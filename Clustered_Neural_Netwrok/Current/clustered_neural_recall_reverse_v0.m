%==========================================================================
%***************FUNCTION: clustered_neural_recall**************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% no_simulated_instances: The number of patterns considered in the recall phase.
% err_bits: Number of initial noisy nodes.
% max_noise_amp: The maximum amplitude of noise during the recall process
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
% recall_algorithm_option: The recall algorithm identifier (0 for winner-take-all, 1 for the original bit flipping and 2 for the simplified bit flipping
% try_max: The maximum number of iterations in outside-cluster recall algorithm
% db_file_in: The address of the input database
% db_name_in: The name of the input database
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a clustered neural associative
% memory and performs the recall phase. More specifically, 
% the function first reads the weight matrix found in the learning phase. 
% Then the function performs the recall phase introduced in our NIPS 2012
% paper, i.e. perform the recall algorithm within each cluster one after
% another. If the recall process was successful, the state of neurons is
% maintained, and reverted back to their original version otherwise. This
% process is repeated try_max times, after which a recall error is declared
% if the output of the algorithm is not equal to the original pattern.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


% function clustered_neural_recall_reverse(N,L,no_simulated_instances,err_bits,max_noise_amp,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,try_max,db_file_in,db_name_in)

%=============================INITIALIZATION===============================

N = 576;
alpha0=0.95;
beta0 = 0.75;
theta0= 0.008;
recall_algorithm_option = 1;
try_max = 20;
dataset_zero_thr = 0.05;
no_simulated_instances = 1000;
err_bits = 4;
max_noise_amp = 1;
gamma_BFO = .9;
gamma_BFS = 0.92;
no_of_clusters = 250;                                        % The total number of clusters we would like to have
db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_1/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_1_Train_Layer_2.mat';
db_name_in = 'final_database_vectorized';

%---------------------------Load the Database------------------------------
load(db_file_in);
eval(['dataset_recall = ',db_name_in,';']);   
[dataset_size,pattern_length] = size(dataset_recall);

dataset_recall = dataset_recall .*(abs(dataset_recall)>dataset_zero_thr);
dataset_recall = ones(dataset_size,pattern_length).*((dataset_recall>dataset_zero_thr)-(dataset_recall<-dataset_zero_thr));
dataset_recall = abs(dataset_recall);
% a = sqrt(sum(dataset_recall'.*dataset_recall'));
% dataset_recall = dataset_recall./(a'*ones(1,1024));
% dataset_recall = dataset_recall.*(abs(dataset_recall)>.03);
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
dataset_recall = [ones(dataset_size,1),dataset_recall];
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 

error_count = 0;
bit_error_count = 0;

y_min = -10;
y_max = 100;
fixed_set = ones(1,N+1);
%--------------------------------------------------------------------------

%-------------------------Load the Weight Matrix---------------------------
% fid = fopen(weight_matrix_file, 'r');               
% if (fid > -1)                            
%     W = fscanf(fid, '%f',[N,inf]);                                
%     W = W';                                
%     fclose(fid);   
    W = W_global;
    [no_of_constraints,~] = size(W);        
    threshold_zero = .005;
for i = 1:no_of_constraints
    W(i,:) = soft_threshold(W(i,:),threshold_zero);
end
% else    
%     error('Invalid input matrix');
% end

% for i = 1:no_of_constraints
%     W(i,:) = soft_threshold(W(i,:),.005)';
% end
%--------------------------------------------------------------------------

%--------------------Create the Sub-folder If Necessary--------------------
slash_flag = 0;
for i = length(db_file_in):-1:1
    if (strcmp(db_file_in(i),'/'))
        if (slash_flag == 0)
            break;
        else
            slash_flag = slash_flag+1;
        end
    end
end


destination_folder = [db_file_in(1:i-1),'/Recall_Results/Zero_One/N_',num2str(N)];
destination_folder_indices = [db_file_in(1:i-1),'/Learn_Results/Zero_One/N_',num2str(N),'/Cluster_Indices'];
if (~exist(destination_folder,'dir'))
    mkdir(destination_folder);
end
%--------------------------------------------------------------------------

%==========================================================================



%%
%===============================MAIN LOOP==================================
for net_simul_itr = 1:no_simulated_instances                % Simulate the error correction procedure for the given ensemble.                                
        
    %-----------------------------Generate Noise---------------------------
    nois = zeros(1,N+1);                                  % This is the noise added to the whole pattern of length N*L_in        
    pp = 1+floor((N-1)*rand(1,err_bits));                
    for h = 1:err_bits        
        nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
    end        
    nois(1) = 0;
    %----------------------------------------------------------------------                            
                
    %-----------------------Generate the Pattern---------------------------
    mu = 1+floor((dataset_size-1)*rand);            % Pick a pattern index at random                   
    pattern = dataset_recall(mu,:);
    %----------------------------------------------------------------------                                
                                                                                                                
    %------------Initialize the Network with a Noisy SubPattern------------
    x = pattern+nois;                                       % Initialize the network with a nosiy version of the subpattern               
    x = max(x,0);
    x = min(x,1);
    nois = x-pattern;
    x = [x];
    pattern = [pattern];
    sum(abs(nois))
    %----------------------------------------------------------------------
            
        
    %------------------------Perform the Recall Step-----------------------
    success_flag = 0;    
    try_itr = 0;
    
    while ((try_itr < try_max) && (success_flag == 0))
        success_flag = 1;
        try_itr = try_itr+1;
        deg = sum(abs(sign(W)));
        [~,schedule] = sort(deg,'descend'); %randperm(N)
        schedule=randperm(N+1);
        for l = 1:N+1
%             load([destination_folder_indices,'/cluster_indices_N_',num2str(N),'_L_',num2str(no_of_clusters)...
%             ,'_cluster_index_',num2str(l),'.mat'])
            cc = schedule(l);
            if cc == 1
                continue
            end
            fixed_set = zeros(1,N+1);            
            fixed_set(cc) = 1;
            
            noise_effective = (x-pattern).*sign(sum(abs(sign(W_cluster))));
            if (sum(abs(noise_effective))== 0)
                111;
            end
            ss = sign(sum(abs(sign(W))));
            x_temp = zeros(1,sum(sign(sum(abs(sign(W_cluster))))));
            pattern_temp = zeros(1,sum(sign(sum(abs(sign(W_cluster))))));
            ikk = 0;
            for iji = 1:N+1                
                if (ss(iji)>0)
                    ikk = ikk + 1;
                    pattern_temp(ikk) = pattern(iji);
                    x_temp(ikk) = x(iji);
                end
            end    
            %--------------------------------------------------------------

            
            %---------------------Iterate Until Convergence----------------            
%             [x_out] = recall_step_clustered_reverse(W_cluster,x,N,gamma_BFO,gamma_BFS,0,1,recall_algorithm_option,fixed_set);    
            [x_out] = recall_step_clustered_scheduled(W,x,N,gamma_BFO,gamma_BFS,0,1,1,fixed_set);    
            
            ikk = 0;
            x_out_temp = zeros(1,sum(sign(sum(abs(sign(W_cluster))))));
            for iji = 1:N+1                
                if (ss(iji)>0)
                    ikk = ikk + 1;
                    x_out_temp(ikk) = x_out(iji);
                end
            end  
            %--------------------------------------------------------------
            
            %----------Verify the Success of the Recall Algorithm----------
%             if (norm(x_out_temp - pattern_temp)<.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                                       
                
             if (norm(W*x_out') -norm(W*x')<-5e-2)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                                       
                x(cc) =x_out(cc);            
                mean(W*x')
                sum(abs(x-pattern)>.0001)
            end
            
            if (norm(x-pattern) < .01)
                break;
            else                    
                success_flag = 0; 
            end
            %--------------------------------------------------------------
            
        end
    end
    %----------------------------------------------------------------------
                                                          
                
    %--------------------Calculate Recall Error Count----------------------
    if (norm(x-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.
        error_count = error_count+1;
        bit_error_count = bit_error_count+sum(abs(x-pattern)>1e-1);
%         sum(abs(sign(x-pattern)));        
    else
        111;
    end
    %----------------------------------------------------------------------
    
    %--------------------------Display Progress----------------------------
    if (mod(net_simul_itr, 2) == 0)
        error_count/net_simul_itr
    end
     %----------------------------------------------------------------------

    
end    
        
            
    
%------------------Transform Error Count to Error Rate---------------------
PER = error_count/net_simul_itr;
BER = bit_error_count/net_simul_itr/N;
%--------------------------------------------------------------------------


    
%----------------Store the Bit and Pattern Error Rates-----------------
if (recall_algorithm_option == 0)
    fid = fopen(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
    fprintf(fid,'\n');
    fclose(fid);
elseif (recall_algorithm_option == 1)
    fid = fopen(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
    fprintf(fid,'\n');
    fclose(fid);
elseif (recall_algorithm_option == 2)
    fid = fopen(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
    fprintf(fid,'\n');
    fclose(fid);
else
    error('Unknown recall algorithm');
end
%----------------------------------------------------------------------

    
