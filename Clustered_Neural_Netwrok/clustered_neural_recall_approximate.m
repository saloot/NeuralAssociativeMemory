%==========================================================================
%***************FUNCTION: clustered_neural_recall**************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% no_clusters: The number of clusters if we have multiple levels (no_clusters = 1 for single level)
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
% index: The index of the simulation setup among various random scenarios
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


% function clustered_neural_recall_approximate(cluster_size,const_learn,no_clusters,db_file_in_orig,db_name_in,...
%     db_file_projected,db_name_projected,no_simulated_instances,err_bits,max_noise_amp,alpha0,beta0,theta0,...
%     gamma_BFO,gamma_BFS,recall_algorithm_option,try_max,Q,learn_itr_max,simulation_set)



%=============================INITIALIZATION===============================

%---------------------------Load the Training Dataset----------------------
db_file_in = ['/scratch/amir/Databases/Caltech-101/Caltech101_Silhouettes/caltech101_silhouettes_28.mat'];
db_name_in = ['X'];
load(db_file_in);
eval(['dataset_recall = ',db_name_in,';']);   
clear(db_name_in)
[dataset_size,N] = size(dataset_recall);
%--------------------------------------------------------------------------

%------------------------Compute the Projected Dataset---------------------
% if (no_of_PC > 0)
%     [COEFF,SCORE,latent] = princomp(dataset_recall);
%     dataset_projected = SCORE(:,1:no_of_PC)*COEFF(:,1:no_of_PC)';
%     dataset_recall_orig = dataset_recall;
%     dataset_recall = dataset_projected;
% end
y_min = min(min(dataset_recall));
y_max = max(max(dataset_recall));
no_simulated_instances = 1000;


Q = 1;
alpha0 = 0.95;
beta0 = 1;

learn_itr_max = 250;
%--------------------------------------------------------------------------

%--------------------------Quantize if Necessary---------------------------
if (Q > 0)    
%     dataset_recall = round(Q*tanh(dataset_recall))/Q;
end
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
thr1 = max(max(abs(dataset_recall)));
dataset_recall = [thr1*ones(dataset_size,1),dataset_recall];
%--------------------------------------------------------------------------


%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 

error_count = 0;
bit_error_count = 0;
error_count_target = 0;
bit_error_count_target = 0;
learning_error = zeros(1,dataset_size);
recall_error = zeros(1,dataset_size);
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


destination_folder_recall = [db_file_in(1:i-1),'/Recall_Results/Clustered_Version'];
destination_folder_learn = [db_file_in(1:i-1),'/Learn_Results/Clustered_Version'];
if (~exist(destination_folder_recall,'dir'))
    mkdir(destination_folder_recall);
end
%--------------------------------------------------------------------------

%==========================================================================
theta0 = 0.01;
simulation_set = 1;
no_clusters = 140;
try_max = 80;
cluster_size = 420;
no_clusters = 40;
const_learn = min(cluster_size/2,30);
const_learn = 2* cluster_size;
const_learn = 160;
max_noise_amp = 1;
W_total_tot = [];
W_global_tot = [];
no_const_tot = [];
cluster_range = [];
index_pattern_tot = [];
varphi = .75;
psi = 0.0025;
err_bits =20;
for l = 1:no_clusters                                
        
        %-----------------------Store the Connectivity Matrix----------------------
        file_name = [destination_folder_learn,'/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...            
            num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];            
        fid = fopen(file_name,'r');        
        
        if (fid == -1)
            fclose(fid);
            file_name = [destination_folder_learn,'/Partial_Convergence/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...            
                num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];                    
            fid = fopen(file_name,'r');        
        end
        if (fid > -1)               
            fclose(fid);
            load(file_name)
            W_global_tot = [W_global_tot;W_global];
            W_total_tot = [W_total_tot;W_total];
            [mi,~] = size(W_total);
            no_const_tot = [no_const_tot,mi];
            cluster_range = [cluster_range,l];
            index_pattern_tot = [index_pattern_tot;index_pattern_neurons];
        end
        
end
%%
% W_total_tot = soft_threshold_matrix(W_total_tot,2*psi+.001);

err_bit_effective = [];
BER_out = [];
SNR_in = [];
SNR_out = [];
no_converged_pats = length(index_converged_patterns);

%===============================MAIN LOOP==================================
for net_simul_itr = 1:no_simulated_instances                % Simulate the error correction procedure for the given ensemble.                                
        
    %-----------------------------Generate Noise---------------------------
    nois = zeros(1,N);                                  % This is the noise added to the whole pattern of length N*L_in        
    pp = 1+floor((N-1)*rand(1,err_bits));                
    for h = 1:err_bits        
        nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
    end  
    nois = [0,nois];
    %----------------------------------------------------------------------                            
                
    %-----------------------Generate the Pattern---------------------------
    
    mu = 1+floor((no_converged_pats-1)*rand);            % Pick a pattern index at random                   
    mu = index_converged_patterns(mu);
    pattern = dataset_recall(mu,:);    
%     projected_pattern = dataset_recall(mu,:);
%     learning_error(mu) = norm(abs((pattern-projected_pattern(2:end))));    
    %----------------------------------------------------------------------                                
                                                                                                                
    %------------Initialize the Network with a Noisy SubPattern------------
    x = pattern+nois;                                       % Initialize the network with a nosiy version of the subpattern               
    x = max(min(x,y_max),y_min);
    nois_effective = x - pattern;
    SNR_in = [SNR_in, norm(x)/norm(nois_effective)];
%     sum(abs(sign(x - pattern)))
    err_bit_effective = [err_bit_effective,sum(abs(sign(x - pattern)))];
    %----------------------------------------------------------------------
            
        
    %------------------------Perform the Recall Step-----------------------
    success_flag = 0;    
    try_itr = 0;

    
    while ((try_itr < try_max) && (success_flag == 0))
        success_flag = 1;
        try_itr = try_itr+1;
        consts_so_far = 0;
        counter = 0;
        for l = 1:length(cluster_range)
            consts = no_const_tot(l);
            
            W_total = W_total_tot(1+consts_so_far:consts_so_far + consts,:);
            consts_so_far = consts_so_far + consts;
            
            if (consts < 20)
                continue;
            end
            counter = counter + 1;
            index_pattern_neurons = index_pattern_tot(l,:);
            x_temp = [thr1,x(1+index_pattern_neurons)];
            pattern_temp = [thr1,pattern(1+index_pattern_neurons)];             
            nois_temp = x_temp - pattern_temp;
            if (norm(nois_temp))
                111;
            end
            %--------------------------------------------------------------    
            
            %---------------------Iterate Until Convergence----------------            
            [x_out] =recall_step_clustered_v2(W_total,x_temp,varphi,psi,y_min,y_max,Q);    
            %--------------------------------------------------------------
            
            
            
            %----------Verify the Success of the Recall Algorithm----------
            fraction_satisfied_consts = sum(abs(W_total*x_out')< psi)/consts;
            if (fraction_satisfied_consts > varphi)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                       
                if (norm(x_temp-x_out))
                111;
            end 
                x(1+index_pattern_neurons) = x_out(2:end);                 
            else            
                success_flag = 0;                
            end
            %--------------------------------------------------------------
            
        end
    end
    %----------------------------------------------------------------------
                                                          
                
    %--------------------Calculate Recall Error Count----------------------
    x_target = dataset_recall(mu,:);    
    recall_error(mu) = sum(abs(sign(x_target(2:end)-x(2:end))));
    
    if (sum(recall_error(mu))>0)
        error_count_target = error_count_target+1;
        bit_error_count_target = bit_error_count_target+recall_error(mu);
        SNR_out = norm(x_target)/norm(x_target -x);
        BER_out = [BER_out,sum(recall_error(mu))];
    else
        SNR_out = [SNR_out,1e8+ 1000];
        BER_out = [BER_out,0];
    end

    %----------------------------------------------------------------------
    
    %--------------------------Display Progress----------------------------
    if (mod(net_simul_itr, 100) == 0)
%         error_count_target/net_simul_itr
%         BER_out
        sum((BER_out./err_bit_effective(1:end))<1)
                
        sum((BER_out./err_bit_effective(1:end))>1)
    end
     %----------------------------------------------------------------------

    
end    
        
            
    
%------------------Transform Error Count to Error Rate---------------------
PER = error_count_target/net_simul_itr;
BER = sum(bit_error_count_target)/net_simul_itr/N;
%--------------------------------------------------------------------------


    
%----------------Store the Bit and Pattern Error Rates-----------------
if (recall_algorithm_option == 0)
    fid = fopen(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(no_clusters),'/clustered_neural_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
    fprintf(fid,'\n');
    fclose(fid);
else
    error('Unknown recall algorithm');
end
%----------------------------------------------------------------------

    
