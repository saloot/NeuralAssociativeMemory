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


% function clustered_neural_recall(cluster_size,const_learn,no_clusters,db_file_in_orig,db_name_in,...
%     db_file_projected,db_name_projected,no_simulated_instances,err_bits,max_noise_amp,alpha0,beta0,theta0,...
%     gamma_BFO,gamma_BFS,recall_algorithm_option,try_max,Q,learn_itr_max,simulation_set)



%=============================INITIALIZATION===============================

%---------------------------Load the Training Dataset----------------------
load(db_file_in_orig);
eval(['dataset_learn_orig = ',db_name_in,';']);   
[dataset_size,N] = size(dataset_learn);
% load(db_file_projected);
% eval(['dataset_recall = ',db_name_projected,';']);   
y_min = min(min(dataset_recall));
y_max = max(max(dataset_recall));
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
destination_folder_learn = [db_file_in(1:i-1),'//Learn_Results/Zero_One/Clustered_Version'];
if (~exist(destination_folder_recall,'dir'))
    mkdir(destination_folder_recall);
end
%--------------------------------------------------------------------------

%--------------------------Quantize if Necessary---------------------------
if (Q > 0)    
    dataset_learn = round(Q*tanh(sigma1*dataset_learn_orig-mean1));
    dataset_learn = dataset_learn;    
else
    dataset_learn = dataset_learn_orig;   
end
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
dataset_learn = [ones(dataset_size,1),dataset_learn];
%--------------------------------------------------------------------------


%==========================================================================



%%
%===============================MAIN LOOP==================================
for net_simul_itr = 1:no_simulated_instances                % Simulate the error correction procedure for the given ensemble.                                
        
    %-----------------------------Generate Noise---------------------------
    nois = zeros(1,N);                                  % This is the noise added to the whole pattern of length N*L_in        
    pp = 1+floor((N-1)*rand(1,err_bits));                
    for h = 1:err_bits        
        nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
    end        
    %----------------------------------------------------------------------                            
                
    %-----------------------Generate the Pattern---------------------------
    mu = 1+floor((dataset_size-1)*rand);            % Pick a pattern index at random                   
    pattern = dataset_learn(mu,:);    
    projected_pattern = dataset_recall(mu,:);
    learning_error(mu) = sum(abs(sign(pattern(2:end)-dataset_recall(mu,2:end))));
    
    %----------------------------------------------------------------------                                
                                                                                                                
    %------------Initialize the Network with a Noisy SubPattern------------
    x = pattern(2:end)+nois(2:end);                                       % Initialize the network with a nosiy version of the subpattern               
    %----------------------------------------------------------------------
            
        
    %------------------------Perform the Recall Step-----------------------
    success_flag = 0;    
    try_itr = 0;
    
    W_global_tot = [];
    for l = 1:no_clusters                        
        load(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(l)]);            
        save_flag = save_place(end,2);                        
        if (save_flag == 2)                           
            destination_folder_save = [destination_folder_learn,'/Partial_Convergence'];               
        else            
            destination_folder_save = destination_folder_learn;
        end        
        %-----------------------Store the Connectivity Matrix----------------------
        file_name = [destination_folder_save,'/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...            
        num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];            
        fid = fopen(file_name,'r');        
        if (fid > -1)               
            load(file_name)
            W_global_tot = [W_global_tot;W_global];
            
        end
    end
    
    W_tot = soft_threshold_matrix(W_global_tot,.03);                
    W_tot = W_tot./(sqrt(sum(W_tot'.*W_tot'))'*ones(1,size(W_tot,2)));                
    W_global_deg_dist = sum(abs(sign(W_tot)));
        
    %---------------------Remove Zero Columns--------------------------
    [val,ind] = sort(W_global_deg_dist);
    for i = 1:length(val)
        if (val(i) > 6)
            break;
        else
            W_tot(:,ind(i)) = W_tot(:,ind(i)) + .05*random_vector(size(W_tot,1),6-val(i))';
        end
    end
%         W_global = W_global(:,ind(i:end));        
    W_tot = W_tot./(sqrt(sum(W_tot'.*W_tot'))'*ones(1,size(W_tot,2)));        

    
    
    while ((try_itr < try_max) && (success_flag == 0))
        success_flag = 1;
        try_itr = try_itr+1;
        consts_so_far = 1;
        for l = 1:no_clusters
            
            
            load(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(l)]);
            save_flag = save_place(end,2);            
            if (save_flag == 2)           
                destination_folder_save = [destination_folder_learn,'/Partial_Convergence'];   
            else    
                destination_folder_save = destination_folder_learn;
            end
            %-----------------------Store the Connectivity Matrix----------------------
            file_name = [destination_folder_save,'/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
            num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];
            fid = fopen(file_name,'r');
            if (fid > -1)
                load(file_name)
            else
                error('Invalid cluster file')                            
            end
            
            consts = size(W_total,1);
            
            
            W = [W_total(:,1),W_tot(consts_so_far:consts_so_far+consts-1,index_pattern_neurons+1)];
            consts_so_far = consts_so_far + consts;
            x_temp = [8,x(index_pattern_neurons)];                                    
            pattern_temp = [8,pattern(1+index_pattern_neurons)]; 
            projected_pattern_temp =[8,projected_pattern(1+index_pattern_neurons)];
            %--------------------------------------------------------------    
            
            %---------------------Iterate Until Convergence----------------            
            [x_out] = recall_step_clustered(W,x_temp,gamma_BFO,gamma_BFS,y_min,y_max,recall_algorithm_option,max_y_threshold);    
            %--------------------------------------------------------------
            
            %----------Verify the Success of the Recall Algorithm----------
            if (norm(x_out-pattern_temp)<.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                       
                for iji = 1:length(index_l)
                    x(index_l(iji)) = x_out(iji); 
                end
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
    if (sum(recall_error)>0)
        error_count_target = error_count_target+1;
        bit_error_count_target = bit_error_count_target+recall_error;
    end
    
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
    fid = fopen(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(no_clusters),'/clustered_neural_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
    fprintf(fid,'\n');
    fclose(fid);
elseif (recall_algorithm_option == 1)
    fid = fopen(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(no_clusters),'/clustered_neural_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
    fprintf(fid,'\n');
    fclose(fid);
elseif (recall_algorithm_option == 2)
    fid = fopen(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(no_clusters),'/clustered_neural_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
    fprintf(fid,'\n');
    fclose(fid);
else
    error('Unknown recall algorithm');
end
%----------------------------------------------------------------------

    
