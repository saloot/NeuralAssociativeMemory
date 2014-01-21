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


% function clustered_classifier_reverse(N,K,L,no_simulated_instances,err_bits,max_noise_amp,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,try_max,index)

%=============================INITIALIZATION===============================


%----------------------------Input Parameters------------------------------
target_class = 4;                           % The class that input patterns belong to
no_of_classes = 5;                          % The total number of classes there is
alpha0 = 0.95;                              % The optimization update step in the learning algorithm
beta0= 0.75;                                % The sparsity penalty coefficient in the learning algorithm
theta0 = 0.011;                              % The sparsity threshold in the learning algorithm
N = 576;                                    % The size of the network
test_flag = 1;                              % Indictaes if we are using the test dataset or the train dataset
success_count = 0;                          % The number of successful classification
dataset_zero_thr = .05;                     % The threshold below which the entires in the normalize dataset are considered as zero and above which are considered as +1/-1    
gamma_BFO = .95;
gamma_BFS = 0.95;

db_name_in = 'final_database_vectorized';   % The name of the input dataset containing the target class

if (test_flag)        
    db_file_in = (['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Test_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_',num2str(target_class),'_Test_Layer_2.mat']);        
else    
    db_file_in = (['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_',num2str(target_class),'_Train_Layer_2.mat']);    
end
%--------------------------------------------------------------------------

%---------------------------Load the Input Database------------------------
load(db_file_in);
eval(['CIFAR_DB=',db_name_in,';']);
[dataset_size,pattern_length] = size(CIFAR_DB);

a = sqrt(sum(CIFAR_DB'.*CIFAR_DB'));
CIFAR_DB = CIFAR_DB./(a'*ones(1,pattern_length));                           % Normalize the input dataset
CIFAR_DB  = (CIFAR_DB>dataset_zero_thr)-(CIFAR_DB<-dataset_zero_thr);       % Change the dataset to a -1/0/+1 dataset
%--------------------------------------------------------------------------

%-----------------------Load the Learned Matrices--------------------------
for j = 1:no_of_classes        
    %------------------Read the Already Learned Constraints----------------
    fid = fopen(['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(j),'/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Learn_Results/Zero_One/N_',num2str(N),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'r');                
    if (fid > -1)                            
        W = fscanf(fid, '%f',[N,inf]);                                
        W = W';                                
        fclose(fid);                                                                            
        [m,~] = size(W);                       
    else        
        m = 0;
    end    
    %----------------------------------------------------------------------
                
    eval(['W_class',num2str(j),'=W;']);
        
end

%-------------------------Other Initializations----------------------------
% mkdir(['/scratch/amir/Clustered_Neural/Classify_Results'],['N_',num2str(N)]);        % Create a specific folder for the current N and K
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%==========================================================================



%%
%===============================MAIN LOOP==================================
for itr = 1:dataset_size                                    % Simulate over the patterns in the target dataset
        
    
    %-----------------------Generate the Pattern---------------------------
    mu = 1+floor(rand*dataset_size);        
    pattern = CIFAR_DB(mu,:);    
    %----------------------------------------------------------------------                                
                                                                                                                
    max_sum = 0;
    min_norm = 10000;
    
    %---------------------Perform the Classification Step------------------
    success_flag = 0;    
    try_itr = 0;
    
    for j = 1:no_of_classes
        eval(['W_tot = W_class',num2str(j),';']);           % Load the corresponding weight matrix

        %---------------------Iterate Until Convergence----------------            
        [x_out] = recall_step_clustered(W_tot,pattern,N,gamma_BFO,gamma_BFS,1,1,-100,100,2);    
%         [x_out,itr_total] = classifier_faultyrecall_step(W_out,pattern,N,gamma_BFO,gamma_BFS,-100,100,1,.6,const_neur_noise,const_update_threshold)
        %--------------------------------------------------------------
                        
                
        %----------Verify the Success of the Recall Algorithm----------
        if (mean(abs(x_out*W_tot'))<min_norm)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                                       
            pattern_class = j;
            min_norm = mean(abs(x_out*W_tot'));
        end        
        %--------------------------------------------------------------
        
%     while ((try_itr < try_max) && (success_flag == 0))
%         success_flag = 1;
%         try_itr = try_itr+1;
%         for l = 1:L
%             
%             
%             %---------------------Read the Weight Matrix-------------------
%             load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);
%             n = length(index_l);
%             fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt'], 'r');
%             if (fid == -1)
%                 l
%                 error('No weight matrix!');
%             end
%             %['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt']
%             W = fscanf(fid, '%f',[n,inf]);        
%             W = W';      
%             [m,~] = size(W);
%             fclose(fid);
%             
%             if (l == 55)
%                 111;
%             end
%             x_temp = zeros(1,length(index_l));
%             pattern_temp = zeros(1,length(index_l));
%             for iji = 1:length(index_l)
%                 pattern_temp(iji) = pattern(index_l(iji));
%                 x_temp(iji) = x(index_l(iji));
%             end    
%             %--------------------------------------------------------------    
            
            

            
%         end
    end
    
    if (pattern_class==target_class)            
        success_count = success_count+1;
    end
    %----------------------------------------------------------------------                                                                           
 
    
    %--------------------------Display Progress----------------------------
    if (mod(itr, 20) == 0)
        success_count/itr
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

    
