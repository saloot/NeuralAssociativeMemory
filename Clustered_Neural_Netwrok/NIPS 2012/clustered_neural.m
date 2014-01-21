%==========================================================================
%********************FUNCTION: clustered_neural****************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% N_const: Number of constraints which must be learned during the learning phase
% N_const_per_job: Number of constraints learned by each submitted job to the cluster
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
% error_bits: Number of initial noisy nodes. It can be a vector
% max_noise_amp: The maximum amplitude of noise during the recall process
% no_simulated_instances: The number of patterns considered in the recall phase.
% simul_instances_per_job: Number of simulated instances in the recall phase per submitted job to the cluster
% recall_algorithm_option: The recall algorithm identifier (0 for winner-take-all, 1 for the original bit flipping and 2 for the simplified bit flipping
% try_max: The maximum number of iterations in outside-cluster recall algorithm
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a clustered neural associative
% memory and performs the learning and recall phases. More specifically, 
% the function first learns a pre-specified number of constraints. Then the
% function goes directly to the recall phase. 
% Each step is performed by dividing the jobs in parallel and submitting
% them to the cluster. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


% function clustered_neural(N,K,L,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,error_bits, max_noise_amp,no_simulated_instances,recall_algorithm_option,try_max)
N = 40;
K = 20;
L = 50;
alpha0 = 0.95;
beta0 = 0.75;
theta0 = 0.05;
gamma_BFO = 0.95;
gamma_BFS = 0.82;
max_noise_amp = 1;
error_bits = [30:10:200];
no_simulated_instances = 2000;
recall_algorithm_option = 2;
try_max = 10;
%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%--------------------Initialize the SSH Connection-------------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');
%--------------------------------------------------------------------------

%==========================================================================
N_const = 20;
N_const_per_job = N_const;
simul_instances_per_job = no_simulated_instances;


%%
%=============================LEARNING PHASE===============================
for index = 1:index_max                     % Perform the learning for different random simulation scenarios
    for cluster = 1:L                       % Execute the learning algorithm for each cluster
        load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster),'.mat']);
        n = length(index_l);                % Load the number of pattern nodes in the corresponding cluster
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'r');        
        
        %-------------Read the Already Learned Constraints-----------------
        if (fid > -1)    
            W = fscanf(fid, '%f',[n,inf]);        
            W = W';        
            fclose(fid);                                                    
            [m,~] = size(W);        
        else        
            m = 0;       
        end
        %------------------------------------------------------------------
        
        %-------------Learn More Constraints if Needed---------------------
        if (m < N_const)
            fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'r'); 
            
            %-------Take Into Account the Constraints Being Learned--------
            if (fid > -1)
                current_count = max(fscanf(fid, '%d'),0);               
                fclose(fid);
            else
                current_count = 0;
            end
            %--------------------------------------------------------------
            
            if (m + current_count < N_const)
                for a =1:floor((N_const-m-current_count)/N_const_per_job)
           
                    %------------Submit the Job to the Cluster-------------
                    command = ['cd /scratch/amir/Clustered_Neural/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(cluster),'_learn" -v N_var=',num2str(N),',K_var=',num2str(K),',L_var=',num2str(L),',const_var=',num2str(N_const_per_job),',train_set_var=',num2str(cluster),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',index=',num2str(index),' clustered_neural_learn.pbs'];        
                    [channel, result]  =  sshfrommatlabissue(channel,command);             
                    %------------------------------------------------------
    
                    %----------Check the success of the submission---------
                    if (isequal(result, {''})) 
                        display('Unsubmitted learning job!');
                    else
                        display('Learning job submitted successfully!');

                        %--Update the Number of Constraints in the Process-
                        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'w'); 
                        fprintf(fid,'%d',current_count+N_const_per_job);
                        fclose(fid);
                        current_count = current_count + N_const_per_job;
                        %--------------------------------------------------
                    end
                    %------------------------------------------------------
                end
        
                %----------------Submit the Job to the Cluster-------------
                if (N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job>0)
                    command = ['cd /scratch/amir/Clustered_Neural/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(cluster),'_learn" -v N_var=',num2str(N),',K_var=',num2str(K),',L_var=',num2str(L),',const_var=',num2str(N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job),',train_set_var=',num2str(cluster),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',index=',num2str(index),' clustered_neural_learn.pbs'];            
                    [channel, result]  =  sshfrommatlabissue(channel,command);                                         

                    %----------Check the success of the submission---------
                    if (isequal(result, {''}))     
                        display('Unsubmitted learning job!');        
                    else            
                        %--Update the Number of Constraints in the Process-
                        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'w'); 
                        fprintf(fid,'%d',current_count+N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job);
                        fclose(fid);
                        current_count = current_count + N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job;
                        %--------------------------------------------------                        
                        display('Learning job submitted successfully!');            
                    end    
                    %------------------------------------------------------
                end
            end
            %--------------------------------------------------------------
        end
        %------------------------------------------------------------------
    

    end
end
%==========================================================================


%%
%=============================RECALL PHASE=================================

%---------------Constantly Monitor for the End of Learning Phase-----------
for index = 1:index_max
    learn_flag = 1;
    
    %--------------Check if the Learning Phase Is Finished-----------------
    for ll = 1:L                               
        %------------Read the Current Weight Matrix from the File----------        
        load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster),'.mat']);
        n = length(index_l);                % Load the number of pattern nodes in the corresponding cluster
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'r');        
        
        %-------------Read the Already Learned Constraints-----------------
        if (fid > -1)    
            W = fscanf(fid, '%f',[n,inf]);        
            W = W';        
            fclose(fid);                                                    
            [m,~] = size(W);        
        else        
            m = 0;       
        end
        %------------------------------------------------------------------        
        
        if (m < N_const*2/3)                        
            learn_flag = 0;
        end
    end
    %----------------------------------------------------------------------
        
    
    %----------------Do the Recall Phase if Learning Is Done---------------
    if (learn_flag == 1)                                  
        for noise_itr = 1:length(error_bits)                
            err_bits = error_bits(noise_itr);                                       % Determine the number of erroneous bits.        
        
            for net_simul_itr = 1:no_simulated_instances/simul_instances_per_job   % Simulate the error correction procedure for the given ensemble.                                                    
                %-------------Submit Recall Jobs to the Cluster------------                                            
                command = ['cd /scratch/amir/Clustered_Neural/Submit_Cluster;qsub -N "adaptive_neur_N_'...
                    ,num2str(N),'K',num2str(K),'C',num2str(l),'_recall" -v N_var=',num2str(N),...
                    ',K_var=',num2str(K),',L_var=',num2str(L),',max_inst=',num2str(simul_instances_per_job),...
                    ',err_bit=',num2str(err_bits),',noise_amp=',num2str(max_noise_amp),',alpha0=',num2str(alpha0)...
                    ,',beta0=',num2str(beta0),',theta0=',num2str(theta0),',gamma_BFO=',num2str(gamma_BFO),...
                    ',gamma_BFS=',num2str(gamma_BFS),',algorithm_option=',num2str(recall_algorithm_option),...
                    ',try_max=',num2str(try_max),',index=',num2str(index),' clustered_neural_recall.pbs'];                            
                [channel, result]  =  sshfrommatlabissue(channel,command);     
                                                               
                if (isequal(result, {''}))                             
                    error('Unsubmitted recall job!');                                                             
                else                    
                    display('Recall job submitted successfully!');                    
                end            
                %----------------------------------------------------------                                    
            end            
        end                
    end
    %----------------------------------------------------------------------        
    
end
%--------------------------------------------------------------------------

%==========================================================================

%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
%==========================================================================