%==========================================================================
%*********************FUNCTION: faulty_neural******************************
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
% pattern_neur_noise_vector: The vector of maximum amount of noise that pattern neurons will have
% const_neur_noise_vector: The vector of maximum amount of noise that constraint neurons will have
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a clustered neural associative
% memory and performs the learning and recall phases. More specifically, 
% the function first learns a pre-specified number of constraints. Then the
% function goes directly to the recall phase. 
% The recall phase is done using faulty neurons, i.e. neurons that have
% some inherent noise in making decisions. 
% Each step is performed by dividing the jobs in parallel and submitting
% them to the cluster. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


% function faulty_neural(N,K,L,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,error_bits, max_noise_amp,no_simulated_instances,simul_instances_per_job,recall_algorithm_option,try_max,pattern_neur_noise_vector,const_neur_noise_vector)

%%
%=============================INITIALIZATION===============================

%--------------------------Set Default Values------------------------------

max_noise_amp = 3;
no_simulated_instances = 2000;
simul_instances_per_job = 2000;
try_max = 40;


if ~exist('error_bits','var')
    error_bits = [1:8:41];
end

if ~exist('N_in','var')
    N = 40;
end

if ~exist('K_in','var')
    K = 20;
end

if ~exist('L_in','var')
    L = 50;
end

if ~exist('alpha0','var')
    alpha0 = 0.95;
end

if ~exist('beta0','var')
    beta0 = 0.75;
end

if ~exist('theta0','var')
    theta0 = 0.05;
end

if ~exist('gamma_BFO','var')
    gamma_BFO = 0.95;
end

if ~exist('gamma_BFS','var')
    gamma_BFS = 0.95;
end

if ~exist('recall_algorithm_option','var')
    recall_algorithm_option = 1;
end

if ~exist('pattern_neur_noise_vector','var')
    pattern_neur_noise_vector = [0,.1,.3,.4,.6];
end

if ~exist('const_neur_noise_vector','var')
    const_neur_noise_vector = [0,.05,.2,.3];
end
%--------------------------------------------------------------------------

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%==========================================================================
N_const = min(40,N-K);
N_const_per_job = N_const;



%%
%=============================LEARNING PHASE===============================
for index = 1:index_max                     % Perform the learning for different random simulation scenarios
    for cluster = 1:L                       % Execute the learning algorithm for each cluster
        load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster),'.mat']);
        n = length(index_l);                % Load the number of pattern nodes in the corresponding cluster
        fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'r');        
        
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
            fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'r'); 
            
            %-------Take Into Account the Constraints Being Learned--------
            if (fid > -1)
                current_count = max(fscanf(fid, '%d'),0);               
                fclose(fid);
            else
                current_count = 0;
            end
            %--------------------------------------------------------------
            
            if (m + current_count < N_const)
                
                %--------------------Initialize the SSH Connection-----------------
                channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');
                %------------------------------------------------------------------
        
                %------------------Keep Track of SSH Created Channels--------------
%                 load('/scratch/amir/ssh_connections.mat');
%                 connection_list = [connection_list;channel];
%                 save('/scratch/amir/ssh_connections.mat','connection_list')
                %------------------------------------------------------------------
                
                for a =1:floor((N_const-m-current_count)/N_const_per_job)
           
                    %------------Submit the Job to the Cluster-------------
                    command = ['cd /scratch/amir/Fault_Tolerant_Decoding/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(cluster),'_learn" -v N_var=',num2str(N),',K_var=',num2str(K),',L_var=',num2str(L),',const_var=',num2str(N_const_per_job),',cluster_index=',num2str(cluster),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',index=',num2str(index),' faulty_neural_learn.pbs'];        
                    [channel, result]  =  sshfrommatlabissue(channel,command);             
                    %------------------------------------------------------
    
                    %----------Check the success of the submission---------
                    if (isequal(result, {''})) 
                        display('Unsubmitted learning job!');
                    else
                        display('Learning job submitted successfully!');

                        %--Update the Number of Constraints in the Process-
%                         fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'w'); 
%                         fprintf(fid,'%d',current_count+N_const_per_job);
%                         fclose(fid);
                        current_count = current_count + N_const_per_job;
                        %--------------------------------------------------
                    end
                    %------------------------------------------------------
                end
        
                %----------------Submit the Job to the Cluster-------------
                if (N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job>0)
                    command = ['cd /scratch/amir/Fault_Tolerant_Decoding/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(cluster),'_learn" -v N_var=',num2str(N),',K_var=',num2str(K),',L_var=',num2str(L),',const_var=',num2str(N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job),',cluster_index=',num2str(cluster),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',index=',num2str(index),' faulty_neural_learn.pbs'];            
                    [channel, result]  =  sshfrommatlabissue(channel,command);                                         

                    %----------Check the success of the submission---------
                    if (isequal(result, {''}))     
                        display('Unsubmitted learning job!');        
                    else            
                        %--Update the Number of Constraints in the Process-
%                         fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'w'); 
%                         fprintf(fid,'%d',current_count+N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job);
%                         fclose(fid);
%                         current_count = current_count + N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job;
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
%--------------------Initialize the SSH Connection-----------------                   
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');          
%------------------------------------------------------------------
                            
%------------------Keep Track of SSH Created Channels--------------            
% load('/scratch/amir/ssh_connections.mat');                    
% connection_list = [connection_list;channel];                    
% save('/scratch/amir/ssh_connections.mat','connection_list')                        
%------------------------------------------------------------------
counter = 0;
for index = 1:index_max
    learn_flag = 1;
    
    %--------------Check if the Learning Phase Is Finished-----------------
    for ll = 1:L                               
        %------------Read the Current Weight Matrix from the File----------                
        load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster),'.mat']);
        n = length(index_l);                % Load the number of pattern nodes in the corresponding cluster
        fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'r');        
        
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
        
        if (m < N_const/2)                        
            learn_flag = 0;
        end
        counter = counter + 1;
    end
    %----------------------------------------------------------------------
        
    
    %----------------Do the Recall Phase if Learning Is Done---------------    
   
    if (learn_flag == 1)                                  
        for noise_itr = 1:length(error_bits)                
            err_bits = error_bits(noise_itr);                                       % Determine the number of erroneous bits.        
            for jj = 1:length(pattern_neur_noise_vector);
                pattern_neur_noise = pattern_neur_noise_vector(jj);
                for kk = 1:length(const_neur_noise_vector);
                    const_neur_noise = const_neur_noise_vector(kk);
                    for net_simul_itr = 1:no_simulated_instances/simul_instances_per_job   % Simulate the error correction procedure for the given ensemble.                                                    
                        %-------------Submit Recall Jobs to the Cluster------------                                            
                        command = ['cd /scratch/amir/Fault_Tolerant_Decoding/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(l),'_recall" -v N_var=',num2str(N),',K_var=',num2str(K),',L_var=',num2str(L),',max_inst=',num2str(simul_instances_per_job),',err_bit=',num2str(err_bits),',noise_amp=',num2str(max_noise_amp),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',gamma_BFO=',num2str(gamma_BFO),',gamma_BFS=',num2str(gamma_BFS),',algorithm_option=',num2str(recall_algorithm_option),',try_max=',num2str(try_max),',index=',num2str(index),',pattern_neur_noise=',num2str(pattern_neur_noise),',const_neur_noise=',num2str(const_neur_noise),' faulty_neural_recall.pbs'];                            
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
        end                
    end
    %----------------------------------------------------------------------        
    counter
end
%--------------------------------------------------------------------------

%==========================================================================

%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
%==========================================================================