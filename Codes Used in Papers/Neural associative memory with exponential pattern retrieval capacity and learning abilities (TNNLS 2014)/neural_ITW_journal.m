%==========================================================================
%***********************FUNCTION: neural_journal_v1******************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% N_const_per_job: Number of constraints learned by each submitted job to the cluster
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% error_bits: Number of initial noisy nodes. It can be a vector
% max_noise_amp: The maximum amplitude of noise during the recall process
% no_simulated_instances: The number of patterns considered in the recall phase.
% simul_instances_per_job: Number of simulated instances in the recall phase per submitted job to the cluster
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a neural network and performs the
% learning and recall phases of a bipartite neural associative memory. More 
% specifically, the function first learns a pre-specified number of 
% constraints. Then the function goes directly to the recall phase. 
% Each step is performed by dividing the jobs in parallel and submitting
% them to the cluster. 
% The learning algorithm is based on the approach mentioned in our journal
% paper. The recall process considers three methods: The winner-take-all,
% the original bit flipping (mentioned in our ITW 2011 paper) and the
% simplified bit-flipping (proposed in ISIT 2012).
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function neural_ITW_journal(N,K,alpha0,beta0,theta0,N_const_per_job,error_bits,max_noise_amp,gamma_BFO,no_simulated_instances,simul_instances_per_job,gamma_BFS,mca_flag)

%%
%=============================INITIALIZATION===============================

% N = 400;
% K = 200;
% alpha0 = 0.75;
% beta0 = 1;
% theta0 = .031;
% N_const_per_job = 20;
% error_bits = [15:5:50];
% max_noise_amp = 1;
% gamma_BFO = 1;
% gamma_BFS = 0.99;
% no_simulated_instances = 2000;
% simul_instances_per_job = no_simulated_instances;
% mca_flag = 0;

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_parameters_v1_N_',num2str(N),'_K_',num2str(K),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
N_const = N-K;
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%--------------------Initialize the SSH Connection-----------------            
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');
%------------------------------------------------------------------

%==========================================================================
    

%%
%=============================LEARNING PHASE===============================
for index = 1:index_max                                                 % Perform the learning for different random simulation scenarios
    
    %---------------Read the Already Learned Constraints-------------------
    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(index),'.txt'], 'r');        
    if (fid > -1)    
        W = fscanf(fid, '%f',[N,N_const]);        
        W = W';        
        fclose(fid);                            
        [m,~] = size(W);        
    else        
        m = 0;       
    end
    %----------------------------------------------------------------------
    
    %---------------Learn More Constraints if Needed-----------------------
    if (m < N_const)
        fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(index),'.txt'], 'r'); 
        %---------Take Into Account the Constraints Being Learned----------
        if (fid > -1)
            current_count = max(fscanf(fid, '%d'),0);            
            fclose(fid);
        else
            current_count = 0;
        end
        %------------------------------------------------------------------
        
        
        if (m + current_count < N_const)
        
        
        %------------------Keep Track of SSH Created Channels--------------
%         load('/scratch/amir/ssh_connections.mat');
%         connection_list = [connection_list;channel];
%         save('/scratch/amir/ssh_connections.mat','connection_list')
        %------------------------------------------------------------------            
        
            for a =1:floor((N_const-m-current_count)/N_const_per_job)
           
                %--------------Submit the Job to the Cluster---------------
                command = ['cd cluster/ITW_Journal/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(index),'" -v N_var=',num2str(N),',K_var=',num2str(K),',const_var=',num2str(N_const_per_job),',train_set_var=',num2str(index),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',mca_flag=',num2str(mca_flag),' neural_jou_cluster_learn.pbs'];        
                [channel, result]  =  sshfrommatlabissue(channel,command);             
                %----------------------------------------------------------
    
                %-----------Check the success of the submission------------
                if (isequal(result, {''})) 
                    display('Unsubmitted learning job!');
                else
%                     %----Update the Number of Constraints in the Process---
%                     fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(index),'.txt'], 'w'); 
%                     fprintf(fid,'%d',current_count+N_const_per_job);
%                     fclose(fid);
                    current_count = current_count + N_const_per_job;
%                     %------------------------------------------------------
                    display('Learning job submitted successfully!');
                end
                %----------------------------------------------------------
            end
        
            %----------------Submit the Job to the Cluster-----------------
            if (N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job>0)
                command = ['cd cluster/ITW_Journal/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(index),'" -v N_var=',num2str(N),',K_var=',num2str(K),',const_var=',num2str(N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job),',train_set_var=',num2str(index),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',mca_flag=',num2str(mca_flag),' neural_jou_cluster_learn.pbs'];            
                [channel, result]  =  sshfrommatlabissue(channel,command);                 
                    
                %------------Check the success of the submission-----------
                if (isequal(result, {''}))     
                    display('Unsubmitted learning job!');        
                else        
                    %----Update the Number of Constraints in the Process---
%                     fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(index),'.txt'], 'w'); 
%                     fprintf(fid,'%d',current_count+N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job);
%                     fclose(fid);
                    current_count = current_count + N_const-m-current_count-floor((N_const-m-current_count)/N_const_per_job)*N_const_per_job;
                    %------------------------------------------------------
                    display('Learning job submitted successfully!');            
                end    
                %----------------------------------------------------------
            end
            %--------------------------------------------------------------
        end        
    end
    %----------------------------------------------------------------------  
end
%==========================================================================


%%
%=============================RECALL PHASE=================================

%---------Perform the Recall Phase for the Finished Learning Jobs----------
for l =1:index_max                                % For all simulation scenarios do...
    
    %-------------Read the Current Weight Matrix from the File-------------        
    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(l),'.txt'], 'r');        
    if (fid > -1)
        W = fscanf(fid, '%f',[N,N_const]);
        W = W';
        fclose(fid);                        
        [m,~] = size(W);        
    else        
        m = 0;
    end
    %----------------------------------------------------------------------
    
    %----------------Do the Recall Phase if Learning Is Done---------------
    if (m >= N_const)                                                                
        for noise_itr = 1:length(error_bits)                        
            err_bits = error_bits(noise_itr);       % Determine the number of erroneous bits.                    
            for net_simul_itr = 1:ceil(no_simulated_instances/simul_instances_per_job)   % Simulate the error correction procedure for the given ensemble.
                                            
                %-------------Submit Recall Jobs to the Cluster------------                                
                command = ['cd cluster/ITW_Journal/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(l),'" -v N_var=',num2str(N),',K_var=',num2str(K),',max_inst=',num2str(simul_instances_per_job),',err_bit=',num2str(err_bits),',noise_amp=',num2str(max_noise_amp),',ind=',num2str(l),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',gamma_BFO=',num2str(gamma_BFO),',gamma_BFS=',num2str(gamma_BFS),' neural_jou_cluster_recall.pbs'];                            
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
    %------------------------------------------------------------------
    
end   
%--------------------------------------------------------------------------
%==========================================================================


%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
%==========================================================================