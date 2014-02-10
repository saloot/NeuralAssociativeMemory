%==========================================================================
%********************FUNCTION: spatially_neural****************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N_horiz: The number of columns in the 2D patterns
% N_vert: The number of rows in the 2D patterns
% L_horiz: The number of clusters within a neural plane.
% L_vert: The number of neural planes
% N_const_in: Number of constraints which must be learned during the learning phase
% N_const_per_job: Number of constraints learned by each submitted job to the cluster
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
% p_err_init: Initial 'bit' error probability. It can be a vector.
% max_noise_amp: The maximum amplitude of noise during the recall process
% no_simulated_instances: The number of patterns considered in the recall phase.
% simul_instances_per_job: Number of simulated instances in the recall phase per submitted job to the cluster
% recall_algorithm_option: The recall algorithm identifier (0 for winner-take-all, 1 for the original bit flipping and 2 for the simplified bit flipping
% try_max: The maximum number of iterations in the global recall algorithm
% fixed_index: The 'length' of the rectangular window which is used to fix the boundaries during the recall phase. Put 0 for the uncostrained system
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a spatially-coupled neural associative
% memory and performs the learning and recall phases. More specifically, 
% the function first learns a pre-specified number of constraints. Then the
% function goes directly to the recall phase. 
% Each step is performed by dividing the jobs in parallel and submitting
% them to the cluster. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function spatially_neural(N_horiz,N_vert,L_horiz,L_vert,N_const_in,N_const_per_job,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,p_err_init,max_noise_amp,no_simulated_instances,simul_instances_per_job,recall_algorithm_option,try_max,fixed_index)

%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Spatially_Coupled/Initialization_Files/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_parameters_N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'.mat']);           
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
 

%%
%=============================LEARNING PHASE===============================
for l_v = 1:L_vert                          
    for l_h = 1:L_horiz                       
        
        n = N_loc_horiz*N_loc_vert;
        fid = fopen(['/scratch/amir/Spatially_Coupled/Learn_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(l_h),'_cluster_vert_',num2str(l_v),'.txt'], 'r');        
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
        if (m < N_const_in)
            fid = fopen(['/scratch/amir/Spatially_Coupled/Learn_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(l_h),'_cluster_vert_',num2str(l_v),'.txt'], 'r'); 
                    
            %-------Take Into Account the Constraints Being Learned--------
            if (fid > -1)
                current_count = max(fscanf(fid, '%d'),0);               
                fclose(fid);
            else
                current_count = 0;
            end
            %--------------------------------------------------------------
            
            if (m + current_count < N_const_in)
                for a =1:floor((N_const_in-m-current_count)/N_const_per_job)
           
                    %------------Submit the Job to the Cluster-------------
                    command = ['cd /scratch/amir/Spatially_Coupled/Submit_Cluster;qsub -N "spatially_neur_N_h',num2str(N_horiz),'N_v',num2str(N_vert),'_learn" -v N_hor=',...
                        num2str(N_horiz),',N_ver=',num2str(N_vert),',L_horiz=',num2str(L_horiz),',L_vert=',num2str(L_vert),',const_var=',num2str(N_const_per_job)...
                        ,',cluster_horiz=',num2str(l_h),',cluster_vert=',num2str(l_v),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),...
                        ',theta0=',num2str(theta0),',db_name=''"''"''',db_name,'''"''"'',db_file_name=''"''"''',db_file_name,'''"''"'' spatally_neural_learn.pbs'];        
                    [channel, result]  =  sshfrommatlabissue(channel,command);             
                    %------------------------------------------------------
    
                    %----------Check the success of the submission---------
                    if (isequal(result, {''})) 
                        display('Unsubmitted learning job!');
                    else
                        display('Learning job submitted successfully!');

                        %--Update the Number of Constraints in the Process-
                        current_count = current_count + N_const_per_job;
                        %--------------------------------------------------
                    end
                    %------------------------------------------------------
                end
        
                %----------------Submit the Job to the Cluster-------------
                if (N_const_in-m-current_count-floor((N_const_in-m-current_count)/N_const_per_job)*N_const_per_job>0)                    
                    command = ['cd /scratch/amir/Spatially_Coupled/Submit_Cluster;qsub -N "spatially_neur_N_h',num2str(N_horiz),'N_v',num2str(N_vert),'_learn" -v N_hor=',...
                        num2str(N_horiz),',N_ver=',num2str(N_vert),',L_horiz=',num2str(L_horiz),',L_vert=',num2str(L_vert),',const_var=',num2str(N_const_in-m-current_count-floor((N_const_in-m-current_count)/N_const_per_job)*N_const_per_job)...
                        ,',cluster_horiz=',num2str(l_h),',cluster_vert=',num2str(l_v),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),...
                        ',theta0=',num2str(theta0),',db_name=''"''"''',db_name,'''"''"'',db_file_name=''"''"''',db_file_name,'''"''"'' spatally_neural_learn.pbs'];        
                    [channel, result]  =  sshfrommatlabissue(channel,command);                                         

                    %----------Check the success of the submission---------
                    if (isequal(result, {''}))     
                        display('Unsubmitted learning job!');        
                    else            
                        %--Update the Number of Constraints in the Process-
%                         current_count = current_count + N_const_in-m-current_count-floor((N_const_in-m-current_count)/N_const_per_job)*N_const_per_job;
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
%===================CHECK IF LEARNING PHASE IS FINISHED====================

%---------------Constantly Monitor for the End of Learning Phase-----------
learn_flag = 1;
for l_v = 1:L_vert
        
    %--------------Check if the Learning Phase Is Finished-----------------
    for l_h = 1:L_horiz                               
        
        %------------Read the Current Weight Matrix from the File----------                        
        fid = fopen(['/scratch/amir/Spatially_Coupled/Learn_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(l_h),'_cluster_vert_',num2str(l_v),'.txt'], 'r');
              
        if (fid > -1)    
            W = fscanf(fid, '%f',[n,inf]);        
            W = W';        
            fclose(fid);                                                    
            [m,~] = size(W);        
        else        
            m = 0;       
        end
        %------------------------------------------------------------------        
        
        %-------------------If Learning Is Not Finished--------------------
        if (m < N_const_in-5)                        
            learn_flag = 0;
        end
        %------------------------------------------------------------------        
    end
    %----------------------------------------------------------------------
end
%==========================================================================

%%
%================PERFORM RECALL IF LEARNING PHASE IS FINISHED==============
if (learn_flag == 1)                                  
    
    for noise_itr = 1:length(p_err_init)                
        err_bits = p_err_init(noise_itr);                                       % Determine the number of erroneous bits.        
        
        for net_simul_itr = 1:no_simulated_instances/simul_instances_per_job   % Simulate the error correction procedure for the given ensemble.                                                    
            %-------------Submit Recall Jobs to the Cluster------------                                            
            
            command = ['cd /scratch/amir/Spatially_Coupled/Submit_Cluster;qsub -N "spatially_neur_N_h',num2str(N_horiz),'N_v',...
                num2str(N_vert),'_recall" -v N_hor=',num2str(N_horiz),',N_ver=',num2str(N_vert),',L_hor=',num2str(L_horiz),...
                ',L_ver=',num2str(L_vert),',max_inst=',num2str(simul_instances_per_job),...
                ',err_bit=',num2str(err_bits),',noise_amp=',num2str(max_noise_amp),',alpha0=',num2str(alpha0),...
                ',beta0=',num2str(beta0),',theta0=',num2str(theta0),',gamma_BFO=',num2str(gamma_BFO),',gamma_BFS=',num2str(gamma_BFS),...
                ',recall_algorithm_option=',num2str(recall_algorithm_option),',try_max=',num2str(try_max),',fixed=',num2str(fixed_index),...
                ' spatially_recall_fixed.pbs'];                            
            [channel, result]  =  sshfrommatlabissue(channel,command);     
                                                               
            if (isequal(result, {''}))                             
                error('Unsubmitted recall job!');                                                             
            else                    
                display('Recall job submitted successfully!');                    
            end            
            %----------------------------------------------------------                                        
        end                
    end
    %----------------------------------------------------------------------        
    
end
%--------------------------------------------------------------------------

%==========================================================================

%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
%==========================================================================