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


% function clustered_neural(N,K,L,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,error_bits, max_noise_amp,no_simulated_instances,simul_instances_per_job,recall_algorithm_option,try_max)

%%
%=============================INITIALIZATION===============================

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*a))); 

dataset_specifier = 0;
eval(['simulation_set = att']);
no_clusters = 30;
%--------------------------------------------------------------------------

%--------------------------Choose the Training Dataset---------------------
switch dataset_specifier
    case 0                    
        db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Mixed/CIFAR_10_Gray_Mixed_DB.mat'];  
        db_name_in = 'CIFAR_10_Gray_Mixed_DB';
    case 1
        db_file_in = ['/scratch/amir/Databases/CIFAR_10/Color_Mixed/CIFAR_10_Color_Mixed_DB.mat'];  
        db_name_in = 'CIFAR_10_Color_Mixed_DB';
    case 2
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/CIFAR_10_Gray_Mixed_DB_whitened.mat';
        db_name_in = 'CIFAR_10_Gray_Mixed_DB_whitened';
    case 3
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Mixed/Whitened/CIFAR_10_Color_Mixed_DB_whitened.mat';
        db_name_in = 'CIFAR_10_Color_Mixed_DB_whitened';
    case 4
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Mixed/OMP_5_8_50/CIFAR_10_Color_Mixed_OMP_5_8_50.mat';
        db_name_in = 'CIFAR_10_Mixed_OMP_5_8_50';
    case 5
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/Sparse Filtered/Layer_2/Pooled_3_3/Unnormalized/Pooled_Sparse_Filtered_CIFAR_10_Gray_Mixed_DB_Gray_Class_1_Train_Layer_2.mat';
        db_name_in = 'final_database_vectorized';
    case 6
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Class_6/Train_Set/OMP_5/Pooled/pooled_qudrant_class_6_omp_8_50_2500.mat';
        db_name_in = 'pooled_features';        
    case 7
        db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(1),'/Train_Set/OMP_7/CIFAR_10_train_gray_class_',num2str(1),'_OMP_7_dict_50_patch_8_stand.mat'];        
        db_name_in = ['rec_ims'];                    
    case 10
        db_file_in = ['/scratch/amir/Databases/STL_10/Gray_Mixed/STL_10_Gray_Mixed_DB.mat'];  
        db_name_in = 'STL_10_Gray_Mixed_DB';
    case 11
        db_file_in = ['/scratch/amir/Databases/STL_10/Color_Mixed/STL_10_Color_Mixed_DB.mat'];  
        db_name_in = 'STL_10_Color_Mixed_DB';
    case 12
        db_file_in = '/scratch/amir/Databases/STL_10/Gray_Mixed/Whitened/STL_10_Gray_Mixed_DB_whitened.mat';
        db_name_in = 'STL_10_Gray_Mixed_DB_whitened';
    case 13
        db_file_in = '/scratch/amir/Databases/STL_10/Color_Mixed/Whitened/STL_10_Color_Mixed_DB_whitened.mat';
        db_name_in = 'STL_10_Color_Mixed_DB_whitened';
    case 14    
        db_file_in = ['/scratch/amir/Databases/STL_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/STL_10_train_gray_scale_class_',num2str(target_class),'.mat'];            
        db_name_in = ['STL_10_Gray_DB_class_',num2str(target_class)];                        
    otherwise        
        error('Invalid dataset!')   
end
%--------------------------------------------------------------------------

%-------------------------Load the Training Dataset------------------------
load(db_file_in);
eval(['dataset_learn = ',db_name_in,';']);          
eval(['clear ',db_name_in,';']);                   
[dataset_size,~] = size(dataset_learn); 
%--------------------------------------------------------------------------



%------------------Create Subdirectories if Necessary----------------------
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

destination_folder = [db_file_in(1:i-1),'/Learn_Results/Clustered_Version'];
if (~exist(destination_folder,'dir'))        
    mkdir(destination_folder);    
end
%--------------------------------------------------------------------------

%------------------------Initialize the SSH Connection---------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos'); 
%--------------------------------------------------------------------------

%==========================================================================


%%
%=============================LEARNING PHASE===============================
counter = 0;
simulation_parameters = [];
cluster_size_range = [no_of_PC+41];
for ik = 1:size(cluster_size_range,1)
    n = cluster_size_range(ik,1);                    
    for const_learn = [n-no_of_PC]            
        for alpha0 = [.95]                
            for theta0 =[.005]
                for beta0 = [1]
                    for Q = [0,20]
                        %---------------Quantize if Necessary--------------
                        if (Q > 0)
                        a = matrix2vector(dataset_learn);            
                        mean1 = mean(a);
                        a = a-mean1;
                        [d,c] = hist(abs(a),40000);
                        sigma1 = Q/max(abs(min(a)),max(a));
                        

                        else
                            sigma1 = 0;
                            mean1 = 0;
                        end
                        %--------------------------------------------------
                        
                        %----------Determine the Number of Dominant Princiapl Components-----------
                        [COEFF,SCORE,latent] = princomp(dataset_learn);
                        for i = 1:length(latent)
                            if (sum(latent(1:i))/sum(latent)>.95)                
                                no_of_PC = i;                
                                break;            
                            end    
                        end
                        dataset_projected = SCORE(:,1:no_of_PC)*COEFF(:,1:no_of_PC)';
                        dataset_learn = dataset_projected;
                        %--------------------------------------------------------------------------

                        for learn_itr_max = [2250]
                            current_params = [n,no_clusters,const_learn,alpha0,beta0,theta0,Q,learn_itr_max,no_of_PC];
                            simulation_parameters = [simulation_parameters;current_params];
                            
                            for cluster = 1:no_clusters                       % Execute the learning algorithm for each cluster                                                                
                                counter = counter + 1;
                                fid = fopen([destination_folder,'/job_being_process_cluster_size_',num2str(n),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
                                            num2str(theta0),'_clustere_',num2str(cluster),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.txt'],'r');
    
                                %--Take Into Account the Constraints Being Learned--
                                job_flag = -1;
                                if (fid > -1)
                                    r = fscanf(fid, '%d %d',[2,inf]);
                                    job_flag = r(1);
                                    job_id = r(2);
                                end
                            
                                if (job_flag == 0)
                                    command = ['qstat -x ',num2str(job_id)];
                                    [channel, result]  =  sshfrommatlabissue(channel,command);
                                    ttemp = result{1};
                                    if (isequal(result, {''})) 
                                        job_flag = -1;
                                    else
                                        if (isequal(ttemp(1:21),'qstat: Unknown Job Id'))
                                            job_flag = -1;
                                        end
                                    end
                                elseif (job_flag == 1)
                                    command = ['qstat -x ',num2str(job_id)];
                                    [channel, result]  =  sshfrommatlabissue(channel,command);
                                    ttemp = result{1};
                                    if (isequal(result, {''})) 
                                        job_flag = -1;
                                    else
                                        if (isequal(ttemp(1:21),'qstat: Unknown Job Id'))
                                            job_flag = -1;
                                        end
                                    end
                                end
                            
                                if (job_flag == -1)                                                            
                                    %------------Submit the Job to the Cluster-------------
                                    command = ['cd /scratch/amir/Clustered_Neural/Submit_Cluster;qsub -N "Neur_clustered_',num2str(counter),...
                                        '_learn" -v no_clusters=',num2str(no_clusters),',cluster_size=',num2str(n),',cluster_index=',num2str(cluster),',const_learn=',num2str(const_learn),...
                                        ',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',db_name_in=''"''"''',db_name_in,'''"''"'',db_file_in=''"''"''',db_file_in,...
                                        '''"''"'',learn_itr_max=',num2str(learn_itr_max),',sigma1=',num2str(sigma1),',mean1=',num2str(mean1),...
                                        ',simulation_set=',num2str(simulation_set),',no_of_PCs=',num2str(no_of_PC),',Q=',num2str(Q),' clustered_neural_learn.pbs'];                                                                                
                                    [channel, result]  =  sshfrommatlabissue(channel,command);             
                                    %------------------------------------------------------
                                
    
                                    %----------Check the success of the submission---------
                                    if (isequal(result, {''})) 
                                        display('Unsubmitted learning job!');
                                    else
                                        display('Learning job submitted successfully!');
                                        ttemp = result{1};
                                        job_id = str2num(ttemp(1:6));
                                        fid = fopen([destination_folder,'/job_being_process_cluster_size_',num2str(n),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
                                                num2str(theta0),'_clustere_',num2str(cluster),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.txt'],'w');
                                        fprintf(fid, '%d %d', [0 job_id]);
                                        fclose(fid);
                                    
                                    end
                                elseif (job_flag == 0)
                                    display(['Job has already been submitted with job ID = ',num2str(job_id)]);
                                elseif (job_flag == 1)
                                    display(['Job has already been started with job ID = ',num2str(job_id)]);                               
                                elseif (job_flag == 2)
                                    display('Job has finished successfuly');
                                    pause;
                                else
                                    error('invalid job_flag');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%==========================================================================

%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
%==========================================================================
111
save(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set)],'simulation_parameters','db_file_in','dataset_specifier','db_name_in','no_clusters')
clear 'simulation_set'
clear 'att'
pause

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
        
        if (m < 30)                        
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
                command = ['cd /scratch/amir/Clustered_Neural/Submit_Cluster;qsub -N "adaptive_neur_N_',num2str(N),'K',num2str(K),'C',num2str(l),'_recall" -v N_var=',num2str(N),',K_var=',num2str(K),',L_var=',num2str(L),',max_inst=',num2str(simul_instances_per_job),',err_bit=',num2str(err_bits),',noise_amp=',num2str(max_noise_amp),',alpha0=',num2str(alpha0),',beta0=',num2str(beta0),',theta0=',num2str(theta0),',gamma_BFO=',num2str(gamma_BFO),',gamma_BFS=',num2str(gamma_BFS),',algorithm_option=',num2str(recall_algorithm_option),',try_max=',num2str(try_max),',index=',num2str(index),' clustered_neural_recall.pbs'];                            
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