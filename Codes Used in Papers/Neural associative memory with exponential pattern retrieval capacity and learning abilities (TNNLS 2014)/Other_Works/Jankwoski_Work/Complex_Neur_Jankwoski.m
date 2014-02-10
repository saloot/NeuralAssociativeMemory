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


% function Complex_Neur_Jankwoski

%%
%=============================INITIALIZATION===============================

N = 800;
K = 400;

random_flag = 1;

err_bits_range = '[';

e_max = 10;

for e = 0:e_max   
    err_bits_range = [err_bits_range,num2str(e)];    
    if (e < e_max)                    
        err_bits_range = [err_bits_range,','];                
    else        
        err_bits_range = [err_bits_range,']'];        
    end   
end

max_noise_amp = 1;
no_of_patterns = 200;
no_of_simulated_instance = 5000;
simul_instances_per_job = 5000;

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

%--------------------Initialize the SSH Connection-------------
if (~exist('channel','var'))                
    channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');                                             
end
%--------------------------------------------------------------

%==========================================================================
    

%%
%=============================RECALL PHASE=================================

%---------Perform the Recall Phase for the Finished Learning Jobs----------
for l =1:index_max                                % For all simulation scenarios do...
                
    for itr = 1:ceil(no_of_simulated_instances/simul_instances_per_job)   % Simulate the error correction procedure for the given ensemble.
          
        
        %-------------Submit Recall Jobs to the Cluster------------                                                
        command = ['cd /scratch/amir/ITW_Journal/Submit_Cluster;qsub -N "Jankowski_neur_N_',num2str(N),...
                '_recall" -v N=',num2str(N),',K=',num2str(K),',max_noise_amp=',num2str(max_noise_amp),...
                ',err_bits_range=''"''"''',err_bits_range,'''"''"'',no_of_patterns=',num2str(no_of_patterns),',no_of_simulated_instance=',num2str(simul_instances_per_job),',train_set_index=',num2str(l),...
                ',random_flag=',num2str(random_flag),' jankowski_neural_recall.pbs'];                                            
                        
        [channel, result]  =  sshfrommatlabissue(channel,command);     
                                        
        if (isequal(result, {''}))                                                                         
            error('Unsubmitted recall job!');                                                      
        else            
            display('Recall job submitted successfully!');
        end        
        %----------------------------------------------------------
    end   
end   
%--------------------------------------------------------------------------

%==========================================================================


%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
clear channel
%==========================================================================