%==========================================================================
%***************FUNCTION: neural_recall_ITW_journal************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% no_simulated_instances: The number of patterns considered in the recall phase.
% err_bits: Number of initial noisy nodes.
% max_noise_amp: The maximum amplitude of noise during the recall process
% index_in: The index of the simulation setup among various random scenarios
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a bipartite neural associative
% memory and performs the recall phase. More specifically, 
% the function first reads the weight matrix found in the learning phase. 
% Then the function performs the recall phase with three different
% algorithms the winner-take-all, the original bit flipping (introduced in
% our ITW 2011) paper and the simplified bit-flipping (introduced in our
% ISIT 2012 paper).
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function neural_recall_ITW_journal(N,K,no_simulated_instances,err_bits,max_noise_amp,index_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS)

%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_parameters_v1_N_',num2str(N),'_K_',num2str(K),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
mkdir(['/scratch/amir/ITW_Journal/Recall_Results'],['N_',num2str(N),'_K_',num2str(K)]);        % Create a specific folder for the current N and K

error_count_WTA = 0;
error_count_BFO = 0;
error_count_BFS = 0;
bit_error_count_WTA = 0;
bit_error_count_BFO = 0;
bit_error_count_BFS = 0;

first_itr_error_count_WTA = 0;
first_itr_error_count_BFO = 0;
first_itr_error_count_BFS = 0;
first_itr_bit_error_count_WTA = 0;
first_itr_bit_error_count_BFO = 0;
first_itr_bit_error_count_BFS = 0;

y_min = 0;
y_max = 20;

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------


%-------------------Adjust the Training-Related Materials------------------
load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(index_in),'.mat']);   
%--------------------------------------------------------------------------

%==========================================================================


%%
%=======================READ THE WEIGHT MATRIX=============================

%-------------Load the Weight Matrix from the File--------------    
fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(index_in),'.txt'], 'r');
if (fid == -1)
    error('Learning phase is not complete yet!');
end
W = fscanf(fid, '%f',[N,N_const]);
W = W';
fclose(fid);

%----------Check Whether the Initial Learning Phase Is Done------------
[m1,~] = size(W);
    if (m1 < N_const-10)
        display('Learning phase is not complete yet!');
        m1
    else    
        display(' ');
        display('---------------------------------------------------------');
        display('Learning phase finished successfully');
        display('---------------------------------------------------------');
        display(' '); 
    end
%--------------------------------------------------------------------------

%----------------------Calculate the Special Answers-----------------------
x_s = 1+floor(2*rand(1,N));
b = W*x_s';        
%--------------------------------------------------------------------------

%==========================================================================


for net_simul_itr = 1:no_simulated_instances   % Simulate the error correction procedure for the given ensemble.                                
        
    %-----------------------------Generate Noise---------------------------
    nois = zeros(1,N);          % This is the noise added to the whole pattern of length N*L_in        
    pp = 1+floor((N-1)*rand(1,err_bits));                
    for h = 1:err_bits        
        nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
    end        
    %----------------------------------------------------------------------                            
                
    %----------------------Generate the Pattern Index----------------------
    p = 1+floor((pattern_learn_number-1)*rand);                 % Pick a pattern index at random               
    temp = dec2bin(mu_list(index),K);                          
    message = zeros(1,K);           % Generate the message from the index                                                 
    for j = 1:K                        
        message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;                                                     
    end    
    pattern = message*G;   
%     pattern = x_s+pattern;                        % Generate the pattern from the message and the special answer
    %----------------------------------------------------------------------                                
                                                                                                                
    %------------Initialize the Network with a Noisy SubPattern------------
    x = pattern+nois;                                       % Initialize the network with a nosiy version of the subpattern               
    %----------------------------------------------------------------------
           
    %-----------------------Iterate Until Convergence----------------------
    [x_WTA,x_BFO,x_BFS,results_first_itr] = recall_step_journal(W,x,b,N,m1,gamma_BFO,gamma_BFS,max_noise_amp,err_bits,y_min,y_max);    
    %----------------------------------------------------------------------
           
    %-----------------Calculate Subpatterns Error Rate---------------------
    if (norm(x_WTA-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                       
        error_count_WTA = error_count_WTA+1;
        bit_error_count_WTA = bit_error_count_WTA+sum(abs(sign(x_WTA-pattern)));
    end
    
    if (norm(x_BFO-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.
        error_count_BFO = error_count_BFO+1;
        bit_error_count_BFO = bit_error_count_BFO+sum(abs(sign(x_BFO-pattern)));
    end
    if (norm(x_BFS-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.
        error_count_BFS = error_count_BFS+1;
        bit_error_count_BFS = bit_error_count_BFS+sum(abs(sign(x_BFS-pattern)));
    end
    
    
    if (norm(results_first_itr(1,:)-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                       
        first_itr_error_count_WTA = first_itr_error_count_WTA+1;
        first_itr_bit_error_count_WTA = first_itr_bit_error_count_WTA+sum(abs(sign(results_first_itr(1,:)-pattern)));
    end
    
    if (norm(results_first_itr(2,:)-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.
        first_itr_error_count_BFO = first_itr_error_count_BFO+1;
        first_itr_bit_error_count_BFO = first_itr_bit_error_count_BFO+sum(abs(sign(results_first_itr(2,:)-pattern)));
    end
    if (norm(results_first_itr(3,:)-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.
        first_itr_error_count_BFS = first_itr_error_count_BFS+1;
        first_itr_bit_error_count_BFS = first_itr_bit_error_count_BFS+sum(abs(sign(results_first_itr(3,:)-pattern)));
    end
    %----------------------------------------------------------------------        
    
    %-------------------------Display Progress-----------------------------
     if (mod(net_simul_itr, 1000) == 0)         
         display(' ');
         display(['Error rate WTA = ',num2str(error_count_WTA/net_simul_itr)])
         display(['Error rate BFO = ',num2str(error_count_BFO/net_simul_itr)])
         display(['Error rate BFS = ',num2str(error_count_BFS/net_simul_itr)])        
         display(' ');
     end
     %---------------------------------------------------------------------
    
end    
        
            
    
%------------------Transform Error Count to Error Rate---------------------
PER_WTA = error_count_WTA/net_simul_itr;
PER_BFO = error_count_BFO/net_simul_itr;
PER_BFS = error_count_BFS/net_simul_itr;
BER_WTA = bit_error_count_WTA/net_simul_itr/N;
BER_BFO = bit_error_count_BFO/net_simul_itr/N;
BER_BFS = bit_error_count_BFS/net_simul_itr/N;

first_itr_PER_WTA = first_itr_error_count_WTA/net_simul_itr;
first_itr_PER_BFO = first_itr_error_count_BFO/net_simul_itr;
first_itr_PER_BFS = first_itr_error_count_BFS/net_simul_itr;
first_itr_BER_WTA = first_itr_bit_error_count_WTA/net_simul_itr/N;
first_itr_BER_BFO = first_itr_bit_error_count_BFO/net_simul_itr/N;
first_itr_BER_BFS = first_itr_bit_error_count_BFS/net_simul_itr/N;
%--------------------------------------------------------------------------

%==========================================================================
  
%%
%===========================SAVE THE RESULTS===============================
fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER_WTA,BER_WTA);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER_BFO,BER_BFO);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER_BFS,BER_BFS);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'_first_itr.txt'], 'a+');  
fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,first_itr_PER_WTA,first_itr_BER_WTA);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'_first_itr.txt'], 'a+');  
fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,first_itr_PER_BFO,first_itr_BER_BFO);
fprintf(fid,'\n');
fclose(fid);

fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'_first_itr.txt'], 'a+');  
fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,first_itr_PER_BFS,first_itr_BER_BFS);
fprintf(fid,'\n');
fclose(fid);
%==========================================================================

    
