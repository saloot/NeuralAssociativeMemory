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


function clustered_neural_recall(N,K,L,no_simulated_instances,err_bits,max_noise_amp,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,try_max,index)

%=============================INITIALIZATION===============================


%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);           

% load('/scratch/amir/Clustered_Neural/Databases/Spoken_Words/DB/FFT/wav_db.mat');
% load('/scratch/amir/Clustered_Neural/Databases/Spoken_Words/DB/FFT/wav_db_quant_mapped.mat');
% [pattern_learn_number,~] = size(wav_db);
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
mkdir(['/scratch/amir/Clustered_Neural/Recall_Results'],['N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)]);        % Create a specific folder for the current N and K
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 

error_count = 0;
bit_error_count = 0;

y_min = -10;
y_max = 100;
%--------------------------------------------------------------------------

%==========================================================================



%%
%===============================MAIN LOOP==================================
for net_simul_itr = 1:no_simulated_instances                % Simulate the error correction procedure for the given ensemble.                                
        
    %-----------------------------Generate Noise---------------------------
    nois = zeros(1,N_tot);                                  % This is the noise added to the whole pattern of length N*L_in        
    pp = 1+floor((N_tot-1)*rand(1,err_bits));                
    for h = 1:err_bits        
        nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
    end        
    %----------------------------------------------------------------------                            
                
    %-----------------------Generate the Pattern---------------------------
    mu = 1+floor((pattern_learn_number-1)*rand);            % Pick a pattern index at random               
    temp = dec2bin(mu_list(1,mu),KK);                  
    randindex = mu_list(2:KK+1,mu);    
    message = zeros(1,K_tot);                               % Generate the message from the index                        
           
    for j = 1:KK    
        message(randindex(j)) = (z_max-z_min)*(temp(j) - 48)+z_min;                                     
    end            
    pattern = message*G_tot;  
%     pattern = wav_db(mu,:);
    %----------------------------------------------------------------------                                
                                                                                                                
    %------------Initialize the Network with a Noisy SubPattern------------
    x = pattern+nois;                                       % Initialize the network with a nosiy version of the subpattern               
    %----------------------------------------------------------------------
            
        
    %------------------------Perform the Recall Step-----------------------
    success_flag = 0;    
    try_itr = 0;
    
    while ((try_itr < try_max) && (success_flag == 0))
        success_flag = 1;
        try_itr = try_itr+1;
        for l = 1:L
            
            
            %---------------------Read the Weight Matrix-------------------
            load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);
            n = length(index_l);
            fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt'], 'r');
            if (fid == -1)
                l
                error('No weight matrix!');
            end
            %['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt']
            W = fscanf(fid, '%f',[n,inf]);        
            W = W';      
            [m,~] = size(W);
            fclose(fid);
            
            if (l == 55)
                111;
            end
            x_temp = zeros(1,length(index_l));
            pattern_temp = zeros(1,length(index_l));
            for iji = 1:length(index_l)
                pattern_temp(iji) = pattern(index_l(iji));
                x_temp(iji) = x(index_l(iji));
            end    
            %--------------------------------------------------------------    
            
            %---------------------Iterate Until Convergence----------------            
            [x_out] = recall_step_clustered(W,x_temp,n,gamma_BFO,gamma_BFS,max_noise_amp,err_bits,y_min,y_max,recall_algorithm_option);    
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

    
