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



%%
%=============================INITIALIZATION===============================

if (~exist('initialization_done','var'))    % If already not initialized by the GUI..

    %------------------------Simulation Variables--------------------------
    if (~exist('initialization_done','var'))    % If already not initialized by the GUI..
    K = 150;                            % Number of message bits
    N = 300;                            % Number of pattern neurons in the network
    Q = 8;                              % Number of quantization levels
    const_to_learn = 150;               % Number of contraints to learn over the patterns
    random_dataset_flag = 1;            % If 1, this flag tells the code to use the dataset generated by the file "neural_initialization.m". If 0, it will read the dataset from the file specified by the user
    
    no_simulated_instances = 1000;      % The number of patterns that are going to be denoised during the recall phase
    max_noise_amp = 1;                  % Maximum value of integer-valued noise added to each bit
    err_bits_range = [0:10];            % The number of bits that will be corrupted initially for the recall phase
    gamma_BFO = 0.95;                   % The update threshold for the Original Bit flipping algorithm
    gamma_BFS = 0.95;                   % The update threshold for the Simplified Bit flipping algorithm
    theta0 = 0.02;                      % The initial sparisty threshold
    alpha0 = 0.9;                       % The initial learning rate
    beta0 = 0.8;                        % The sparsity penalty
    nu = 0.025;                         % update threshold for the constraint neurons during the recall phase    
    index_in = 1;                       % Index of the random graph in the considered ensemble (for random_dataset_flag =1)
    %----------------------------------------------------------------------

    %-----------------------Other Initializations--------------------------
    addpath(genpath('./Common_Library'));                                  % Include the library of common functions in the search path
    mkdir(['../Recall_Results'],['N_',num2str(N),'_K_',num2str(K)]);        % Create a specific folder for the current N and K
    
    x_min = 0;                          % The minimum value of the pattern nodes during the recall phase
    x_max = 20;                         % The maximum value of the pattern nodes during the recall phase

    a=clock;                                % Initialize the seed for random number generation with the clock value.
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*a))); 
    %----------------------------------------------------------------------
end

%==========================================================================

%========================LOAD THE LEARNED DATASET==========================

%--------------Load the Random Synthetic Training Dataset------------------
if random_dataset_flag == 1
    fid = fopen(['./Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(index_in),'.mat'], 'r');                        % The path towards the dataset
    if (fid > -1 )
        load(['../Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(index_in),'.mat']);
    else    
        error('I can not find the learning dataset!');    
    end
    
    
    %...........................Create the Dataset.........................
    dataset_learn = zeros(length(mu_list),N);
    [C,~] = size(dataset_learn);
    for mu = 1:length(mu_list)
        
        %~~~~~~~~~~~~~~~~~~~~~Pick a Pattern at Random~~~~~~~~~~~~~~~~~~~~~
        index = 1+floor((length(mu_list)-1)*rand);
        temp = dec2bin(mu_list(index),K);                      
        message = zeros(1,K);           % Generate the message from the index                                             
        for j = 1:K                    
            message(j) = (temp(j) - 48);
        end            
        dataset_learn(mu,:) = message*G;
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end
    %......................................................................
    
%--------------------------------------------------------------------------    

%-----------------Load the Real-Valued Training Dataset--------------------
else
    fid = fopen('../Initialization_Files/Communes-Yes.txt', 'r');                        % The path towards the dataset
    if (fid > -1)           
        dataset_learn = fscanf(fid, '%f',[no_communes,no_votations]); 
        dataset_learn = dataset_learn';
        fclose(fid);                                    
        [C,~] = size(dataset_learn);            
    else    
        error('I can not find the learning dataset!');    
    end
    
    
    %.....................Quantize the Training Dataset....................
    max_dataset = max(max(dataset_learn));    
    dataset_learn = round(Q*dataset_learn/max_dataset);    
    %......................................................................

    %.....................Transform the Dataset to Bits....................
    q = (1+floor(log(Q)/log(2)));
    dataset_learn_binary = zeros(C,N*q);
    for i = 1:C
         dataset_learn_binary(i,:) = reshape(de2bi(dataset_learn(i,:),q),[1,N*q]);
    end
    N = N*q;
    dataset_learn = dataset_learn_binary;
    %......................................................................

end
%--------------------------------------------------------------------------

%==========================================================================


%%
%=======================READ THE WEIGHT MATRIX=============================

%-------------Load the Weight Matrix from the File--------------    
fid = fopen(['../Learn_Results/N_',num2str(N),'_K_',num2str(K),'/W_alpha_',num2str(alpha0),'_theta_',num2str(theta0),'_index_',num2str(index_in),'.txt'], 'r');
if (fid == -1)
    error('Learning phase is not complete yet!');
end
W = fscanf(fid, '%f',[N,const_to_learn]);
W = W';
fclose(fid);

%----------Check Whether the Initial Learning Phase Is Done------------
[m1,~] = size(W);
    if (m1 < const_to_learn-10)
        display('Learning phase is not complete yet!');
    else    
        display(' ');
        display('---------------------------------------------------------');
        display('Learning phase finished successfully');
        display('---------------------------------------------------------');
        display(' '); 
    end
%--------------------------------------------------------------------------

%==========================================================================


%==========================PERFORM THE RECALL PHASE========================
for error_itr = 1:length(err_bits_range)
    err_bits = err_bits_range(error_itr);
    
    error_count_WTA = 0;                % Pattern error rate for the Winnter-Take-ALL (WTA) algorithm
    error_count_BFO = 0;                % Pattern error rate for the Original Bit Flipping (BFO) algorithm
    error_count_BFS = 0;                % Pattern error rate for the Simplified Bit Flipping (BFS) algorithm
    bit_error_count_WTA = 0;            % Bit error rate for the Winnter-Take-ALL (WTA) algorithm
    bit_error_count_BFO = 0;            % Bit error rate for the Original Bit Flipping (BFO) algorithm
    bit_error_count_BFS = 0;            % Bit error rate for the Simplified Bit Flipping (BFS) algorithm

    for net_simul_itr = 1:no_simulated_instances   
        
        %---------------------------Generate Noise-------------------------
        nois = zeros(1,N);          % This is the noise added to the whole pattern of length N*L_in        
        pp = 1+floor((N-1)*rand(1,err_bits));                
        for h = 1:err_bits        
            nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
        end        
        %------------------------------------------------------------------
                
        %--------------------Generate the Pattern Index--------------------
        p = 1+floor((C-1)*rand);                 % Pick a pattern index at random                   
        pattern = dataset_learn(p,:);
        %------------------------------------------------------------------                                

        %----------Initialize the Network with a Noisy SubPattern----------
        x = pattern+nois;                                       % Initialize the network with a nosiy version of the subpattern               
        %------------------------------------------------------------------
           
        %---------------------Iterate Until Convergence--------------------
        [x_WTA,x_BFO,x_BFS] = recall_step(W,x,N,m1,gamma_BFO,gamma_BFS,max_noise_amp,err_bits,x_min,x_max,nu);    
        %------------------------------------------------------------------
           
        %--------------Calculate Subpatterns Error Rate--------------------
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
        %------------------------------------------------------------------
    
        %-----------------------Display Progress---------------------------
         if (mod(net_simul_itr, 1000) == 0)         
             display(' ');
             display(['Number of errors = ',num2str(err_bits)]);
             display(['Pattern error rate WTA = ',num2str(error_count_WTA/net_simul_itr)])
             display(['Pattern error rate BFO = ',num2str(error_count_BFO/net_simul_itr)])
             display(['Pattern error rate BFS = ',num2str(error_count_BFS/net_simul_itr)])        
             display(' ');
         end
         %-----------------------------------------------------------------
    end    

    %----------------Transform Error Count to Error Rate-------------------
    PER_WTA = error_count_WTA/net_simul_itr;
    PER_BFO = error_count_BFO/net_simul_itr;
    PER_BFS = error_count_BFS/net_simul_itr;
    BER_WTA = bit_error_count_WTA/net_simul_itr/N;
    BER_BFO = bit_error_count_BFO/net_simul_itr/N;
    BER_BFS = bit_error_count_BFS/net_simul_itr/N;
    %----------------------------------------------------------------------

  
    %---------------------------SAVE THE RESULTS---------------------------
    fid = fopen(['../Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER_WTA,BER_WTA);
    fprintf(fid,'\n');
    fclose(fid);
    
    fid = fopen(['../Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER_BFO,BER_BFO);
    fprintf(fid,'\n');
    fclose(fid);

    fid = fopen(['../Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER_BFS,BER_BFS);
    fprintf(fid,'\n');
    fclose(fid);
    %----------------------------------------------------------------------
end
%==========================================================================

    
