%==========================================================================
%*****************FUNCTION: Complex_Neur_Jankowski_Learn*******************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function implements the learning phase of a multi-state complex-
% valued neural network proposed by Jankowski et al. in the paper 
% "Complex-Valued Multistate Neural Associative Memory". 

% NOTE: In general, the patterns should be produced randomly. However,
% since we are interested in comparing the performance of their proposed
% algorithm with that of our journal paper, we have chosen a synthetic
% dataset so that the simulations are done on the same dataset. 
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================


%============================INITIALIZATIONS===============================
if (~exist('initialization_done','var'))    
    N = 160;                                % N is the number of neurons in network.
    K = 80;                                 % K is the number of message bits.
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.
    index_max = 50;                         % This is the maximum number of random scenarios generated for simulation
    random_flag = 1;                        % Determines if patterns are drawn from a subspace or generated randomly
    no_of_patterns = 200;                   % The number of patterns in the dataset    
end
dataset_learn = zeros(no_of_patterns,N);

if (random_flag)
    mkdir('../../Learn_Results/Jankowski',['N_',num2str(N),'_Random'])
else
    mkdir('../../Learn_Results/Jankowski',['N_',num2str(N),'_K_',num2str(K)])
end
%==========================================================================



for train_set_index = 1:index_max

%============================LOAD THE DATASET==============================
    %-----------------Load the Patterns from a Subspace--------------------
    if random_flag == 0
        fid = fopen(['../../Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(train_set_index),'.mat'], 'r');                        % The path towards the dataset
        if (fid > -1 )
            load(['../../Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(train_set_index),'.mat']);
        else    
            error('I can not find the learning dataset!');    
        end
        
        %.........................Create the Dataset.......................        
        [C,~] = size(dataset_learn);
        for mu = 1:C
            
            %~~~~~~~~~~~~~~~~~~~Pick a Pattern at Random~~~~~~~~~~~~~~~~~~~
            index = 1+floor((C-1)*rand);
            temp = dec2bin(mu_list(index),K);                      
            message = zeros(1,K);           % Generate the message from the index                                             
            for j = 1:K                    
                message(j) = (temp(j) - 48);
            end            
            dataset_learn(mu,:) = message*G;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end
        %..................................................................
    
    %----------------------------------------------------------------------


    %-------------------------Pick Random Patterns-------------------------
    else
        for mu = 1:no_of_patterns 
           dataset_learn(mu,:) = randi(2,1,N)-1;
        end
    end
    %----------------------------------------------------------------------
    
    S = max(max(abs(dataset_learn)));
    phi0 = 2*pi/(S+1);
    
    
%==========================================================================


%==========================PERFORM THE LEARNING============================    
    W = zeros(N,N);
    for mu = 1:no_of_patterns                        
        x = dataset_learn(mu,:); 
        x = exp(1i*x*phi0);
        W = W + x.'*conj(x);            
    end

    W = W/N;
%==========================================================================

%============================STORE THE RESULTS=============================
    if (random_flag)
        save(['../../Learn_Results/Jankowski/N_',num2str(N)...
            ,'_Random/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat'], 'W','dataset_learn');
    else
        save(['../../Learn_Results/Jankowski/N_',num2str(N)...
            ,'_K_',num2str(K),'/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat'], 'W','dataset_learn');
    end    
%==========================================================================
    
end


