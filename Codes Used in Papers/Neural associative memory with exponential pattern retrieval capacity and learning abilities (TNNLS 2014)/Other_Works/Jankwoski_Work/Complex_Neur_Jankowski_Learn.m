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
N = 800;
K = 400;
try_index_max = 50;
z_max = 1;
z_min = 0;
dataset_learn = zeros(no_of_patterns,N);
no_of_patterns = 200;
random_flag = 1;
%==========================================================================




for train_set_index = 1:try_index_max
    
    %========================LOADING THE DATASET===========================
    
    %-----------------Adjust the Training-Related Materials----------------
    load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),...
        '/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(train_set_index),'.mat']);   
    %----------------------------------------------------------------------



    %----------------------Pick a Pattern at Random------------------------
    if (random_flag)
        for mu = 1:no_of_patterns 
           dataset_learn(mu,:) = randi(2,1,N)-1;
        end
    else
        for mu = 1:no_of_patterns 
            index = mu+1;
            temp = dec2bin(mu_list(index),K);                                  
            message = zeros(1,K);           % Generate the message from the index                                                         
            for j = 1:K                                    
                message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;                                             
            end    
            dataset_learn(mu,:) = message*G; 

        end
    end

    S = max(max(abs(dataset_learn)));
    phi0 = 2*pi/(S+1);
    %======================================================================


    %===========================LEARNING PHASE=============================
    W = zeros(N,N);
    for mu = 1:no_of_patterns                        
        x = dataset_learn(mu,:); 
        x = exp(1i*x*phi0);
        W = W + x.'*conj(x);            
    end

    W = W/N;
    %======================================================================

    %==========================STORE THE RESULTS===========================
    if (random_flag)
        save(['/scratch/amir/ITW_Journal/Learn_Results/Jankowski/N_',num2str(N)...
            ,'_Random/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat'], 'W','dataset_learn');
    else
        save(['/scratch/amir/ITW_Journal/Learn_Results/Jankowski/N_',num2str(N)...
            ,'_K_',num2str(K),'/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat'], 'W','dataset_learn');
    end    
    %======================================================================
    
end


