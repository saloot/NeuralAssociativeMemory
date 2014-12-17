%==========================================================================
%*****************FUNCTION: Complex_Neur_Jankowski_Learn*******************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% const_learn: Number of constraints which must be learned during the learning phase
% cluster_index: The index of clusterS for which the learning should be done
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% index: The index of the simulation setup among various random scenarios
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
    N = 300;                                % N is the number of neurons in network.
    K = 150;                                % K is the number of message bits.    
    random_flag = 1;                        % Determines if patterns are drawn from a subspace or generated randomly
    no_of_patterns = 200;                   % The number of patterns in the dataset
    no_of_simulated_instance = 500;         % The number of patterns that are going to be denoised during the recall phase
    max_noise_amp = 1;                      % Maximum value of integer-valued noise added to each bit
    err_bits_range = [0:10];                % The number of bits that will be corrupted initially for the recall phase
end
addpath(genpath('../../Common_Library'));                         % Include the library of common functions in the search path
max_itr_recall = 20;                        % Maximum number of iterations to perform the recall algorithm
if (random_flag)
    mkdir('../../Recall_Results/Jankowski',['N_',num2str(N),'_Random'])
else
    mkdir('../../Recall_Results/Jankowski',['N_',num2str(N),'_K_',num2str(K)])
end
%==========================================================================


for train_set_index = 1:index_max
    
    display(['Performing the recall phase for graph ',num2str(train_set_index),' of the ensemble.'])

%================LOAD THE DATASET AND THE CONNECTIVITY MATRIX==============
    if (random_flag)
        load(['../../Learn_Results/Jankowski/N_',num2str(N)...
            ,'_Random/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat']);
    else
        load(['../../Learn_Results/Jankowski/N_',num2str(N)...
            ,'_K_',num2str(K),'/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat']);    
    end

    S = max(max(abs(dataset_learn)));
    phi0 = 2*pi/(S+1);
%==========================================================================

%==============================RECALL PHASE================================
    for iki = 1:length(err_bits_range)
        err_bits = err_bits_range(iki);
        
        bit_error_count = 0;
        pattern_error_count = 0;
        
        for iji = 1:no_of_simulated_instance
            mu = 1 + floor(rand*no_of_patterns);
            p = dataset_learn(mu,:);     
        
            nois = zeros(1,N);
            pp = 1+floor((N-1)*rand(1,err_bits));                
            for h = 1:err_bits        
                nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
            end        
            x = p + nois;
            x = exp(1i*x*phi0);
    
            for itr = 1:max_itr_recall
            %   for ij = 1:N
                ind = 1 + floor(rand*N);
                h = W*x.';
                %x(ind) = csign(h(ind)*exp(1i*phi0/2),phi0);                       
                x_old = x;
                x = csign(h*exp(1i*phi0/2),phi0);            
                x = x.';            
                if ( norm(abs(x-x_old))/norm(abs(x)) < 1e-5)
                    break
                end
            %   end
            end
        
            if ( sum((abs(x-exp(p*1i*phi0))>1e-3)) > 0)        
                pattern_error_count = pattern_error_count + 1;
                bit_error_count = bit_error_count + sum((abs(x-exp(p*1i*phi0))>1e-3));
            end    
        end 
    
        %=======================STORE THE RESULTS==========================
        BER = bit_error_count/N/iji;
        PER = pattern_error_count/iji;

        if (random_flag)
            fid = fopen(['../../Recall_Results/Jankowski/N_',num2str(N),...
                '_Random/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'a+');  
        else
            fid = fopen(['../../Recall_Results/Jankowski/N_',num2str(N),'_K_',num2str(K),...
                '/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'a+');  
        end
    
        if (fid > -1)
            fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
            fprintf(fid,'\n');
            fclose(fid);
        else
            error('Can not store the results');
        end    
        %=================================================================
    end
end
%==========================================================================


