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


function Complex_Neur_Jankowski_Recall(N,K,err_bits_range,max_noise_amp,no_of_patterns,no_of_simulated_instance,train_set_index,random_flag)

%============================INITIALIZATIONS===============================
max_itr_recall = 20;


addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
%==========================================================================


%================LOAD THE DATASET AND THE CONNECTIVITY MATRIX==============
if (random_flag)
    load(['/scratch/amir/ITW_Journal/Learn_Results/Jankowski/N_',num2str(N)...
        ,'_Random/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat']);
else
    load(['/scratch/amir/ITW_Journal/Learn_Results/Jankowski/N_',num2str(N)...
        ,'_K_',num2str(K),'/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(train_set_index),'.mat']);    
end

S = max(max(abs(dataset_learn)));
phi0 = 2*pi/(S+1);
%==========================================================================

%==============================RECALL PHASE================================
for iki = 1:(length(err_bits_range)-1)/2
    err_bits = str2num(err_bits_range(2*(iki)));    
    
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
    
    %=========================STORE THE RESULTS============================
    BER = bit_error_count/N/iji;
    PER = pattern_error_count/iji;

    if (random_flag)
        fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/Jankowski/N_',num2str(N),...
            '_Random/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'a+');  
    else
        fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/Jankowski/N_',num2str(N),'_K_',num2str(K),...
            '/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'a+');  
    end
    
    if (fid > -1)
        fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
        fprintf(fid,'\n');
        fclose(fid);
    else
        error('Can not store the results');
    end
    %======================================================================
    
end
%==========================================================================


