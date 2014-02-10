% Based on the paper "Complex-Valued Multistate Neural Associative Memory"

function Complex_Neur_Lee_Learn_Recall(N,K,err_bits_range,max_noise_amp,no_of_patterns,no_of_simulated_instance,train_set_index,random_flag)

z_max = 1;
z_min = 0;
max_itr_recall = 5;
bit_error_count = 0;
pattern_error_count = 0;

%==========================LOADING THE DATASET=============================
%-------------------Adjust the Training-Related Materials------------------
load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(train_set_index),'.mat']);   
%--------------------------------------------------------------------------

addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path

dataset_learn = zeros(no_of_patterns,N);
 %------------------Pick a Pattern at Random--------------------
if (random_flag)
    for mu = 1:no_of_patterns    
        dataset_learn(mu,:) = [randi(2,1,N)-1];
    end
else
    for mu = 1:no_of_patterns 
        index = 1+1 + floor(rand*2^K);
        temp = dec2bin((index),K);                                  
        message = zeros(1,K);           % Generate the message from the index                                                         
        for j = 1:K                                    
            message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;                                             
        end    
        dataset_learn(mu,:) = [message*G];
    end
end

S = max(max(abs(dataset_learn)));
phi0 = 2*pi/(S+1);
%==========================================================================


%============================LEARNING PHASE================================
Sigma = exp(1i*phi0*dataset_learn)';
W = Sigma*pinv(Sigma);
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
        fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/Lee/N_',num2str(N),...
            '_Random/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'a+');  
    else
        fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/Lee/N_',num2str(N),'_K_',num2str(K),...
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
