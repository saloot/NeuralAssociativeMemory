% Based on the paper "Complex-Valued Multistate Neural Associative Memory"

% function Complex_Neur_Muezzinoglu

N = 40;
K = 400;
train_set_index = 1;
z_max = 1;
z_min = 0;
max_itr_recall = 5;
no_of_simulated_instance = 10;
bit_error_count = 0;
pattern_error_count = 0;

%==========================LOADING THE DATASET=============================
%-------------------Adjust the Training-Related Materials------------------
% load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(train_set_index),'.mat']);   
%--------------------------------------------------------------------------

no_of_patterns = min(length(mu_list),6);
dataset_learn = zeros(no_of_patterns,N+1);
 %------------------Pick a Pattern at Random--------------------
for mu = 1:no_of_patterns 
    index = mu+1;
    temp = dec2bin(mu_list(index),K);                                  
    message = zeros(1,K);           % Generate the message from the index                                                         
    for j = 1:K                                    
        message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;                                             
    end    
%     dataset_learn(mu,:) = [message*G,1];
    dataset_learn(mu,:) = [randi(2,1,N)-1,1];
end

S = max(max(abs(dataset_learn)));
phi0 = 2*pi/(S+1);

%==========================================================================



% %=============================LEARNING PHASE===============================
% W = zeros(N,N);
% A = [];
% for mu = 1:no_of_patterns                        
%     x = dataset_learn(mu,:); 
%     y = x;
%     for ij = 1:N
%         y(ij) = x(ij) + 1;
%         if (y(ij)>S-1)
%             y(ij) = 0;
%         end
%         v = [];
%     
%         for i = 1:N+1
%             for j = i+1:N+1
%                 c = cos(phi0*(x(j)-x(i))) - cos(phi0*(y(j)-y(i)));
%                 s = -sin(phi0*(x(j)-x(i))) + sin(phi0*(y(j)-y(i)));
%                 v = [v,c,s];
%             end
%         end
%         A = [A;v];       
% 
%         y(ij) = x(ij) - 1;
%         
%         if (y(ij)<0)
%             y(ij) = S-1;
%         end
%         
%         v = [];    
%         for i = 1:N+1
%             for j = i+1:N+1
%                 c = cos(phi0*(x(j)-x(i))) - cos(phi0*(y(j)-y(i)));
%                 s = -sin(phi0*(x(j)-x(i))) + sin(phi0*(y(j)-y(i)));
%                 v = [v,c,s];
%             end
%         end
%         A = [A;v];       
%         
%     end
% end
% %--------------------------------------------------------------
    

% 
% SOLVE Aq> 0 
q = linprog(zeros(1,N*(N+1)),A,zeros(1,2*no_of_patterns*N));

W_hat = zeros(N+1,N+1);
for i = 1:N
    for j = i+1:N
        W_hat(i,j) = q(2*j-3) + 1i * q(2*j-2);
        W_hat(j,i) = q(2*j-3) - 1i * q(2*j-2);
% x = exp(1i*x*phi0);
%     W = W + x.'*conj(x);
% W = W/N;
% %==========================================================================
% 
% 
% 
% %==============================RECALL PHASE================================
% for mu = 1:no_of_patterns
%     p = dataset_learn(mu,:);     
%     
%     nois = zeros(1,N);
%     x = p + nois;
%     x = exp(1i*x*phi0);
%     
%     for itr = 1:max_itr_recall
%         for ij = 1:N
%             ind = 1 + floor(rand*N);
%             h = W*x.';
%             x(ind) = csign(h(ind)*exp(1i*phi0/2),phi0);
% %             x = x.';
%         end
%     end
%         
%     if ( sum(sign(abs(x-exp(p*1i*phi0)))) > 0)        
%         pattern_error_count = pattern_error_count + 1;
%         bit_error_count = bit_error_count + sum(sign(abs(x-exp(p*1i*phi0))));
%     end    
%     
% end
% %==========================================================================
% 
% 
% 
% %==========================PLOTTING THE RESULTS============================
% BER = bit_error_count/N/no_of_patterns
% PER = pattern_error_count/no_of_patterns
% %==========================================================================