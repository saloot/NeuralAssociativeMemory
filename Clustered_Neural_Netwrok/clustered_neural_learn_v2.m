%==========================================================================
%********************FUNCTION: clustered_neural_learn**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% const_learn: Number of constraints which must be learned during the learning phase
% cluster_index: The cluster for which the learning should be done
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% index: The index of the simulation setup among various random scenarios
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function implements the idea mentioned in our journal paper to learn 
% a sparse vector W which is orthogonal to a set of given patterns. 

% The code starts by reading patterns at random and adjusting the weights 
% accordingly such that in the end we have Wx = 0, for all data vectors x 
% in the training set. The learning approach which is very similar (almost 
% identical to) the one proposed in the paper "The subspace learning 
% algorithm as a formalism for pattern recognition and neural networks" 
% plus a penalty term to make the results sparse. 

% The code stops once a maximum number of iterations is passed or
% convergence is reached. Once finished, the learnt vectors are added to
% the end of the appropriatefile. 
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================

function clustered_neural_learn_v2(N,K,L,const_learn,cluster_index,alpha0,beta0,theta0,index,quantization_flag)


%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
% load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
% mkdir(['/scratch/amir/Clustered_Neural/Learn_Results'],['N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)]);        % Create a specific folder for the current N and K

b=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*b))); 

learn_itr_max = 4000;
cos_phi_thr = 0.1;
max_y_threshold = .05;
min_W = .005;
fraction_of_convergence = .99;
%--------------------------------------------------------------------------

%-------------------Construct the Pattern Generator Matrix-----------------
% load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster_index),'.mat']);
% G = [];
% for i = 1:length(index_l)
%     G = [G,G_tot(:,index_l(i))];              
% end    

%------------------------------Quantization--------------------------------
% partition = [0.0001,.00015,.0002,.00025,.0003,.00035,.0004,.00045,.0005,.0006,.0007,.0008,.0009,.001,.002,.004,.006,.008,.01];
% partition = [0.0001,.00015,.0002,.00025,.0003,.00035,.0004,.00045,.0005,.0006,.0007,.0008,.0009,.001,.002,.0025,.003,.0035,.004,.0045,.005,.0055,.006,.0065,.007,.0075,.008,.0085,.009,.0095,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1];
% partition = [0.0001,.00015,.0002,.00025,.0003,.00035,.0004,.00045,.0005,.0006,.0007,.0008,.0009,.001,.002,.0025,.003,.0035,.004,.0045,.005,.0055,.006,.0065,.007,.0075,.008,.0085,.009,.0095,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.12,.14,.16,.18]
% codebook = zeros(1,length(partition)+1);
% codebook(2:length(partition)+1) = partition;
%--------------------------------------------------------------------------
% n = length(index_l);

% load('/scratch/amir/Clustered_Neural/Databases/Spoken_Words/DB/FFT/dataset_learn.mat');
% load('/scratch/amir/Clustered_Neural/Databases/Spoken_Words/DB/FFT/wav_db_quant.mat');

% if (quantization_flag == 2)
%     load('/scratch/amir/Databases/Spoken_Words/DB_new/FFT/wav_db_fft_quantized_indices_coarse.mat');
%     dataset_learn = wav_db_FFT_quantized_indices;
% elseif (quantization_flag == 1)
%     load('/scratch/amir/Databases/Spoken_Words/DB_new/FFT/wav_db_fft_quantized.mat');
%     dataset_learn = wav_db_FFT_quantized;
% else
%     load('/scratch/amir/Databases/Spoken_Words/DB_new/FFT/wav_db_fft.mat');
%     dataset_learn = wav_db_fft;
% end
% dataset_learn = dataset_learn/0.0001;\

[pattern_learn_number,N_tot] = size(dataset_learn);
%--------------------------------------------------------------------------

%-------------Update the Number of Constraints in the Process--------------

% fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r'); 
% if (fid > -1)
%     current_count = max(fscanf(fid, '%d'),0);            
%     fclose(fid);    
% else    
    current_count = 0;            
% end        
% 
% fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r');
% if (fid > -1)           
%     W = fscanf(fid, '%f',[n,inf]);                
%     W = W';            
%     fclose(fid);                                    
%     [m,~] = size(W);            
% else    
%     m = 0;    
% end
    
m = 0;
N_const_max = 100;
n = 1024;
d_max = 18;
const_to_learn = min(const_learn,N_const_max-m-current_count);
% const_to_learn = min(const_learn,N-(K/2)-m);


% fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'w'); 
% fprintf(fid,'%d',current_count+const_to_learn);
% fclose(fid);
%--------------------------------------------------------------------------

%==========================================================================
    
% W_tot = [];
%%
%=============================LEARNING PHASE===============================
for itr = 1: const_to_learn

    alpha1 = alpha0;
    theta1 = theta0;
    %-----------------------Read Current Learned Matrix--------------------
%     fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r');
%     if (fid > -1)           
%         W_tot = fscanf(fid, '%f',[n,inf]);                
%         W_tot = W_tot';            
%         fclose(fid);                                    
        [m,~] = size(W_tot);            
%     else    
%         m = itr;    
%     end
    
%     if (m>4)
%         index_w = zeros(1,n);
        if (m >= 2)
            s = sum(abs(sign(W_tot)));
        else
            s = zeros(1,n);
        end
        
%         for i = 1:n
%             pp = rand;
%             if (s(i)/d_max <= pp)
%                 index_w(i) = 1;
%             end
%         end
%     else
        index_w = ones(1,n);
%     end
%     if (m > 4)
%         index_w = (sum(abs(sign(W_tot))) < 6);
%         if (sum(index_w) == 0)
%             index_w = (sum(abs(sign(W_tot)))/m < 0.25);
%         end
%     else
%         index_w = ones(1,n);
%     end
    
    n_new = sum(index_w);
    
    if (n_new/n < .3)
        theta1 = 0;
        alpha1 = alpha0/3;
%         max_y_thr = .005;
    else
        theta1 = theta0*n_new/n;
%         max_y_thr = .001;
    end
    %----------------------------------------------------------------------

    %------------------Randomly Initialize the Weight Matrix---------------
    W = random_vector(n_new,2*round(log(n_new)));
    W = W';
    W = rand(n_new,1);
    W = W/norm(W);
    %----------------------------------------------------------------------
    
    %---------------------Adjust Simulation Parameters--------------------- 
    theta = theta1;   
    theta_itr = 1;
    cost_per_itr_av = [];            
    alph_itr = 0;    
    %----------------------------------------------------------------------

    %--------------------------Main Learning Loop--------------------------           
    for learn_itr = 1:learn_itr_max                         % Repeat the steps below several times in order to enhance the quality of the learning phase        
        if (mod(learn_itr,50)==0)
%             W = random_vector(n_new,2*round(log(n_new)));
%             W = W';
            W = rand(n_new,1);
            W = W/norm(W);
            theta_itr = theta_itr + 1;                
            theta = theta1/(2+theta_itr);
            alpha1 = alpha0/theta_itr;
        end
        temp_cost2 =0;
        temp_cost = 0;        
        max_cos_phi = 0;
        max_y = 0;
        
%         alph0 = max(50*alpha1/(50+(learn_itr)),0.0005);
        if (mod(learn_itr,202) == 1)            
            alph_itr = alph_itr + 1;
            alph0 = max(alpha1/alph_itr,0.0005);        
        end
        continue_flag = 0;
        if (learn_itr>1)
            last_cost_one_itr = cost_one_itr;
        end
        cost_one_itr = 0;        
        cost_short_term =0;
        for ii = 1:pattern_learn_number
            
            %------------------Pick a Pattern at Random--------------------
            mu = 1+floor((pattern_learn_number-1)*rand); 
%             temp = dec2bin(mu_list(1,mu),KK);      
%             randindex = mu_list(2:KK+1,mu);
%             message = zeros(1,K_tot);                       % Generate the message from the index                        
% 
%             for j = 1:KK
%                 message(randindex(j)) = (z_max-z_min)*(temp(j) - 48)+z_min;                             
%             end    

%             x = message*G;
            x_temp2 = dataset_learn(mu,:);
%             x_temp2 = [];
%             for i = 1:length(index_l)
%                 x_temp2 = [x_temp2,x_temp(index_l(i))];
%             end
            deg=[];
            x = [];
            for i = 1:n
                if (index_w(i) > 0)
                    x = [x,x_temp2(i)];
                    deg=[deg,s(i)];
                end
            end
            %--------------------------------------------------------------
        
            %---------------------Update the Weight Vector-----------------
            y = x*W;                 % Calculate the projection of x on the vector W(j,:).
            
            if (norm(x)>0)
                alph = alph0/norm(x)^2;
            else
                alph = 0;
            end
            if (norm(W) > .001)                                
%                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,0.05);W.*(1-(tanh(100*W(j,:).^2)).^2)
%                  W = soft_threshold(W - alph*y*( x' - (y*W/(norm(W))^2)),theta0);%-beta0*W.*(1-(tanh(100*W.^2)).^2);
%                  W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta)-.05*W.*(1-(tanh(100*W.^2)).^2).*(1+tanh(100*(deg-d_max)))';
%                 W = soft_threshold(W - alph*y*(x'-y*W)+ alph*(1-norm(W)^2)*W,theta); 
%                  
            end            
            
            W = soft_threshold(W,theta);                       % Set entries less than 0.0001 to zero.
            %--------------------------------------------------------------
        
%             beta0 = max(beta0 + 0.00001*(sum(abs(sign(W)))-theta0),0);
            
            %------------------Check for Numerical Errors------------------
            if (sum(isnan(W))>0)    
%                 W = random_vector(n_new,2*round(log(n_new)));
%                 W = W';
                W = rand(n_new,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
                display('Error Nan');        
                break;                                    
            end
            
            if (norm(W) < .001)                
                display('Error zero!');
%                 W = random_vector(n_new,2*round(log(n_new)));
%                 W = W';
                    W = rand(n_new,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
            end
            %--------------------------------------------------------------

            %-----------------Update the Simulation Costs------------------
            cost_one_itr = cost_one_itr + abs(x*W);
            cost_short_term = cost_short_term+norm(x*W);
            temp_cost = temp_cost + y^2;
            temp_cost2 = temp_cost2+y^2;
            if (norm(x)>0)
                cos_phi = abs(y)/(norm(W)*norm(x));
            else
                cos_phi = 1;
            end
            
            if (abs(x*W)>max_y)
                max_y = abs(x*W);
            end
                        
            if ((cos_phi> max_cos_phi)&(norm(x)>0))
                max_cos_phi = cos_phi;                
            end                                        
            %--------------------------------------------------------------
                                                           
            
            
        end 
        
        
        %--------------------Update the Simulation Costs-------------------
        cost_per_itr_av = [cost_per_itr_av,cost_short_term];
        cost_short_term = 0;
        %------------------------------------------------------------------                

        W_temp = W;
        W2 = zeros(n,1);
        j = 1;
        for i = 1:n
            if (index_w(i) > 0)
                W2(i) = W_temp(j);
                j = j+1;
            end
        end
        
%         W_temp = W2;
%         W2 = zeros(N_tot,1);        
%         for i = 1:length(index_l)                
%             W2(index_l(i)) = W_temp(i);            
%         end

        
        
        if (sum(abs(dataset_learn*W2)<max_y_threshold)>fraction_of_convergence*pattern_learn_number)
            W2=soft_threshold(W2,min_W);         
        end
        
        ww = W2;
        ww(~ww) = inf;
        
        %--------------------------Display Progress------------------------
        if (mod(learn_itr,2)==0)
            theta_itr = theta_itr + 1;
            theta = 4*theta1/(4+theta_itr);            
            
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);
            display(['Converged count = ',num2str(sum((abs(dataset_learn*W2)<max_y_threshold)))]);
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);
            display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
            display(['Maximum projection = ',num2str(max(abs(dataset_learn*W2)))]);  
            display(['Minum of W = ',num2str(min(abs(ww)))]);  
            display('-----------------------------------');
            111;
            display(' ');             
        end           
        %--------------------------------------------------------------
         
            
        %----------------------Check Convergence---------------------------
%         if ((max_cos_phi < cos_phi_thr)&&(max_y<min(max_y_thr,theta)))            
%         if (theta>0) 
%             a = min(max_y_thr,theta);            
%         else
%             a = max_y_thr;
%         end
%         if (max_y<a)
%             if (surely_flag == 1)
%                 break;
%             else
%                 surely_flag = 1;
%             end
%         end
        ww = W2;
        ww(~ww) = inf;
        if (quantization_flag<2)
            if ((sum((abs(dataset_learn*W2))<max_y_threshold)>=fraction_of_convergence*pattern_learn_number)&&(min(abs(ww))>min_W)&&(norm(W)>.1))
%                 if (surely_flag == 1)
                    break;
%                 else
%                     surely_flag = 1;
%                 end
            end            
        else   
%             if ((sum((abs(dataset_learn*W2))<max_y_threshold)>=fraction_of_convergence*pattern_learn_number)&&(min(abs(ww))>min_W)&&(norm(W)>.1)&&(max(abs(dataset_learn*W2))<min_W-.01))
              if ((sum((abs(dataset_learn*W2))<max_y_threshold)>=fraction_of_convergence*pattern_learn_number)&&(min(abs(ww))>min_W)&&(norm(W)>.1))
%                 if (surely_flag == 1)
                    break;
%                 else
%                     surely_flag = 1;
%                 end
            end
        end
        %------------------------------------------------------------------
        
        

    end  
    %======================================================================


    %==========================SAVE THE RESULTS============================    
%     if ((max_y<max_y_thr)&&(sum(abs(sign(W)))<=n/2)&&(sum(abs(sign(W)))>=4))
%     if (theta>0) 
%         a = min(max_y_thr,theta);        
%     else
%         a = max_y_thr;
%     end
%     if (max_y<a)
    if ((sum((abs(dataset_learn*W2))<max_y_threshold)>=fraction_of_convergence*pattern_learn_number)&&(min(abs(ww))>min_W)&&(norm(W)>.1))
        
        
        W=soft_threshold(W,min_W);         
%         W_temp = W;
%         W = zeros(n,1);
%         j = 1;
%         for i = 1:n
%             if (index_w(i) > 0)
%                 W(i) = W_temp(j);
%                 j = j+1;
%             end
%         end
        W_temp = W;
        W2 = zeros(n,1);
        j = 1;
        for i = 1:n
            if (index_w(i) > 0)
                W2(i) = W_temp(j);
                j = j+1;
            end
        end
        
%         W_temp = W2;
%         W2 = zeros(N_tot,1);        
%         for i = 1:length(index_l)                
%             W2(index_l(i)) = W_temp(i);            
%         end

        ww = W2;
        ww(~ww) = inf;
        W_tot = [W_tot;W2'];
                
        
        display(['-------Convergence!!!--------']);        
        display(['Converged count = ',num2str(sum((abs(dataset_learn*W2)<max_y_threshold)))]);        
        display(['Cost in last iteration = ',num2str(cost_one_itr)]);
        display(['Norm-2 of W = ',num2str(norm(W))]);
        display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
        display(['Maximum projection = ',num2str(max(abs(dataset_learn*W2)))]);  
        display(['Minum of W = ',num2str(min(abs(ww)))]);  
        display('-----------------------------------');
        display(' ');   
        
%         fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
%         fprintf(fid, '%f \t',W');
%         fprintf(fid, '\n');
%         fclose(fid);
%     
% 
%         fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
%         fprintf(fid, '%d \t',learn_itr);
%         fprintf(fid, '\n');
%         fclose(fid);
% 
%         fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_cost_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
%         fprintf(fid, '%f \t',log(cost_per_itr_av));
%         fprintf(fid, '\n');
%         fclose(fid);
%     
%         fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r'); 
%         if (fid > -1)
%             current_count = fscanf(fid, '%d');
%             fclose(fid);
%         else
%             current_count = 0;
%             display('error!');
%         end
%     
%         fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'w'); 
%         fprintf(fid,'%d',max(current_count-1,0));
%         fclose(fid);
    end
    %======================================================================
    
end
111;

   