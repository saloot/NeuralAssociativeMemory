%==========================================================================
%********************FUNCTION: clustered_neural_learn**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% const_learn: Number of constraints which must be learned during the learning phase
% train_set_index: The index of the simulation setup among various random scenarios
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
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

function neural_learn_ITW_journal(N,K,const_learn,train_set_index,alpha0,beta0,theta0,mca_flag)

%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_parameters_v1_N_',num2str(N),'_K_',num2str(K),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
if (mca_flag == 0)
    max_y_thr = .001;
else
    max_y_thr = .01;
end


addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
mkdir(['/scratch/amir/ITW_Journal/Learn_Results'],['N_',num2str(N),'_K_',num2str(K)]);        % Create a specific folder for the current N and K

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%-------------------Adjust the Training-Related Materials------------------
load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(train_set_index),'.mat']);   
%--------------------------------------------------------------------------

%-------------Update the Number of Constraints in the Process--------------
if (mca_flag == 0)
    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'r'); 
       
    if (fid > -1)
        current_count = max(fscanf(fid, '%d'),0);            
        fclose(fid);    
    else    
        current_count = 0;            
    end        

    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'r');            
    if (fid > -1)           
        W = fscanf(fid, '%f',[N,N-K]);                
        W = W';            
        fclose(fid);                                    
        [m,~] = size(W);            
    else    
        m = 0;    
    end
    
    const_to_learn = min(const_learn,N-K-m-current_count);

    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'w'); 
    fprintf(fid,'%d',current_count+const_to_learn);
    fclose(fid);
    current_count = current_count + const_to_learn;
else
    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'r'); 
       
    if (fid > -1)
        current_count = max(fscanf(fid, '%d'),0);            
        fclose(fid);    
    else    
        current_count = 0;            
    end        

    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'r');            
    if (fid > -1)           
        W = fscanf(fid, '%f',[N,N-K]);                
        W = W';            
        fclose(fid);                                    
        [m,~] = size(W);            
    else    
        m = 0;    
    end
    
    const_to_learn = min(const_learn,N-K-m-current_count);

    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'w'); 
    fprintf(fid,'%d',current_count+const_to_learn);
    fclose(fid);
    current_count = current_count + const_to_learn;    
end
%--------------------------------------------------------------------------

%==========================================================================




%=============================LEARNING PHASE===============================
for itr = 1: const_to_learn

    %------------------Randomly Initialize the Weight Matrix---------------
    W = random_vector(N,4*round(log(N)));
    W = W';
%     W = rand(N,1);
%     W = W/norm(W);   
    %----------------------------------------------------------------------
    
    %--------------------Adjust Simulation Parameters----------------------
    theta = theta0;       
    total_cost2 = [];
    theta_itr = 1;    
    %----------------------------------------------------------------------

    %--------------------------Main Learning Loop------------------------------                
    for learn_itr = 1:learn_itr_max         % Repeat the steps below several times in order to enhance the quality of the learning phase        
        temp_cost2 =0;
        temp_cost = 0;
        temp_cost3 = 0;
        max_y = 0;
        
        alph0 = max(50*alpha0/(50+(learn_itr)),0.0005);
        continue_flag = 0;
        for mu = 1:length(mu_list)
            
            %------------------Pick a Pattern at Random--------------------
            index = 1+floor((length(mu_list)-1)*rand);
            temp = dec2bin(mu_list(index),K);                      
            message = zeros(1,K);           % Generate the message from the index                                             
            for j = 1:K                    
                message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;                             
                
            end            
            x = message*G;                                                    
            %--------------------------------------------------------------
            
            %-----------------------------Generate Noise---------------------------
            if (mca_flag > 0)
                nois = zeros(1,N);          % This is the noise added to the whole pattern of length N*L_in        
                pp = rand;
                if (pp<.001)
                    pp = 1+floor((N-1)*rand(1,mca_flag));                
                    for h = 1:mca_flag        
                        nois(pp(h)) =1;%(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
                    end     
                end
                x = x+nois;
            end
            %----------------------------------------------------------------------                            
            
            
            %--------------------Update the Weight Vector------------------
            y = x*W;                 % Calculate the projection of x on the vector W(j,:).
            
            if (norm(x)>.001)                
                alph = alph0/norm(x)^2;
            end
            if (norm(W) > .001)
                W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
                % W = soft_threshold(W - alph*y*( x' - (y*W/(norm(W))^2)),theta);%-beta0*W.*(1-(tanh(100*W.^2)).^2);
                
                % W = soft_threshold(W - alph*y*(x'*(2*norm(W)^2-1) - y*W),theta);
                % W = soft_threshold(W - alph*y*(x'-y*W)+ alph*(1-norm(W)^2)*W,theta);
                % W = soft_threshold(W - alph*y*( x' - (y*W/(norm(W))^2)),theta);%-beta0*W.*(1-(tanh(100*W(j,:).^2)).^2);
            end
            W = soft_threshold(W,0.0001);                       % Set entries less than 0.0001 to zero.
            %--------------------------------------------------------------
        
            
            %-----------------Check for Numerical Errors-------------------
            if (sum(isnan(W))>0)    
                display('Error Nan');        
                break;                                    
            end
            
            if (norm(W) < .001)                
                display('Error zero!');
            end
            %--------------------------------------------------------------

            %-----------------Update the Simulation Costs------------------
            temp_cost = temp_cost + y^2;
            temp_cost2 = temp_cost2+y^2;
            temp_cost3 = temp_cost3 + y^2;
            if (abs(y)> max_y)
                max_y = abs(y);
            end                                        
            %--------------------------------------------------------------                                                            
            
            %-----------------------Display Progress-----------------------                        
            if (mod(mu,200)==0)
                total_cost2 = [total_cost2,temp_cost2];
                temp_cost2 = 0;
            end
            if (mod(mu,10000)==0)                    
                if (temp_cost3<1e-3)
                    continue_flag = 1;                   
                else                    
                    temp_cost3 = 0;
                end
                norm(G*W)
                temp_cost
                norm(W)        
                sum(abs(sign(W)))                                
            end
            if (mod(mu,20000)==0)
                111;
            end
            %--------------------------------------------------------------
          
            if (continue_flag == 1)
                break;
            end
            %-------------------------Update theta-------------------------
            if (mod(mu,1000)==0)
                theta_itr = theta_itr + 1;
                theta = theta0/theta_itr;
            end
            %--------------------------------------------------------------
        end 
        
    
%         theta = max(theta0/(1+learn_itr),0);       
    
        %-------------------Update Learnig Cost--------------------------------    
            
        %----------------------Check Convergence---------------------------
        if (mca_flag == 0)
            if (max_y < max_y_thr)        
                break;
            end
        else
            if (temp_cost/mu<.1*max_y_thr)
                break;
            end
        end
        %------------------------------------------------------------------      
        
         display('Another learn_itr has passed!');
         
    end  
    
%     display('Another learn_itr has passed!');
    %======================================================================


    %========================SAVE THE RESULTS==============================    
    if (mca_flag ==0)
        if (max_y <max_y_thr)        
            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'a');
            fprintf(fid, '%f \t',W');
            fprintf(fid, '\n');
            fclose(fid);
    

            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'a');
            fprintf(fid, '%d \t',learn_itr);
            fprintf(fid, '\n');
            fclose(fid);

            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Learn_cost_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'a');
            fprintf(fid, '%f \t',log(total_cost2));
            fprintf(fid, '\n');
            fclose(fid);
    
            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'r'); 
            if (fid > -1)
                current_count = fscanf(fid, '%d');
                fclose(fid);
            else
                display('error!');
            end
    
            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'w'); 
            fprintf(fid,'%d',current_count-1);
            fclose(fid);
        end
    else
        if (temp_cost/mu<.1*max_y_thr)
            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'a');
            fprintf(fid, '%f \t',W');
            fprintf(fid, '\n');
            fclose(fid);
    

            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Learn_itr_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'a');
            fprintf(fid, '%d \t',learn_itr);
            fprintf(fid, '\n');
            fclose(fid);

            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Learn_cost_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'a');
            fprintf(fid, '%f \t',log(total_cost2));
            fprintf(fid, '\n');
            fclose(fid);
    
            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'r'); 
            if (fid > -1)
                current_count = fscanf(fid, '%d');
                fclose(fid);
            else
                display('error!');
            end
    
            fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/learn_count_mca_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(train_set_index),'.txt'], 'w'); 
            fprintf(fid,'%d',current_count-1);
            fclose(fid);
        end
    end
    %==========================================================================
end

