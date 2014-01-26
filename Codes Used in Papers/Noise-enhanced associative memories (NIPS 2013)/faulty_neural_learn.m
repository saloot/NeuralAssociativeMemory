%==========================================================================
%**********************FUNCTION: faulty_neural_learn***********************
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

function faulty_neural_learn(N,K,L,const_learn,cluster_index,alpha0,beta0,theta0,index)

%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
mkdir(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results'],['N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)]);        % Create a specific folder for the current N and K

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%-------------------Construct the Pattern Generator Matrix-----------------
load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster_index),'.mat']);
G = [];
for i = 1:length(index_l)
    G = [G,G_tot(:,index_l(i))];              
end    
[~,n]= size(G);
%--------------------------------------------------------------------------

%-------------Update the Number of Constraints in the Process--------------
fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r'); 
if (fid > -1)
    current_count = max(fscanf(fid, '%d'),0);            
    fclose(fid);    
else    
    current_count = 0;            
end        

fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r');
if (fid > -1)           
    W = fscanf(fid, '%f',[n,inf]);                
    W = W';            
    fclose(fid);                                    
    [m,~] = size(W);            
else    
    m = 0;    
end
    
N_const_max = N-K;
d_max = 7;
const_to_learn = min(const_learn,N_const_max-m-current_count);
% const_to_learn = min(const_learn,N-(K/2)-m);


fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'w'); 
fprintf(fid,'%d',current_count+const_to_learn);
fclose(fid);
%--------------------------------------------------------------------------


min_W = 0.02;
max_y_threshold = 0.01;
d_min_row = ceil(N_const_max/3);
%==========================================================================
    

%%
%=============================LEARNING PHASE===============================
for itr = 1: const_learn

    %-----------------------Read Current Learned Matrix--------------------
    fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r');
    if (fid > -1)           
        W_tot = fscanf(fid, '%f',[n,inf]);                
        W_tot = W_tot';            
        fclose(fid);                                    
        [m,~] = size(W_tot);            
    else    
        m = 0;    
    end
    
    if (m>4)
        index_w = zeros(1,n);
        s = sum(abs(sign(W_tot)));
        for i = 1:n
            pp = rand;
            if (s(i)/d_max <= pp)
                index_w(i) = 1;
            end
        end
    else
        index_w = ones(1,n);
    end
    
    n_new = sum(index_w);
    
    s = sum(abs(sign(W_tot)))/d_max;
    if (n_new <=d_min_row)
        [ss,ind] = sort(s,2);
        for i = n_new+1:d_min_row
            index_w(ind(i)) = 1;
        end
    end
    
    n_new = sum(index_w);
            
    if (n_new/n < .3)
        theta1 = 0;
        alpha1 = alpha0/3;        
    else
        theta1 = theta0*n_new/n;        
    end    
    %----------------------------------------------------------------------    
    
    %------------------Randomly Initialize the Weight Matrix---------------
%     W = random_vector(N,2*round(log(N)));
    W = random_vector(n_new,2*round(log(n_new)));
    W = W';
    W = rand(n_new,1);
    W = W/norm(W);
%     W = W';
    %----------------------------------------------------------------------
    
    %---------------------Adjust Simulation Parameters---------------------
    theta = theta0;   
    theta_itr = 1;   
    cost_per_itr_av = [];            
    alph_itr = 0;    
    alpha1 = alpha0;
    verify_flag = 0;
    %----------------------------------------------------------------------

    %--------------------------Main Learning Loop--------------------------           
    for learn_itr = 1:learn_itr_max                         % Repeat the steps below several times in order to enhance the quality of the learning phase        
        
        max_y = 0;
        cost_one_itr = 0;        
        cost_short_term =0;
        if (mod(learn_itr,200) == 1)            
            alph_itr = alph_itr + 1;
            alph0 = max(alpha1/alph_itr,0.0005);        
        end
        
        for ii = 1:pattern_learn_number
            
            %------------------Pick a Pattern at Random--------------------
            mu = 1+floor((pattern_learn_number-1)*rand); 
            temp = dec2bin(mu_list(1,mu),KK);      
            randindex = mu_list(2:KK+1,mu);
            message = zeros(1,K_tot);                       % Generate the message from the index                        

            for j = 1:KK
                message(randindex(j)) = (z_max-z_min)*(temp(j) - 48)+z_min;                             
            end    
            
            x_temp = message*G;
            x = [];
            for i = 1:n
                if (index_w(i) > 0)
                    x = [x,x_temp(i)];
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
%                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
            
%                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,0.05);W.*(1-(tanh(100*W(j,:).^2)).^2)
%                  W = soft_threshold(W - alph*y*( x' - (y*W/(norm(W))^2)),theta0);%-beta0*W.*(1-(tanh(100*W.^2)).^2);
                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
%                 W = soft_threshold(W - alph*y*(x'-y*W)+ alph*(1-norm(W)^2)*W,theta); 
            end            
            
            W = soft_threshold(W,0.001);                       % Set entries less than 0.0001 to zero.
            %--------------------------------------------------------------
        
            
             %------------------Check for Numerical Errors------------------
            if (sum(isnan(W))>0)    
                display('Error Nan');        
                break;                                    
            end
            
            if (norm(W) < .001)                
                display('Error zero!');
            end
            %--------------------------------------------------------------

            %-----------------Update the Simulation Costs------------------
            cost_one_itr = cost_one_itr + abs(x*W);
            cost_short_term = cost_short_term+norm(x*W);
%             temp_cost = temp_cost + y^2;
%             temp_cost2 = temp_cost2+y^2;
            
            if (abs(x*W)>max_y)
                max_y = abs(x*W);
            end                                    
            %--------------------------------------------------------------
                                                
        end 
         
        
            
        %-------------------------Update theta-------------------------            
        if (mod(learn_itr,3)==0)                
            theta_itr = theta_itr + 1;                
            theta = theta1/theta_itr;            
        end        
        %--------------------------------------------------------------
        
        %--------------------Update the Simulation Costs-------------------
        cost_per_itr_av = [cost_per_itr_av,cost_short_term];
        cost_short_term = 0;
        ww = W;
        ww(~ww) = inf;
        W=soft_threshold(W,min_W);         
               
        ww = W;
        ww(~ww) = inf;                                
        %------------------------------------------------------------------                
         
        
        %--------------------------Display Progress------------------------
        if (mod(learn_itr,1)==0)                        
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);            
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);
            display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
            display(['Maximum projection = ',num2str(max_y)]);  
            display(['Minum of W = ',num2str(min(abs(ww)))]);  
            display('-----------------------------------');
            111;
            display(' ');             
        end           
        %--------------------------------------------------------------
        
            
        %----------------------Check Convergence---------------------------
        if (max_y < max_y_threshold)    
            if (verify_flag == 1)
                break;
            else
                verify_flag = 1;
            end
        else
            verify_flag = 0;
        end
        %------------------------------------------------------------------
                
        if (mod(learn_itr,5)==0)                        
            111;
        end

    end  
    %======================================================================


    %==========================SAVE THE RESULTS============================    
    W_temp = W;        
    W2 = zeros(N,1);
    j = 1;        
    for i = 1:N
        if (index_w(i) > 0)               
            W2(i) = W_temp(j);                
            j = j+1;            
        end        
    end
    W = W2;
        
    %------Verify one More Time-----------
    max_y = 0;
    for ii = 1:pattern_learn_number
                        
        %------------------Pick a Pattern at Random--------------------
        mu = 1+floor((pattern_learn_number-1)*rand);             
        temp = dec2bin(mu_list(1,mu),KK);                  
        randindex = mu_list(2:KK+1,mu);            
        message = zeros(1,K_tot);                       % Generate the message from the index                        

        for j = 1:KK
            message(randindex(j)) = (z_max-z_min)*(temp(j) - 48)+z_min;                             
        end    

        x = message*G;                        
        %--------------------------------------------------------------
        
        %---------------------Update the Weight Vector-----------------
        y = x*W;                 % Calculate the projection of x on the vector W(j,:).
            
            
        if (abs(x*W)>max_y)
            max_y = abs(x*W);
        end                                    
        %--------------------------------------------------------------
            
                                                            
    end 
    %-------------------------------------
    
    if (max_y < 0.005)
        fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
        fprintf(fid, '%f \t',W');
        fprintf(fid, '\n');
        fclose(fid);
    end

    fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
    fprintf(fid, '%d \t',learn_itr);
    fprintf(fid, '\n');
    fclose(fid);

    fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_cost_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
    fprintf(fid, '%f \t',log(cost_per_itr_av));
    fprintf(fid, '\n');
    fclose(fid);
    
    fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r'); 
    if (fid > -1)
        current_count = fscanf(fid, '%d');
        fclose(fid);
    else
        display('error!');
    end
    
    fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'w'); 
    fprintf(fid,'%d',current_count-1);
    fclose(fid);
    %======================================================================
    
end

   