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

function clustered_neural_learn_v1(N,K,L,const_learn,cluster_index,alpha0,beta0,theta0,index)

%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);           
%--------------------------------------------------------------------------



%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
mkdir(['/scratch/amir/Clustered_Neural/Learn_Results'],['N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)]);        % Create a specific folder for the current N and K

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%-------------------Construct the Pattern Generator Matrix-----------------
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster_index),'.mat']);
G = [];
for i = 1:length(index_l)
    G = [G,G_tot(:,index_l(i))];              
end    
[~,n]= size(G);
%--------------------------------------------------------------------------

%==========================================================================
    

%%
%=============================LEARNING PHASE===============================
for itr = 1: const_learn

    %------------------Randomly Initialize the Weight Matrix---------------
%     W = random_vector(N,2*round(log(N)));
    W = rand(n,1);
%     W = W';
    %----------------------------------------------------------------------
    
    %---------------------Adjust Simulation Parameters---------------------
    theta = theta0;   
    theta_itr = 1;
    total_cost2 = [];        
    %----------------------------------------------------------------------

    %--------------------------Main Learning Loop--------------------------           
    for learn_itr = 1:learn_itr_max                         % Repeat the steps below several times in order to enhance the quality of the learning phase        
        temp_cost2 =0;
        temp_cost = 0;        
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
            
            if (norm(x)>0)
                alph = max(alpha0/norm(x)^2/learn_itr,0);            
            end
            if (norm(W) > .001)                                
                W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
            end            
            
            W = soft_threshold(W,0.0001);                       % Set entries less than 0.0001 to zero.
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
            temp_cost = temp_cost + y^2;
            temp_cost2 = temp_cost2+y^2;
            if (abs(y)> max_y)
                max_y = abs(y);
            end                                        
            %--------------------------------------------------------------
            
                        
            %-------------------------Update theta-------------------------
            if (mod(mu,20000)==0)
                theta_itr = theta_itr + 1;
                theta = theta0/theta_itr;
            end
            %--------------------------------------------------------------
            
            
            %-----------------------Display Progress-----------------------
            if (mod(mu,10000)==0)    
                total_cost2 = [total_cost2,temp_cost2];
                temp_cost2 = 0;                
                temp_cost
                norm(W)        
                sum(abs(sign(W)))   
                max_y
            end
            %--------------------------------------------------------------
        end 
         
        
            
        %----------------------Check Convergence---------------------------
        if (max_y < 0.005)            
            break;
        end
        %------------------------------------------------------------------

    end  
    %======================================================================


    %==========================SAVE THE RESULTS============================    
    if (max_y < 0.005)
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
        fprintf(fid, '%f \t',W');
        fprintf(fid, '\n');
        fclose(fid);
    end

    fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
    fprintf(fid, '%d \t',learn_itr);
    fprintf(fid, '\n');
    fclose(fid);

    fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_cost_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'a');
    fprintf(fid, '%f \t',log(total_cost2));
    fprintf(fid, '\n');
    fclose(fid);
    
    fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'r'); 
    if (fid > -1)
        current_count = fscanf(fid, '%d');
        fclose(fid);
    else
        display('error!');
    end
    
    fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster_index),'_index_',num2str(index),'.txt'], 'w'); 
    fprintf(fid,'%d',current_count-1);
    fclose(fid);
    %======================================================================
    
end

   