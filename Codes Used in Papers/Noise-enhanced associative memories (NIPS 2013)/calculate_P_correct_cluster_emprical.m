%==========================================================================
%***************FUNCTION: calculate_P_correct_cluster_emprical*************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% pattern_neur_noise: The amount of internal noise at pattern neuorns side
% const_neur_noise: The amount of internal noise at pattern neuorns side
% no_simulated_instances: The number of iterations to be conducted to calculate the rate of success
% gamma: The update threshold in the majority voting algorithm for pattern neurons
% index: The index of the simulation ensemble
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% success_rate: The number of times the algorithm has secceeded ine liminating external errors
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function reads the result of the learning phase for a clustered
% neural associative memory and determines the degree distribution of the
% network from a clustered point of view, i.e. the case that all the
% constraint nodes in a cluster are contracted to one node.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================






function [success_rate] = calculate_P_correct_cluster_emprical(N,K,L,alpha0,beta0,theta0,pattern_neur_noise,const_neur_noise,no_simulated_instances,gamma,index)
% function [success_rate] = calculate_P_correct_cluster_emprical(err_bits)
% N = 40;
% K = 20;
% L = 50;
% alpha0 = 0.95;
% beta0 = 0.75;
% theta0 = 0.05;
% 
% const_neur_noise = 0;
% no_simulated_instances = 1000;
% gamma = 0.95;
gamma_BFO = gamma;
gamma_BFS = gamma;
y_min = -10;
y_max = 100;
recall_algorithm_option = 1;
max_noise_amp = 1;
% index = 1;
% err_bits = 4;

% 
for err_bits = 1:4
    success_count = zeros(1,L);
err_bits
%==============================INITIALIZATION==============================

%--------------------------Simulation Parameters---------------------------   

const_update_threshold = const_neur_noise + 1e-3;       % The update threshold for the constraint neurons
%--------------------------------------------------------------------------                            

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),...
    '_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),...
    '_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);           
%--------------------------------------------------------------------------

%==========================================================================


%============================MAIN LOOP=====================================
for net_simul_itr = 1:no_simulated_instances                % Simulate the error correction procedure for the given ensemble.                                       
                   
    %-----------------------Generate the Pattern---------------------------
    mu = 1+floor((pattern_learn_number-1)*rand);            % Pick a pattern index at random               
    temp = dec2bin(mu_list(1,mu),KK);                  
    randindex = mu_list(2:KK+1,mu);    
    message = zeros(1,K_tot);                               % Generate the message from the index                        
           
    for j = 1:KK    
        message(randindex(j)) = (z_max-z_min)*(temp(j) - 48)+z_min;                                     
    end            
    pattern = message*G_tot;  
    %----------------------------------------------------------------------                                                                                                                                                           
        
    %------------------------Perform the Recall Step-----------------------                
    for l = 1:L
                       
        %-----------------------Read the Weight Matrix---------------------
        load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),...
            '_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),...
            '_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);
        
        n = length(index_l);                                % n is the number of pattern neurons in the cluster
        
        fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),...
            '_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),...
            '_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt'], 'r');
        
        if (fid == -1)            
            error('No weight matrix!');
        end
            
        W = fscanf(fid, '%f',[n,inf]);        
        W = W';              
        fclose(fid);            
        %------------------------------------------------------------------    
            
        %--------------------------Generate Noise--------------------------
        nois = zeros(1,n);                                  % This is the noise added to the whole pattern of length N*L_in                    
        pp = 1+floor((n-1)*rand(1,err_bits));                            
        for h = 1:err_bits                
            nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));                        
        end        
        %------------------------------------------------------------------
                                    
      
        %---------------Form the Subpatterns for Each Cluster--------------
        pattern_temp = zeros(1,length(index_l));
        for iji = 1:length(index_l)
            pattern_temp(iji) = pattern(index_l(iji));            
        end    
        x_temp = pattern_temp + nois;
        %------------------------------------------------------------------
                       
        %-----------------------Iterate Until Convergence------------------
        [x_out,~] = faulty_recall_step(W,x_temp,n,gamma_BFO,gamma_BFS,y_min,y_max,recall_algorithm_option,pattern_neur_noise,const_neur_noise,const_update_threshold);                                         
        %------------------------------------------------------------------
                        
        %-----------Verify the Success of the Recall Algorithm-------------
        if (norm(x_out-pattern_temp)<.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                                       
            success_count(l) = success_count(l) + 1;            
        end
        %------------------------------------------------------------------
        
    end
end
%==========================================================================


%======================STORE THE RESULTS===================================
success_rate = success_count/net_simul_itr/index;
save(['./Simulation_Results/P_Correct_Empirical/P_C_',num2str(err_bits),'_upsilon_',num2str(pattern_neur_noise),'_nu_',num2str(const_neur_noise),'_gamma_',num2str(gamma_BFO),'.mat'],'success_rate');
%==========================================================================    
end