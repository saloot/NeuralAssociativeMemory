%==========================================================================
%*****************FUNCTION: calculate_P_correct_cluster********************
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
% gamma: The update threshold in the majority voting algorithm for pattern neurons
% no_errors: The number of external errors in each cluster
% index: The index of the simulation ensemble
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pc_tot: The probability that EACH cluster correct one external error
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function reads the result of the learning phase for a clustered
% neural associative memory and calculates an upper bound on the
% probability that each cluster could correct no_errors external error.
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================

function Pc_tot = calculate_P_correct_cluster(N,K,L,alpha0,beta0,theta0,pattern_neur_noise,const_neur_noise,gamma,no_errors,index)
    constr_update_thereshold = const_neur_noise + 1e-3;
    %--------------------Identify Network Parameters-----------------------
    [min_W_distribution,mean_W_distributiond_max,~,~,~,~,no_const_distribution] = find_network_characteristics(N,K,L,alpha0,beta0,theta0,index);
    %----------------------------------------------------------------------
    
    for l = 1:L
        
        %------------------Within Loop Initliazations----------------------
        if (const_neur_noise>0)
            Pai2_1 = max(0,(const_neur_noise-max(0,min_W_distribution(l)-constr_update_thereshold))/(2*const_neur_noise));
        else
            Pai2_1 = 0;
        end
        Pai2_1_distribution(l) = Pai2_1;
        
        load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),...
            '_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),...
            '_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);        
        
        n = length(index_l);                        % Determine the number of pattern nodes in the corresponding cluster
        
        temp2 = 0;
        %------------------------------------------------------------------
        
        %-------Determine the Degree Distribution of Current Cluster-------
        [~,deg_dist_column,lambda,~,~,~,~,~,~,~,~,~,~,~] = read_cluster_degree_local_v2(N,K,L,alpha0,beta0,theta0,l);
        %------------------------------------------------------------------
        
        %---------Average Over the Degree of the Corrupted Neuron----------
        for ii = 1:length(deg_dist_column)        
            if (lambda(ii)==0)
                continue;
            end
            
            %-------------------In-Loop Initializations--------------------
            d1 = deg_dist_column(ii);                                       % d1 is the degree of the noisy node
            P_neighbor = d1/no_const_distribution(l);                       % The probability of having a common neighbor with the noisy node
            temp = 0;                                                            
            %--------------------------------------------------------------
            
            %--Find Part of the Probabality of Correcting Two Ext. Errors--                    
            if (no_errors == 2)                        
                d_max = max(deg_dist_column.*sign(lambda));                                    
                d_mean = (deg_dist_column*lambda');                        
                P_common_neighbor = d_max/no_const_distribution(l);                        
                P_common_neighbor = d_mean/no_const_distribution(l);                        
                P2_2 = P2_2_calculation_e2(gamma,pattern_neur_noise,d1,0,P_common_neighbor);                                        
            end            
            %--------------------------------------------------------------
            
            %-----Average Over the Degree of the Non-Corrupted Neuron------
            for i = 1:length(deg_dist_column)        
                if (lambda(i) > 0)
                    d = deg_dist_column(i);
                    
                    %---Find the Probabality of Correcting One Ext. Error--
                    if (no_errors == 1)
                        P2_2 = P2_2_calculation(gamma,pattern_neur_noise,d1,Pai2_1);            % P2_2 is the probability that the corrupted pattern neuron is NOT corrected
                        P2_1 = P2_1_calculation(gamma,pattern_neur_noise,d,Pai2_1,P_neighbor);  % P2_1 is the (average) probability that a non-corrupted pattern neuron getsS corrected
                    %------------------------------------------------------
                    
                    %--Find the Probabality of Correcting Two Ext. Errors--
                    elseif (no_errors == 2)                        
                        P_common_neighbor = d_max/no_const_distribution(l);
                        P_common_neighbor = d_mean/no_const_distribution(l);                        
                        P2_1 = P2_1_calculation(gamma,pattern_neur_noise,d,0,min(2*P_common_neighbor,1));                        
                    %------------------------------------------------------
                    
                    
                    %--Find the Probabality of Correcting More Ext. Errors-
                    else
                        P2_1 = 1;
                        P2_2 = 1;
                    end
                    %------------------------------------------------------
                                        
                    %---Find the Prbability of Eliminating An Ext. Error---
                    P_cor_temp = ((1-P2_1)^(n-1))*(1-P2_2);
                    temp = temp + lambda(i) *(1-P_cor_temp);
                    %------------------------------------------------------
            
                end
            end
            %--------------------------------------------------------------
            
            %---------------Calculate the Weighted Average-----------------
            temp2 = temp2 + lambda(ii) * temp;
            %--------------------------------------------------------------
        end
        Pc_tot(l) = 1-temp2;
        %------------------------------------------------------------------
        
    end 