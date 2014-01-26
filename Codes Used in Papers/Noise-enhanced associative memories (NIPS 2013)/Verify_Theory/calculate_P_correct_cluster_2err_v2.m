%==========================================================================
%**************FUNCTION: calculate_P_correct_cluster_2err_v2***************
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
% index: The index of the simulation ensemble
% no_try: The number of times we try to correct external errors
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pc_tot: The probability that EACH cluster correct one external error
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function reads the result of the learning phase for a clustered
% neural associative memory and calculates an upper bound on the
% probability that each cluster could correct two external error.

% NOTE: The different between this code and the earlier version
% (calculate_P_correct_cluster) is the fact that this version is more time
% consuming as it avergaes over all combinations. Hence, it should be more
% accurate. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

% N = 40;
% K = 20;
% L = 50;
% alpha0 = 0.95;
% beta0 = 0.75;
% theta0 = 0.05;
% pattern_neur_noise = 0.5;
% const_neur_noise = 0;
% gamma = .95;
% index = 1;
% no_try = 20;

function [Pc_tot,P2_1_tot,P2_2_tot] = calculate_P_correct_cluster_2err_v2(N,K,L,alpha0,beta0,theta0,pattern_neur_noise,const_neur_noise,gamma,index,no_try)

%%   
%==============================INITIALIZATION==============================
Pe_tot = zeros(1,L);                            % The probability that each cluster could NOT correct two errors
P2_1_tot = zeros(1,L);                          % The (average) probability that in each cluster, the non-corrupted pattern neurons gets corrupted
P2_2_tot = zeros(1,L);                          % The (average) probability that in each cluster, the corrupted pattern neurons do not get corrected

%--------------------Identify Network Parameters-----------------------
[min_W_distribution,mean_W_distributiond_max,~,~,~,~,no_const_distribution] = find_network_characteristics(N,K,L,alpha0,beta0,theta0,index);
                                                                              
%----------------------------------------------------------------------

%==========================================================================


%%
%==============================MAIN LOOP===================================
for l = 1:L
        
    %--------------------Within Loop Initliazations------------------------
    load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),...
        '_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),...
        '_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);                
    n = length(index_l);                        % Determine the number of pattern nodes in the corresponding cluster                
    %----------------------------------------------------------------------
               
    %---------Determine the Degree Distribution of Current Cluster---------
    [~,deg_dist_column,lambda,~,~,~,~,~,~,~,~,~,~,~] = read_cluster_degree_local_v2(N,K,L,alpha0,beta0,theta0,l);        
    %----------------------------------------------------------------------                       
    
    
    %-----------------Average Over the First Corrupted Neuron--------------
    temp_d1 = 0;  
    for ii = 1:length(deg_dist_column)                    
        if (lambda(ii)==0)        
            continue;            
        end        
        d1 = deg_dist_column(ii);                           % d1 is the first noisy node            
        P_common_neighbor_noisy = d1/no_const_distribution(l);
                        
        %-------------Average Over the Second Corrupted Neuron-------------
        temp_d2 = 0;        
        for ij = 1:length(deg_dist_column)                        
            if ( (lambda(ij)==0)  )                    
                continue;                
            end            
            d2 = deg_dist_column(ij);                       % d2 is the second noisy node
            
            %-------------Average Over the Non-Corrupted Neurons-----------
            temp_d3 = 0;  
            for ik = 1:length(deg_dist_column)                                            
                if ( (lambda(ik)==0) )                        
                    continue;                                            
                end                
                d3 = deg_dist_column(ik);                    % d3 is the non-noisy node
                                    
                %---------------------Average Over d_bar-------------------
                P2_2_1 = 0;                    
                P2_2_2 = 0;
                    
                P2_1 = 0;
                temp_4 = 0;
                P_cor_dbar = [];
                for d_bar = 0:min(d1,d2)
                    P_neighbor = ((P_common_neighbor_noisy)^d_bar)* ((1-P_common_neighbor_noisy)^(min(d1,d2)-d_bar));                            
                    [temp_dbar2_1_no_update,temp_dbar2_1_wrong_dir,temp_dbar2_2_no_update,temp_dbar2_2_wrong_dir] = P2_2_calculation_e2_v2(gamma,pattern_neur_noise,d1,d2,0,d_bar);                                                        
                    
                    P_e_get_corrected1 = 1-(temp_dbar2_1_no_update+temp_dbar2_1_wrong_dir);
                    P_e_get_corrected2 = 1-(temp_dbar2_2_no_update+temp_dbar2_2_wrong_dir);
                            
                    Already_calculated_flag(ii,ij) = temp_dbar;                                                
                    P2_2_1 = P2_2_1 + nchoosek(min(d1,d2),d_bar) * P_neighbor * temp_dbar2_1_no_update;
                    P2_2_2 = P2_2_2 + nchoosek(min(d1,d2),d_bar) * P_neighbor * temp_dbar2_2_no_update;
                    
                    P_common_neighbor_non_noisy = (d1+d2-d_bar)/no_const_distribution(l);                        
                  
                    temp_dbar1 = P2_1_calculation(gamma,pattern_neur_noise,d3,0,min(P_common_neighbor_non_noisy,1));
                    P2_1 = P2_1 + nchoosek(min(d1,d2),d_bar) * P_neighbor * temp_dbar1;                                                                                                                       
                        
                    P_c1_temp = ((1-temp_dbar1)^(n-2))*P_e_get_corrected1*P_e_get_corrected2 + ...
                                ((1-temp_dbar1)^(n-2))*(1-P_e_get_corrected1)*P_e_get_corrected2 + ...
                                ((1-temp_dbar1)^(n-2))*P_e_get_corrected1*(1-P_e_get_corrected2) + ...
                                (n-2)*temp_dbar1 * ((1-temp_dbar1)^(n-3))*P_e_get_corrected1*P_e_get_corrected2;
                            
                    q = (1-temp_dbar2_1_no_update)*(1-temp_dbar2_2_no_update)*((1-temp_dbar1)^(n-2)) + ...
                        P_e_get_corrected1 * (1-temp_dbar2_2_no_update)*(n-2)*temp_dbar1*((1-temp_dbar1)^(n-3)) + ...
                        P_e_get_corrected2 * (1-temp_dbar2_1_no_update)*(n-2)*temp_dbar1*((1-temp_dbar1)^(n-3)) + ...
                        P_e_get_corrected1 * P_e_get_corrected2*(n-2)*(n-3)*.5*(temp_dbar1)^2*((1-temp_dbar1)^(n-4));

                    P_cor_temp = P_c1_temp * (1-q^no_try)/(1-q);
                            
                    P_cor_dbar = [P_cor_dbar,nchoosek(min(d1,d2),d_bar) * P_neighbor*P_cor_temp];
                    temp_4 = temp_4 + nchoosek(min(d1,d2),d_bar) * P_neighbor*(1-P_cor_temp);
                end                    
                %----------------------------------------------------------                    
            
                P_cor_temp = ((1-P2_1)^(n-2))*(1-P2_2_1)^2;                   
                temp_d3 = temp_d3 + lambda(ik) *temp_4;          
                
            end
            %--------------------------------------------------------------                
            
            temp_d3            
            temp_d2 = temp_d2 + lambda(ij) *temp_d3;                      
        end
        %------------------------------------------------------------------
            
        temp_d2
        temp_d1 = temp_d1 + lambda(ii) *temp_d2;                  
    end
    Pe_tot(l) = temp_d1
    %----------------------------------------------------------------------                               
end

Pe_tot = Pe_tot.^no_try;
Pc_tot = 1-Pe_tot;
%==========================================================================