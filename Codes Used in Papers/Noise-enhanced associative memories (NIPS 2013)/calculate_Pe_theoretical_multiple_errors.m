%==========================================================================
%******************FUNCTION: calculate_Pe_theoretical**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% error_bits_in: The number of external errors to the network
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% const_neur_noise: The maximum amount of noise a constraint neuron will "suffer" from
% no_errors: The maximum number of errors each cluster tries to correct
% gamma: The update threshold in the original bit-flipping recall algorithm 
% index: The index of the simulation setup among various random scenarios
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe_tot: The upper bound on BER of the inter-cluster error correcing algorithm
% peeling_thr: The threshold on the number of correctable external errors returned by the "generalized" peeling algorithm
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates a theoretical upper bound on the probability of
% error in a clustered neural networks built out of faulty neurons. The
% bound is explained in the technical report of our ISIT 2013 paper.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


% function [Pe_tot,peeling_thr] = calculate_Pe_theoretical_multiple_errors(N,K,L,alpha0,beta0,theta0,error_bits_in,pattern_neur_noise,const_neur_noise,gamma,no_errors,index)

%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_'...
    ,num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);            
%--------------------------------------------------------------------------

%==========================================================================


%======IF NOT CALCULATE BEFORE, CALCULATE THE BOUND FOR EACH CLUSTER=======
P_c = [];
for e = 1:no_errors    
    save_file = ['./Simulation_Results/P_Correct_Empirical/P_C_',num2str(e),'_upsilon_',num2str(pattern_neur_noise),'_nu_0_gamma_0.95.mat'];
    if (~exist(save_file,'file'))                        
        P_c = [P_c;calculate_P_correct_cluster_emprical(N,K,L,alpha0,beta0,theta0,pattern_neur_noise,const_neur_noise,no_simulated_instances,gamma,index)];
    else
        load(['./Simulation_Results/P_Correct_Empirical/P_C_',num2str(e),'_upsilon_',num2str(pattern_neur_noise),'_nu_0_gamma_0.95.mat']);
        P_c = [P_c;success_rate];
    end    
end
%==========================================================================


%===================CALCULATE THE GLOBAL BOUND ON BER======================
[~,~,lambda,rho] = read_cluster_degree(N,K,L,alpha0,beta0,theta0);
p_e_range = sort(error_bits_in)/N_tot;

% P_c = 5*P_c;
% Pe_tot = recurions_P_e_v3(P_c,lambda,rho,p_e_range);
Pe_tot = recurions_P_e(P_c,lambda,rho,p_e_range,4);
%==========================================================================


%=================CALCULATE THE THRESHOLD OF THE PEELING ALGO==============
e = 1;
load(['./Simulation_Results/P_Correct_Empirical/P_C_',num2str(e),'_vs_upsilon'])    
avg_edge = e*p_c_ave;
for e = 2:4
    load(['./Simulation_Results/P_Correct_Empirical/P_C_',num2str(e),'_vs_upsilon'])    
    avg_edge = avg_edge+e*p_c_ave; 
end

load(['./Simulation_Results/P_Correct_Empirical/P_C_',num2str(4),'_vs_upsilon'])    
for i = 1:length(upsilon_range)
    if (abs(upsilon_range(i) - pattern_neur_noise)<1e-5)
        avg_edge_1 = avg_edge(i);
        break
    end
end
[~,peeling_thr]=recurions_P_e_peeling(avg_edge_1,lambda,rho,sort(p_e_range));
%==========================================================================

