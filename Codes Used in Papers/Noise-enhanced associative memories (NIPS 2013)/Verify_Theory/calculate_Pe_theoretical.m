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
% gamma: The update threshold in the original bit-flipping recall algorithm 
% index: The index of the simulation setup among various random scenarios
% no_of_try: The number of times the algorithms tries to eliminate external errors
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% f: The upper bound on BER of the inter-cluster error correcing algorithm
% Pe_tor: The upper bound on the probability of mistake for each clsuter
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates a theoretical upper bound on the probability of
% error in a clustered neural networks built out of faulty neurons. The
% bound is explained in our ISIT 2013 paper.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function [f,Pe_tot] = calculate_Pe_theoretical(N,K,L,alpha0,beta0,theta0,error_bits_in,pattern_neur_noise,const_neur_noise,gamma,index,no_of_try)


%=============================INITIALIZATION===============================

%--------------------------Set Default Values------------------------------
if ~exist('error_bits_in','var')
    error_bits_in = [0:10:50];
end

if ~exist('no_errors','var')
    no_errors = 1;
end

if ~exist('N_in','var')
    N_in = 40;
end

if ~exist('K_in','var')
    K_in = 20;
end

if ~exist('L_in','var')
    L_in = 50;
end

if ~exist('alpha0','var')
    alpha0 = 0.95;
end

if ~exist('beta0','var')
    beta0 = 0.75;
end

if ~exist('theta0','var')
    theta0 = 0.05;
end

if ~exist('gamma','var')
    gamma = 0.95;
end

if ~exist('index','var')
    index = 1;
end

if ~exist('recall_algorithm_option','var')
    recall_algorithm_option = 1;
end

if ~exist('pattern_neur_noise_range','var')
    pattern_neur_noise_range = 0;
end

if ~exist('const_neur_noise_range','var')
    const_neur_noise_range = 0;
end
%--------------------------------------------------------------------------

%---------------------------Simulation Parameters--------------------------
constr_update_thereshold = const_neur_noise + 1e-3;                         % The update threshold for the constraint neurons    
Pai2_1_distribution = zeros(1,L);                                           % The probability that a constrain neuron sends an incorrect message for each cluster    
Pe_tot = zeros(1,L);                                                        % Total probability of error for each cluster
saved_Pe_fil = ['./Simulation_Results/Global_Upper_Bound/Pe_cluster_new_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_epsilon_',num2str(pattern_neur_noise),'_delta_',num2str(const_neur_noise),'.mat'];
%--------------------------------------------------------------------------

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_'...
    ,num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);            
%--------------------------------------------------------------------------

%==========================================================================


%======IF NOT CALCULATE BEFORE, CALCULATE THE BOUND FOR EACH CLUSTER=======
if (~exist(['./Simulation_Results/Global_Upper_Bound/Pe_cluster_loose2_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_epsilon_',num2str(pattern_neur_noise),'_delta_',num2str(const_neur_noise),'.mat'],'file'))
    Pe_tot = 1-calculate_P_correct_cluster(N,K,L,alpha0,beta0,theta0,pattern_neur_noise,const_neur_noise,gamma,no_errors,index);
    save(saved_Pe_fil,'Pe_tot');
else
    load(saved_Pe_fil);
end
%==========================================================================


%=====================CALCULATE THE GLOBAL BOUNT ON BER====================
% Pe_bound = mean(Pe_tot);
Pe_bound = max(Pe_tot);
    
[~,~,lambda,rho] = read_cluster_degree(N,K,L,alpha0,beta0,theta0);
p_e_range = sort(error_bits_in/N_tot);
f = recurions_P_e(1-Pe_bound,lambda,rho,p_e_range,no_of_try);  
%==========================================================================

%==========================SAVE THE RESULTS================================
clear error_bits_in
clear delta

%==========================================================================
