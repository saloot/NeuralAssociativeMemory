%==========================================================================
%******************FUNCTION: read_journal_results**************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N_in: The number of pattern nodes in the graph
% K_in: The dimension of the subspace of the pattern nodes
% L_in: The number of clusters if we have multiple levels (L = 1 for single level)
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% processed_error_bits: The list of processed number of initial erroneous nodes
% processed_PER_WTA: The list of processed Pattern Error Rates for the Winner-Take-All algorithm
% processed_PER_BFO: The list of processed Pattern Error Rates for the original bit-flipping algorithm
% processed_PER_BFS: The list of processed Pattern Error Rates for the simplified bit-flipping algorithm
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a clustered neural associative
% memory and reads the result of recall phase from the appropriate files. 
% The results will then be plotted and compared with theoretical values. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

%%
%=============================INITIALIZATION===============================

%-------------------------Simulation Variables-----------------------------
if (~exist('initialization_done','var'))    % If already not initialized by the GUI..
    K = 150;                            % Number of message bits
    N = 300;                            % Number of pattern neurons in the network
    Q = 8;                              % Number of quantization levels
    const_to_learn = 150;               % Number of contraints to learn over the patterns
    random_dataset_flag = 1;            % If 1, this flag tells the code to use the dataset generated by the file "neural_initialization.m". If 0, it will read the dataset from the file specified by the user
    index_in = 1;                       % Index of the random graph in the considered ensemble (for random_dataset_flag =1)
    no_simulated_instances = 1000;      % The number of patterns that are going to be denoised during the recall phase
    max_noise_amp = 1;                  % Maximum value of integer-valued noise added to each bit
    err_bits_range = [0:10];            % The number of bits that will be corrupted initially for the recall phase
    gamma_BFO = 0.95;                   % The update threshold for the Original Bit flipping algorithm
    gamma_BFS = 0.95;                   % The update threshold for the Simplified Bit flipping algorithm
    theta0 = 0.02;                      % The initial sparisty threshold
    alpha0 = 0.9;                       % The initial learning rate
    beta0 = 0.8;                        % The sparsity penalty
    nu = 0.025;                         % update threshold for the constraint neurons during the recall phase
    index_max = 1;                      % This is the maximum number of random scenarios generated for simulation
end
no_of_patterns_range = [50,100,200];
addpath(genpath('../Other_Works'));                          % Include the library of common functions in the search path
%--------------------------------------------------------------------------

%%
%==============================INITIALIZATION==============================
figure_legend = [];                         % This variable will be field with the legend for the graphs
hist_x_axis_itr = [0:20:1000];             % This is bin places for total learning iterations histogram
hist_x_axis_sparsity = [0:.05:1];           % This is bin places for total sparsity percentage

hist_x_itr = [];                            % This is output bin places produced BY the hist command for total learning iterations (to be used in bar plots)
hist_x_sparsity = [];                       % This is output bin places produced BY the hist command for sparsity percentage (to be used in bar plots)

hist_out_itr = [];                          % This is output values produced bY the hist command for total learning iterations (to be used in bar plots)
hist_out_sparsity = [];                     % This is output values produced bY the hist command for sparsity percentage (to be used in bar plots)

initialization_done = 1                     % Tell the other codes that the initialization has been done
%==========================================================================

%%
%============================PROCESS THE RESULTS===========================

%----------------If Not Already Read the Resultsof Our Work----------------
if (~exist('processed_BER_BFS','var'))
    plot_flag = 0
    run read_results
end
%--------------------------------------------------------------------------
    
%---------------------Execute Jankwoski's Algorithm------------------------
run Complex_Neur_Jankwoski
%--------------------------------------------------------------------------

%------------------------Execute Lee's Algorithm---------------------------
run Complex_Neur_Lee
%--------------------------------------------------------------------------

%==========================================================================

%==============================PLOT THE RESULTS============================
plot(sort(processed_error_bits_BFO),sort(processed_BER_BFO),'b-','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);    
hold on
plot(sort(processed_error_bits_Jankowski_non_random),sort(processed_BER_Jankowski_non_random),'b-.','LineWidth',2,'Color','red');
plot(sort(processed_error_bits_Jankowski_random),sort(processed_BER_Jankowski_random),'b-*','LineWidth',2,'Color','blue');

set(l,'Interpreter','latex')           
xlhand = get(gca,'xlabel');            
ylhand = get(gca,'ylabel');            
set(xlhand,'string','$\epsilon$','fontsize',30)            
set(ylhand,'string','Final BER','fontsize',30)
legend('Our work','Subspace patterns', 'Random patterns')
%==========================================================================
    

