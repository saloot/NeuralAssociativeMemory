%==========================================================================
%********************************READ ME***********************************
%==========================================================================

%--------------------------------Summary-----------------------------------
% This piece of code initializes the necessary parameters for the constraint
% enforcing neural network and stores them on appropriate files so that
% other functions can have access to them. The advantage of storing them on
% hard disk is that we can run a parallel version of the code which makes
% it much more faster.
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================

%=============================INITIALIZATION===============================
clc
if (~exist('initialization_done','var'))    % If already not initialized by the GUI...
    %-------------------------Network Parameters---------------------------
    N = 300;                                % N is the number of neurons in network.
    K = 100;                                 % K is the number of message bits.
    L = 4;                                  % L is the number of clusters of the subpatterns
    N_const = N-K;                          % N_cost represents the number of constraints.
    N_tot = N*L;                            % This is the total length of the pattern
    N_const_tot = N_tot-K;                  % This is the total number of constraints
    deg_row_W = 2*round(log(N));     % This is the average degree of each row for the weight matrix.
    deg_row_W_tot = 2*round(log(N_const_tot)); % This is the average degree of each column for the total weight matrix.
    deg_column_G = 2;                       % This is the average degree of each column for the pattern generator matrix.
    deg_row_G = 6;                          % This is the average degree of each row for the pattern generator matrix.
    %----------------------------------------------------------------------

    %-------------------------Simulation Parameters------------------------    
    learn_itr_max = 100000;                 % This is the number of times that the learning phase is repeated for the patterns in the training set       
    pattern_learn_number = min(10000,10*(2^K));           % This is the number of patterns used in the learning process.
    %----------------------------------------------------------------------

    
    %--------------------------Neural Parameters---------------------------
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.   
    no_of_patterns = 2^K;                   % This is the number of patterns we are going to memorize.
   
    
    
    y = zeros(1,N_const);                   % y is the state of constraint nodes.
    
    alph = 1/(deg_column_G*N*deg_row_W);% alph is the learning rate.    
    alph_tot = 1/(deg_column_G*N_tot*deg_row_W_tot);% alph is the learning rate.    
    betta = 0;                               % betta is a variable that controls computing the novelty part of a pattern with respect to the current basis vectors. 
    epsilon = 1e-4;                         % epsilon is the desired accuracy of the algorithm. 
    sigm = 100;                            % sigm controls the sparsity of the network.
    lambda_tot = 0.001;                             % lambda is step size for sparisty constraint.         
    lambda = 0.001;                             % lambda is step size for sparisty constraint.         
    %----------------------------------------------------------------------
        
end

%----------------------------Other Initializations-------------------------
N_orig = N;                                 % N_orig keeps the original value of N.
K_orig = K;                                 % To save the original value as K is updated during the algorithm.
%--------------------------------------------------------------------------

%==========================================================================


%%
%===================CREATE THE PATTERN GENERATOR MATRIX====================
G_tot = [];                                 % G_tot is the generator matrix for all the patterns.
G = [];
for l = 1:L
    GG = bipartite(N,K,deg_column_G,deg_row_G);  % G is the generator matrix used to generate the patterns from binary vectors.
    G = [G;GG];
    G_tot = [G_tot,GG];    
end
%==========================================================================

%%
%===================CREATE THE TRAINING PATTERN LIST=======================            
for i = 1:pattern_learn_number
    mu = 1+floor( (no_of_patterns-1)*rand);             % If it is the first learning iteration, generate a random pattern index...
    mu_list(i) = mu;                             % ...and store the index in a list
end
%==========================================================================

%%
%=================STORE NECESSARY VALUES IN A FILE=========================

%----------------------Save Initialized Variables--------------------------
save(['neural_learn_parameters_v3_N_',num2str(N),'_K_',num2str(K),'.mat']);                    % Store all the variables in the appropriate file. 
%--------------------------------------------------------------------------

%==========================================================================

