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

%%
%=============================INITIALIZATION===============================
if (~exist('initialization_done','var'))    % If already not initialized by the GUI...
    %-------------------------Network Parameters---------------------------
    N = 1600;                                % N is the number of neurons in network.
    K = 800;                                % K is the number of message bits.
    N_const = N-K;                          % N_cost represents the number of constraints.        
    deg_column_G = 15;
    deg_row_G = N*deg_column_G/K;    
    %----------------------------------------------------------------------

    %-------------------------Simulation Parameters------------------------    
    learn_itr_max = 10000;                             % This is the number of times that the learning phase is repeated for the patterns in the training set       
    index_max = 50;                                     % This is the maximum number of random scenarios generated for simulation
    pattern_learn_number = min(100000,2^K);             % This is the number of patterns used in the learning process.
    KK = round(log(pattern_learn_number)/log(2));
    %----------------------------------------------------------------------

    
    %--------------------------Neural Parameters---------------------------
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.                             
    %----------------------------------------------------------------------
        
end

%----------------------------Other Initializations-------------------------
mkdir(['/scratch/amir/ITW_Journal/Initialization_Files'],['N_',num2str(N),'_K_',num2str(K)]);    % Create a specific folder for the current N and K
addpath('/home1/amir/cluster/Common_Library');                              % Include the library of common functions in the search path

a=clock;                                                                    % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%==========================================================================


%%
%============CREATE THE PATTERN GENERATOR MATRIX & TRAINING SETS===========
for index = 1:index_max

    %----------------Create the Pattern Generating Matrix------------------
    G = bipartite(N,K,deg_column_G,deg_row_G);
    %----------------------------------------------------------------------

    %-----------------Create the Training List Indices---------------------
    mu_list = 1+floor( (pattern_learn_number-1)*rand(1,pattern_learn_number));     
    %----------------------------------------------------------------------
    
    %-----------------------Save Training Sets-----------------------------
    save(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(index),'.mat'],'G','mu_list');        
    %----------------------------------------------------------------------
end

%%
%=================STORE NECESSARY VALUES IN A FILE=========================
%----------------------Clear Unnecessary Variables-------------------------
clear train_set X train_set_local;
%--------------------------------------------------------------------------

%----------------------Save Initialized Variables--------------------------
save(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_parameters_v1_N_',num2str(N),'_K_',num2str(K),'.mat']);                    % Store all the variables in the appropriate file. 
%--------------------------------------------------------------------------
%==========================================================================

