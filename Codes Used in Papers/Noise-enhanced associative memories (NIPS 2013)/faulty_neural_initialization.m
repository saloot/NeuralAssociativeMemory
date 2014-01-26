%==========================================================================
%***************FUNCTION: faulty_neural_initialization*********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% deg_cluster_min The minimum number of clusters that a pattern node is connected to
% index_max: The maximum number of random setups constructed for simulations
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This piece of code initializes the necessary parameters for the clustered
% neural associative memory and stores them on appropriate files so that
% other functions can have access to them. The advantage of storing them on
% hard disk is that we can run a parallel version of the code which makes
% it much more faster.

%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function faulty_neural_initialization(N,K,L,deg_cluster_min,index_max)

%=============================INITIALIZATION===============================
if (~exist('initialization_done','var'))    % If already not initialized by the GUI...
    
    %-------------------------Network Parameters---------------------------    
    N_const = N-K;                          % N_cost represents the number of constraints.
    N_tot = N*(L)/deg_cluster_min;                            % This is the total length of the pattern
    K_tot = K;       
    deg_column_G_tot = 6;                       % This is the average degree of each column for the global pattern generator matrix.
    deg_row_G_tot = N_tot*deg_column_G_tot/K_tot;                          % This is the average degree of each row for the global pattern generator matrix.
    %----------------------------------------------------------------------

    %-------------------------Simulation Parameters------------------------    
    learn_itr_max = 1000;                 % This is the number of times that the learning phase is repeated for the patterns in the training set       
    pattern_learn_number = min(100000,2^K_tot);           % This is the number of patterns used in the learning process.    
    %----------------------------------------------------------------------

    
    %--------------------------Neural Parameters---------------------------
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.       
    %----------------------------------------------------------------------       
end

%----------------------------Other Initializations-------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
mkdir(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files'],['N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)]);        % Create a specific folder for the current N and K

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%==========================================================================



%%
%==================CREATE THE PATTERN GENERATOR MATRIX=====================
for index = 1:index_max
    G_tot = bipartite(N_tot,K_tot,deg_column_G_tot,deg_row_G_tot);              % Generate the pattern generator matrix
    clus_index_node = zeros(N_tot,deg_cluster_min);                             % Contains the cluster each pattern node is connected to

    %-------Determine Which Clusters Each Pattern Node Is Conneced to------
    for i = 1:N_tot                 
        temp = randperm(L);
        clus_index_node(i,:) = temp(1:deg_cluster_min);             
    end    
    %----------------------------------------------------------------------
    
    %-----Determine Which Pattern Nodes Are Connected to Each Clusters-----
    for l = 1:L
        index_l = [];           
        for i = 1:N_tot        
            for j = 1:deg_cluster_min                    
                if (clus_index_node(i,j) == l)                                                   
                    index_l = [index_l,i];                                      % Containt the pattern nodes in cluster l
                    break;                    
                end                
            end            
        end
        %-----Save the Neighborhood of Each Cluster in a Separate File-----
        save(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/faulty_neural_parameters_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat'],'index_l');
        %------------------------------------------------------------------
        
    end    
    %----------------------------------------------------------------------

    %-------------Generate the Training Data Set for Each Setup------------
    KK=round(log(pattern_learn_number)/log(2));
    mu_list = zeros(KK+1,pattern_learn_number);
    for i = 1:pattern_learn_number
        randindex = randperm(K_tot);
        mu_list(1,i) = 1+floor((pattern_learn_number-1)*rand);
        mu_list(2:KK+1,i) = randindex(1:KK);      
    end
    %----------------------------------------------------------------------
    
    %----------------------Save Initialized Variables----------------------
    clear temp
    clear index_l 
    save(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/faulty_neural_parameters_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);                    % Store all the variables in the appropriate file. 
    %----------------------------------------------------------------------
        
end
%==========================================================================

