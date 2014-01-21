%==========================================================================
%**************FUNCTION: clustered_neural_initialization*******************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% mean_cluster_size: The average number of pattern nodes in a cluster
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

function clustered_neural_initialization(N,L,mean_cluster_size)

%=============================INITIALIZATION===============================
        
%---------------------------Network Parameters-----------------------------

%--------------------------------------------------------------------------


%----------------------------Other Initializations-------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
mkdir(['/scratch/amir/Clustered_Neural/Initialization_Files'],['N_',num2str(N),'_L_',num2str(L)]);        % Create a specific folder for the current N and K

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%==========================================================================



%%
%==================CREATE THE PATTERN GENERATOR MATRIX=====================

%-------Determine Which Pattern Nodes Are Connected to Each Clusters-------    
for l = 1:L       
    cluster_size = max(40,mean_cluster_size + floor(sqrt((mean_cluster_size))*randn));
    index_pattern_neurons = randperm(N);
    index_pattern_neurons = index_pattern_neurons(1:cluster_size);                       
                   
    %-----Save the Neighborhood of Each Cluster in a Separate File-----            
    save(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_L_',num2str(L),...
    '/Clustered_parameters_N_',num2str(N),'_L_',num2str(L),'_cluster_index_',num2str(l),'.mat'],'index_pattern_neurons');
    %------------------------------------------------------------------
        
end        
%==========================================================================

