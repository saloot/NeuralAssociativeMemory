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
if (~exist('initialization_done','var'))    
    N = 300;                                % N is the number of neurons in network.
    K = 150;                                % K is the number of message bits.
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.
    index_max = 50;                         % This is the maximum number of random scenarios generated for simulation
    random_flag = 1;                        % Determines if patterns are drawn from a subspace or generated randomly
    no_of_patterns = 200;                   % The number of patterns in the dataset
    no_of_simulated_instance = 500;         % The number of patterns that are going to be denoised during the recall phase
    max_noise_amp = 1;                      % Maximum value of integer-valued noise added to each bit  
    err_bits_range = [0:10];                % The number of bits that will be corrupted initially for the recall phase
    addpath(genpath('../../Common_Library'));                          % Include the library of common functions in the search path
    
    a=clock;                                % Initialize the seed for random number generation with the clock value.
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
end

figure_legend = [];                         % This variable will be field with the legend for the graphs
hist_x_axis_itr = [0:20:1000];             % This is bin places for total learning iterations histogram
hist_x_axis_sparsity = [0:.05:1];           % This is bin places for total sparsity percentage

hist_x_itr = [];                            % This is output bin places produced BY the hist command for total learning iterations (to be used in bar plots)
hist_x_sparsity = [];                       % This is output bin places produced BY the hist command for sparsity percentage (to be used in bar plots)

hist_out_itr = [];                          % This is output values produced bY the hist command for total learning iterations (to be used in bar plots)
hist_out_sparsity = [];                     % This is output values produced bY the hist command for sparsity percentage (to be used in bar plots)

%==========================================================================

%%
%============================PROCESS THE RESULTS===========================
if (random_flag)
    processed_error_bits_Jankowski_random = [];
    processed_PER_Jankowski_random = [];
    processed_BER_Jankowski_random = [];
    processed_count_Jankowski_random = [];
else
    processed_error_bits_Jankowski_non_random = [];
    processed_PER_Jankowski_non_random = [];
    processed_BER_Jankowski_non_random = [];
    processed_count_Jankowski_non_random = [];
end    

if (random_flag)
    fid = fopen(['../../Recall_Results/Jankowski/N_',num2str(N),...
                '_Random/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'r');
else
    fid = fopen(['../../Recall_Results/Jankowski/N_',num2str(N),'_K_',num2str(K),...
                 '/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'r');  
end
    
if (fid > -1)
    results_Jankowski = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
    fclose(fid);                 
else
    error('Undefined Jankowski input file!')
end            
unprocessed_error_bits_Jankowski = results_Jankowski(2,:);
unprocessed_PER_Jankowski = results_Jankowski(6,:);
unprocessed_BER_Jankowski = results_Jankowski(10,:);

for i = 1:length(unprocessed_error_bits_Jankowski)
    
    if (random_flag)
        processed_error_bits_Jankowski = processed_error_bits_Jankowski_random;
    else
        processed_error_bits_Jankowski = processed_error_bits_Jankowski_non_random;
    end
    
    processed_flag = 0;
    for j = 1:length(processed_error_bits_Jankowski)
        if (unprocessed_error_bits_Jankowski(i) == processed_error_bits_Jankowski(j))
            processed_flag = 1;
            break;
        end
    end
    
    if (processed_flag == 0)
        if (random_flag)
            processed_error_bits_Jankowski_random = [processed_error_bits_Jankowski_random,unprocessed_error_bits_Jankowski(i)];
            processed_BER_Jankowski_random = [processed_BER_Jankowski_random,unprocessed_BER_Jankowski(i)];
            processed_PER_Jankowski_random = [processed_PER_Jankowski_random,unprocessed_PER_Jankowski(i)];
            processed_count_Jankowski_random = [processed_count_Jankowski_random,1];
        else
            processed_error_bits_Jankowski_non_random = [processed_error_bits_Jankowski_non_random,unprocessed_error_bits_Jankowski(i)];
            processed_BER_Jankowski_non_random = [processed_BER_Jankowski_non_random,unprocessed_BER_Jankowski(i)];
            processed_PER_Jankowski_non_random = [processed_PER_Jankowski_non_random,unprocessed_PER_Jankowski(i)];
            processed_count_Jankowski_non_random = [processed_count_Jankowski_non_random,1];
        end
    else
        if (random_flag)
            processed_BER_Jankowski_random(j) = processed_BER_Jankowski_random(j) + unprocessed_BER_Jankowski(i);
            processed_PER_Jankowski_random(j) = processed_PER_Jankowski_random(j) + unprocessed_PER_Jankowski(i);
            processed_count_Jankowski_random(j) = processed_count_Jankowski_random(j) + 1;
        else                    
            processed_BER_Jankowski_non_random(j) = processed_BER_Jankowski_non_random(j) + unprocessed_BER_Jankowski(i);
            processed_PER_Jankowski_non_random(j) = processed_PER_Jankowski_non_random(j) + unprocessed_PER_Jankowski(i);
            processed_count_Jankowski_non_random(j) = processed_count_Jankowski_non_random(j) + 1;
        end
        
    end
end

if random_flag   
    processed_BER_Jankowski_random = processed_BER_Jankowski_random./processed_count_Jankowski_random;                        
    processed_PER_Jankowski_random = processed_PER_Jankowski_random./processed_count_Jankowski_random;
else
    processed_BER_Jankowski_non_random = processed_BER_Jankowski_non_random./processed_count_Jankowski_non_random;                        
    processed_PER_Jankowski_non_random = processed_PER_Jankowski_non_random./processed_count_Jankowski_non_random;
end
%==========================================================================