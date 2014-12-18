%==========================================================================
%***********************FUNCTION: neural_journal_v1******************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% N_const_per_job: Number of constraints learned by each submitted job to the cluster
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% error_bits: Number of initial noisy nodes. It can be a vector
% max_noise_amp: The maximum amplitude of noise during the recall process
% no_simulated_instances: The number of patterns considered in the recall phase.
% simul_instances_per_job: Number of simulated instances in the recall phase per submitted job to the cluster
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a neural network and performs the
% learning and recall phases of a bipartite neural associative memory. More 
% specifically, the function first learns a pre-specified number of 
% constraints. Then the function goes directly to the recall phase. 
% Each step is performed by dividing the jobs in parallel and submitting
% them to the cluster. 
% The learning algorithm is based on the approach mentioned in our journal
% paper. The recall process considers three methods: The winner-take-all,
% the original bit flipping (mentioned in our ITW 2011 paper) and the
% simplified bit-flipping (proposed in ISIT 2012).
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


%%
%=============================INITIALIZATION===============================
if (~exist('initialization_done_by_compare','var'))    
    N = 300;                                % N is the number of neurons in network.
    K = 150;                                % K is the number of message bits.
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.
    index_max = 50;                         % This is the maximum number of random scenarios generated for simulation
    random_flag = 1;                        % Determines if patterns are drawn from a subspace or generated randomly
    no_of_patterns = 200;                   % The number of patterns in the dataset    
    max_noise_amp = 1;                      % Maximum value of integer-valued noise added to each bit  
    err_bits_range = [0:10];                % The number of bits that will be corrupted initially for the recall phase    
    
    a=clock;                                % Initialize the seed for random number generation with the clock value.
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a)));
    plot_flag = 1;                          % Plot the graphs
end
initialization_done_by_jankowski = 1;       % Tell the other codes that the initialization has been done
addpath(genpath('../../Common_Library'));                          % Include the library of common functions in the search path
%==========================================================================   

%%
for random_flag = 0:1
%============================THE LEARNING PHASE============================

%-----------------Check If the Learning Phase Has Been Done----------------
learn_flag = 1;
for l =1:index_max                                
    if (random_flag)
        fid = fopen(['../../Learn_Results/Jankowski/N_',num2str(N)...
            ,'_Random/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(l),'.mat'], 'r');
        if (fid < 0)
            learn_flag = 0;
            break
        else
            fclose(fid);
        end
    else
        fid = fopen(['../../Learn_Results/Jankowski/N_',num2str(N)...
            ,'_K_',num2str(K),'/Learned_Data_Capac_',num2str(no_of_patterns),'_index_',num2str(l),'.mat'], 'r');
        if (fid < 0)
            learn_flag = 0;
            break
        else
            fclose(fid);
        end
    end
end
%--------------------------------------------------------------------------

%----------------Perform the Learning Phase If Necessary-------------------
if (learn_flag == 0)
    run Complex_Neur_Jankowski_Learn
end
%--------------------------------------------------------------------------

%==========================================================================

%%
%=============================THE RECALL PHASE=============================

%-----------------Check If the Learning Phase Has Been Done----------------
recall_flag = 0;
if (random_flag)
    fid = fopen(['../../Recall_Results/Jankowski/N_',num2str(N),...
            '_Random/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'r');
    if (fid > 0)
        recall_flag = 1;
    end
else
    fid = fopen(['../../Recall_Results/Jankowski/N_',num2str(N),'_K_',num2str(K),...
            '/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'r');  
    if (fid > 0)
        recall_flag = 1;
    end
end

%--------------------------------------------------------------------------

%----------------Perform the Learning Phase If Necessary-------------------
if (recall_flag == 0)
    run Complex_Neur_Jankowski_Recall
end
%--------------------------------------------------------------------------

%==========================================================================

%%
%===========================READ RECALL RESULTS============================
run Read_jankowski_results
%==========================================================================
end

%============================PLOT THE RESULTS==============================
if plot_flag
    figure;
    plot(sort(processed_error_bits_Jankowski_non_random),sort(processed_BER_Jankowski_non_random),'b-.','LineWidth',2,'Color','red');
    hold on        
    plot(sort(processed_error_bits_Jankowski_random),sort(processed_BER_Jankowski_random),'b-*','LineWidth',2,'Color','blue');
    legend('Subspace patterns', 'Random patterns')
    set(gca,'Interpreter','latex')           
    xlhand = get(gca,'xlabel');            
    ylhand = get(gca,'ylabel');            
    set(xlhand,'string','e','fontsize',30)            
    set(ylhand,'string','Final BER','fontsize',30)
    title('Jankowski''s work')
end
%==========================================================================
