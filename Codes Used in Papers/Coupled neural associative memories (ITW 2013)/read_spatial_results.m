%==========================================================================
%******************FUNCTION: read_cluster_results**************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N_horiz_in: The number of columns in the 2D patterns
% N_vert_in: The number of rows in the 2D patterns
% L_horiz_in: The number of clusters within a neural plane.
% L_vert_in: The number of neural planes
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
% recall_algorithm_option: The recall algorithm identifier (0 for winner-take-all, 1 for the original bit flipping and 2 for the simplified bit flipping
% try_max: The maximum number of iterations in outside-cluster recall algorithm
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% processed_error_bits: The list of processed number of initial erroneous nodes
% processed_PER: The list of processed Pattern Error Rates
% processed_error_bits_fixed: The list of processed number of initial erroneous nodes for the coupled system in which the boundaries are fixed
% processed_PER_fixed: The list of processed Pattern Error Rates for the coupled system in which the boundaries are fixed
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a spatially-coupled neural associative
% memory and reads the result of recall phase from the appropriate files. 
% The results will then be plotted and compared with theoretical values. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


% function [processed_error_bits,processed_PER,processed_error_bits_fixed,processed_PER_fixed] = read_spatial_results(N_horiz_in,N_vert_in,L_horiz_in,L_vert_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option)
warning off

N_horiz_in = 64;
N_vert_in = 64;
L_horiz_in = 29;
L_vert_in = 29;
alpha0 = 0.5;
beta0 = 0.5;
theta0 = 0.5;
gamma_BFO = 0.95;
gamma_BFS = 0.95;
recall_algorithm_option = 2;

%%
%==============================INITIALIZATION==============================
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
figure_legend = [];                         % This variable will be field with the legend for the graphs
% hist_x_axis_itr = [0:20:1000];             % This is bin places for total learning iterations histogram
%==========================================================================

%%
%============================PROCESS THE RESULTS===========================
for ll = 1:length(N_horiz_in)
    N_horiz = N_horiz_in(ll);
    N_vert = N_vert_in(ll);
    L_horiz = L_horiz_in(ll);
    L_vert = L_vert_in(ll);
    
    %----------------Load the Saved Initialized Parameters---------------------
    load(['/scratch/amir/Spatially_Coupled/Initialization_Files/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_parameters_N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'.mat']);           
    %--------------------------------------------------------------------------

                                
    %---------------------------Read Recall Results------------------------
    if (recall_algorithm_option == 0)
        fid = fopen(['/scratch/amir/Spatially_Coupled/Recall_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'r'); 
        if (fid > -1)
            results = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);
            fclose(fid);
        else
            error('Undefined WTA input file!')    
        end
        
        fid = fopen(['/scratch/amir/Spatially_Coupled/Recall_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_results_fixed_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'r'); 
        if (fid > -1)
            results_fixed = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);
            fclose(fid);
        else
            error('Undefined BFO input file!')    
        end
    elseif (recall_algorithm_option == 2)
        fid = fopen(['/scratch/amir/Spatially_Coupled/Recall_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'r'); 
        if (fid > -1)
            results = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);
            fclose(fid);
        else
            error('Undefined BFS input file!')    
        end
        
        fid = fopen(['/scratch/amir/Spatially_Coupled/Recall_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_results_fixed_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'r'); 
        if (fid > -1)
            results_fixed = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);
            fclose(fid);
        else
            error('Undefined BFO input file!')    
        end
        
    elseif (recall_algorithm_option == 1)
        fid = fopen(['/scratch/amir/Spatially_Coupled/Recall_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'r'); 
        if (fid > -1)
            results = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);
            fclose(fid);
        else
            error('Undefined BFO input file!')    
        end
        
        fid = fopen(['/scratch/amir/Spatially_Coupled/Recall_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_results_fixed_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'r'); 
        if (fid > -1)
            results_fixed = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);
            fclose(fid);
        else
            error('Undefined BFO input file!')    
        end
        
    else
        error('Unknown recall algorithm');
    end
    %----------------------------------------------------------------------                       
    
    %---------------------Process the Results------------------------------
    unprocessed_error_bits = results(2,:);
    unprocessed_PER = results(6,:);
    unprocessed_BER = results(10,:);
    
    processed_error_bits = [0];
    processed_PER = [0];
    processed_BER = [0];
    processed_count = [1];
    for i = 1:length(unprocessed_error_bits)
        processed_flag = 0;
        for j = 1:length(processed_error_bits)
            if (unprocessed_error_bits(i) == processed_error_bits(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits = [processed_error_bits,unprocessed_error_bits(i)];
            processed_PER = [processed_PER,unprocessed_PER(i)];
            processed_BER = [processed_BER,unprocessed_BER(i)];
            processed_count = [processed_count,1];
        else
            processed_PER(j) = processed_PER(j) + unprocessed_PER(i);
            processed_BER(j) = processed_BER(j) + unprocessed_BER(i);
            processed_count(j) = processed_count(j) + 1;
        end
    end
    processed_PER = processed_PER./processed_count;
    processed_BER = processed_BER./processed_count;
    %----------------------------------------------------------------------                       
    
    %--------------------Process the Fixed Results-------------------------
    unprocessed_error_bits_fixed = results_fixed(2,:);
    unprocessed_PER_fixed = results_fixed(6,:);
    unprocessed_BER_fixed = results_fixed(10,:);
    
    processed_error_bits_fixed = [0];
    processed_PER_fixed = [0];
    processed_BER_fixed = [0];
    processed_count_fixed = [1];
    for i = 1:length(unprocessed_error_bits_fixed)
        processed_flag = 0;
        for j = 1:length(processed_error_bits_fixed)
            if (unprocessed_error_bits_fixed(i) == processed_error_bits_fixed(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits_fixed = [processed_error_bits_fixed,unprocessed_error_bits_fixed(i)];
            processed_PER_fixed = [processed_PER_fixed,unprocessed_PER_fixed(i)];
            processed_BER_fixed = [processed_BER_fixed,unprocessed_BER_fixed(i)];
            processed_count_fixed = [processed_count_fixed,1];
        else
            processed_PER_fixed(j) = processed_PER_fixed(j) + unprocessed_PER_fixed(i);
            processed_BER_fixed(j) = processed_BER_fixed(j) + unprocessed_BER_fixed(i);
            processed_count_fixed(j) = processed_count_fixed(j) + 1;
        end
    end
    processed_PER_fixed = processed_PER_fixed./processed_count_fixed;
    processed_BER_fixed = processed_BER_fixed./processed_count_fixed;
    %----------------------------------------------------------------------                       
    
    
    
     

    %---------------------Constructu the Legend----------------------------
    temp_leg = ['n=',num2str(N_vert*N_horiz),', L=',num2str(L_vert), ', D=',num2str(L_horiz)];    
    for i = length(num2str(N_vert)):length(num2str(max(N_vert_in)))
        temp_leg = [temp_leg,' '];
    end
    for i = length(num2str(L_vert)):length(num2str(max(L_vert_in)))
        temp_leg = [temp_leg,' '];        
    end            
    
    temp_lega = [temp_leg,' - Unconstrained'];                              % 'Uncostrained' refers to a system in which the boundaries are not fixed
    figure_legend = [figure_legend;temp_lega];
    temp_lega = [temp_leg,' - Constrained  '];                              % 'Costrained' refers to a system in which the boundaries are fixed
    figure_legend = [figure_legend;temp_lega];
    %----------------------------------------------------------------------

    %----------------------Display The Results on Graph--------------------
    plot(sort(processed_error_bits),sort(processed_PER),'b','LineWidth',2,'Color',[tanh(ll/length(N_horiz_in)),ll/length(N_horiz_in),1-(ll/length(N_horiz_in))]);
    hold on
    plot(sort(processed_error_bits_fixed),sort(processed_PER_fixed),'b-.','LineWidth',2,'Color',[tanh(ll/length(N_horiz_in)),ll/length(N_horiz_in),1-(ll/length(N_horiz_in))]);
    hold on
%     title('PER for WTA')

%     figure(h_error_BER_WTA);
%     plot(sort(processed_error_bits),sort(BER_tehroy),'b','LineWidth',2,'Color',[tanh(ll/length(N)),ll/length(N),1-(ll/length(N))]);
%     hold on
%     plot(sort(processed_error_bits),sort(processed_BER),'b--','LineWidth',2,'Color',[tanh(ll/length(N)),ll/length(N),1-(ll/length(N))]);    
%     title('BER for WTA')
    %----------------------------------------------------------------------
    
    
end
%==========================================================================

%======================ADD THE LEGEND TO THE FIGURES=======================
% figure(h_error_PER)
legend(figure_legend)
% figure(h_error_BER_BFO)
% legend(figure_legend)
% figure(h_error_BER_WTA)
% legend(figure_legend)
%==========================================================================

%======================SAVE THE PROCESSED RESULTS==========================
% mkdir(['/scratch/amir/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)],'Read_Results');        % Create a specific folder for the current N and K
% save(['/scratch/amir/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Read_Results/final_results_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_alpha0_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '.mat'],'PER_tehroy_single','processed_error_bits',...
%     'processed_PER','processed_BER','lambda','rho','deg_column','deg_row','processed_count','m');
%==========================================================================