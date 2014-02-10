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


N_in = 400;
K_in = 200;
alpha0 = 0.75;
beta0 = 1;
theta0 = .031;
N_const_per_job = 20;
error_bits = [15:5:50];
max_noise_amp = 1;
gamma_BFO = 1;
gamma_BFS = 0.99;
no_simulated_instances = 2000;
simul_instances_per_job = no_simulated_instances;
mca_flag = 0;


% N_in = 800;
% K_in = 400;
% alpha0= 0.95;
% beta0 = 1;
% theta0 = 0.029;

% gamma_BFO = 0.99;
% gamma_BFS = 0.99;

% gamma_BFO = 0.95;
% gamma_BFS = 0.95;

% function [processed_error_bits_BFO,processed_PER_WTA,processed_PER_BFO,processed_PER_BFS] = read_journal_results(N_in,K_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS)
warning off
%%
%==============================INITIALIZATION==============================
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
% h_itr = figure;    
h_sparse=figure;    
for ll = 1:length(N_in)
    N = N_in(ll);
    K = K_in(ll);
    
    load(['/scratch/amir/ITW_Journal/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_parameters_v1_N_',num2str(N),'_K_',num2str(K),'.mat']);           
    addpath(genpath('/home1/amir/cluster/Common_Library'));   
    addpath('/home1/amir/cluster/ITW_Journal/Verify_Theory');   
    
   
    
    %---------------------------Read Recall Results------------------------
    fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '.txt'], 'r'); 
    if (fid > -1)
        results_WTA = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
        fclose(fid);
    else
        error('Undefined WTA input file!')
    
    end

    fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '.txt'], 'r'); 
    if (fid > -1)
        results_BFS = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
        fclose(fid);
    else
        error('Undefined BFS input file!')
    
    end

    fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '.txt'], 'r'); 
    if (fid > -1)
        results_BFO = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
        fclose(fid);
    else
        error('Undefined BFO input file!')
    
    end
    
    
    fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '_first_itr.txt'], 'r'); 
    if (fid > -1)
        results_WTA_first_itr = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
        fclose(fid);
    else
        error('Undefined WTA input file!')
    
    end

    fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '_first_itr.txt'], 'r'); 
    if (fid > -1)
        results_BFS_first_itr = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
        fclose(fid);
    else
        error('Undefined BFS input file!')
    
    end

    fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/N_',num2str(N),'_K_',num2str(K),'/neural_journal_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '_first_itr.txt'], 'r'); 
    if (fid > -1)
        results_BFO_first_itr = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
        fclose(fid);
    else
        error('Undefined BFO input file!')
    
    end
    %----------------------------------------------------------------------                       
        
    %---------------------Process the Results------------------------------
    unprocessed_error_bits_WTA = results_WTA(2,:);
    unprocessed_PER_WTA = results_WTA(6,:);
    unprocessed_BER_WTA = results_WTA(10,:);
    unprocessed_error_bits_BFO = results_BFO(2,:);
    unprocessed_PER_BFO = results_BFO(6,:);
    unprocessed_BER_BFO = results_BFO(10,:);
    unprocessed_error_bits_BFS = results_BFS(2,:);
    unprocessed_PER_BFS = results_BFS(6,:);
    unprocessed_BER_BFS = results_BFS(10,:);

    processed_error_bits_WTA = [0];
    processed_PER_WTA = [0];
    processed_BER_WTA = [0];
    processed_error_bits_BFO = [0];
    processed_PER_BFO = [0];
    processed_BER_BFO = [0];
    processed_error_bits_BFS = [0];
    processed_PER_BFS = [0];
    processed_BER_BFS = [0];
    processed_count_WTA = [1];
    processed_count_BFO = [1];
    processed_count_BFS = [1];

    unprocessed_error_bits_WTA_first_itr = results_WTA_first_itr(2,:);
    unprocessed_PER_WTA_first_itr = results_WTA_first_itr(6,:);
    unprocessed_BER_WTA_first_itr = results_WTA_first_itr(10,:);
    unprocessed_error_bits_BFO_first_itr = results_BFO_first_itr(2,:);
    unprocessed_PER_BFO_first_itr = results_BFO_first_itr(6,:);
    unprocessed_BER_BFO_first_itr = results_BFO_first_itr(10,:);
    unprocessed_error_bits_BFS_first_itr = results_BFS_first_itr(2,:);
    unprocessed_PER_BFS_first_itr = results_BFS_first_itr(6,:);
    unprocessed_BER_BFS_first_itr = results_BFS_first_itr(10,:);

    processed_error_bits_WTA_first_itr = [0];
    processed_PER_WTA_first_itr = [0];
    processed_BER_WTA_first_itr = [0];
    processed_error_bits_BFO_first_itr = [0];
    processed_PER_BFO_first_itr = [0];
    processed_BER_BFO_first_itr = [0];
    processed_error_bits_BFS_first_itr = [0];
    processed_PER_BFS_first_itr = [0];
    processed_BER_BFS_first_itr = [0];
    processed_count_WTA_first_itr = [1];
    processed_count_BFO_first_itr = [1];
    processed_count_BFS_first_itr = [1];
    
    for i = 1:length(unprocessed_error_bits_WTA)
        processed_flag = 0;
        for j = 1:length(processed_error_bits_WTA)
            if (unprocessed_error_bits_WTA(i) == processed_error_bits_WTA(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits_WTA = [processed_error_bits_WTA,unprocessed_error_bits_WTA(i)];
            processed_PER_WTA = [processed_PER_WTA,unprocessed_PER_WTA(i)];
            processed_BER_WTA = [processed_BER_WTA,unprocessed_BER_WTA(i)];
            processed_count_WTA = [processed_count_WTA,1];
        else
            processed_PER_WTA(j) = processed_PER_WTA(j) + unprocessed_PER_WTA(i);
            processed_BER_WTA(j) = processed_BER_WTA(j) + unprocessed_BER_WTA(i);
            processed_count_WTA(j) = processed_count_WTA(j) + 1;
        end
    end
    processed_PER_WTA = processed_PER_WTA./processed_count_WTA;
    processed_BER_WTA = processed_BER_WTA./processed_count_WTA;

    for i = 1:length(unprocessed_error_bits_BFO)
        processed_flag = 0;
        for j = 1:length(processed_error_bits_BFO)
            if (unprocessed_error_bits_BFO(i) == processed_error_bits_BFO(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits_BFO = [processed_error_bits_BFO,unprocessed_error_bits_BFO(i)];
            processed_PER_BFO = [processed_PER_BFO,unprocessed_PER_BFO(i)];
            processed_BER_BFO = [processed_BER_BFO,unprocessed_BER_BFO(i)];
            processed_count_BFO = [processed_count_BFO,1];
        else
            processed_PER_BFO(j) = processed_PER_BFO(j) + unprocessed_PER_BFO(i);
            processed_BER_BFO(j) = processed_BER_BFO(j) + unprocessed_BER_BFO(i);
            processed_count_BFO(j) = processed_count_BFO(j) + 1;
        end
    end
    processed_PER_BFO = processed_PER_BFO./processed_count_BFO;
    processed_BER_BFO = processed_BER_BFO./processed_count_BFO;

    for i = 1:length(unprocessed_error_bits_BFS)
        processed_flag = 0;
        for j = 1:length(processed_error_bits_BFS)
            if (unprocessed_error_bits_BFS(i) == processed_error_bits_BFS(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits_BFS = [processed_error_bits_BFS,unprocessed_error_bits_BFS(i)];
            processed_PER_BFS = [processed_PER_BFS,unprocessed_PER_BFS(i)];
            processed_BER_BFS = [processed_BER_BFS,unprocessed_BER_BFS(i)];
            processed_count_BFS = [processed_count_BFS,1];
        else
            processed_PER_BFS(j) = processed_PER_BFS(j) + unprocessed_PER_BFS(i);
            processed_BER_BFS(j) = processed_BER_BFS(j) + unprocessed_BER_BFS(i);
            processed_count_BFS(j) = processed_count_BFS(j) + 1;
        end
    end
    processed_PER_BFS = processed_PER_BFS./processed_count_BFS;
    processed_BER_BFS = processed_BER_BFS./processed_count_BFS;
    processed_error_bits = processed_error_bits_BFO;
    
    
    
    
    for i = 1:length(unprocessed_error_bits_WTA_first_itr)
        processed_flag = 0;
        for j = 1:length(processed_error_bits_WTA_first_itr)
            if (unprocessed_error_bits_WTA_first_itr(i) == processed_error_bits_WTA_first_itr(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits_WTA_first_itr = [processed_error_bits_WTA_first_itr,unprocessed_error_bits_WTA_first_itr(i)];
            processed_PER_WTA_first_itr = [processed_PER_WTA_first_itr,unprocessed_PER_WTA_first_itr(i)];
            processed_BER_WTA_first_itr = [processed_BER_WTA_first_itr,unprocessed_BER_WTA_first_itr(i)];
            processed_count_WTA_first_itr = [processed_count_WTA_first_itr,1];
        else
            processed_PER_WTA_first_itr(j) = processed_PER_WTA_first_itr(j) + unprocessed_PER_WTA_first_itr(i);
            processed_BER_WTA_first_itr(j) = processed_BER_WTA_first_itr(j) + unprocessed_BER_WTA_first_itr(i);
            processed_count_WTA_first_itr(j) = processed_count_WTA_first_itr(j) + 1;
        end
    end
    processed_PER_WTA_first_itr = processed_PER_WTA_first_itr./processed_count_WTA_first_itr;
    processed_BER_WTA_first_itr = processed_BER_WTA_first_itr./processed_count_WTA_first_itr;

    for i = 1:length(unprocessed_error_bits_BFO_first_itr)
        processed_flag = 0;
        for j = 1:length(processed_error_bits_BFO_first_itr)
            if (unprocessed_error_bits_BFO_first_itr(i) == processed_error_bits_BFO_first_itr(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits_BFO_first_itr = [processed_error_bits_BFO_first_itr,unprocessed_error_bits_BFO_first_itr(i)];
            processed_PER_BFO_first_itr = [processed_PER_BFO_first_itr,unprocessed_PER_BFO_first_itr(i)];
            processed_BER_BFO_first_itr = [processed_BER_BFO_first_itr,unprocessed_BER_BFO_first_itr(i)];
            processed_count_BFO_first_itr = [processed_count_BFO_first_itr,1];
        else
            processed_PER_BFO_first_itr(j) = processed_PER_BFO_first_itr(j) + unprocessed_PER_BFO_first_itr(i);
            processed_BER_BFO_first_itr(j) = processed_BER_BFO_first_itr(j) + unprocessed_BER_BFO_first_itr(i);
            processed_count_BFO_first_itr(j) = processed_count_BFO_first_itr(j) + 1;
        end
    end
    processed_PER_BFO_first_itr = processed_PER_BFO_first_itr./processed_count_BFO_first_itr;
    processed_BER_BFO_first_itr = processed_BER_BFO_first_itr./processed_count_BFO_first_itr;

    for i = 1:length(unprocessed_error_bits_BFS_first_itr)
        processed_flag = 0;
        for j = 1:length(processed_error_bits_BFS_first_itr)
            if (unprocessed_error_bits_BFS_first_itr(i) == processed_error_bits_BFS_first_itr(j))
                processed_flag = 1;
                break;
            end
        end
    
        if (processed_flag == 0)
            processed_error_bits_BFS_first_itr = [processed_error_bits_BFS_first_itr,unprocessed_error_bits_BFS_first_itr(i)];
            processed_PER_BFS_first_itr = [processed_PER_BFS_first_itr,unprocessed_PER_BFS_first_itr(i)];
            processed_BER_BFS_first_itr = [processed_BER_BFS_first_itr,unprocessed_BER_BFS_first_itr(i)];
            processed_count_BFS_first_itr = [processed_count_BFS_first_itr,1];
        else
            processed_PER_BFS_first_itr(j) = processed_PER_BFS_first_itr(j) + unprocessed_PER_BFS_first_itr(i);
            processed_BER_BFS_first_itr(j) = processed_BER_BFS_first_itr(j) + unprocessed_BER_BFS_first_itr(i);
            processed_count_BFS_first_itr(j) = processed_count_BFS_first_itr(j) + 1;
        end
    end
    processed_PER_BFS_first_itr = processed_PER_BFS_first_itr./processed_count_BFS_first_itr;
    processed_BER_BFS_first_itr = processed_BER_BFS_first_itr./processed_count_BFS_first_itr;
    processed_error_bits_first_itr = processed_error_bits_BFO_first_itr;
    %----------------------------------------------------------------------                       

    %------------------Read Cluster Degree Distribution--------------------
%     [deg_row,deg_column,lambda,rho] = deg_dribution(N,K,alpha0,beta0,theta0,index_max);    
    %----------------------------------------------------------------------                       
    
    %--------------------Select the approprate figure----------------------
    if (ll == 1)
%         h_error_PER_WTA = figure;    
        h_error_PER_BFO = figure;       
    else
        figure(h_error_PER_BFO);
    end
    %----------------------------------------------------------------------
    
    %--------------------------Compare with Theory-------------------------
%     [PER_tehroy_BFO,BER_tehroy_BFO] = Error_prob_ireg(N,N-K,gamma_BFO,deg_column,lambda,processed_error_bits_BFO);    
%     [PER_tehroy_WTA,BER_tehroy_WTA] = Error_prob_ireg(N,N-K,1,deg_column,lambda,processed_error_bits_WTA);
%         
%     PER_tehroy2_BFO = Error_prob_ireg2(N,N-K,gamma_BFO,deg_column,lambda,rho,processed_error_bits_BFO);
%     PER_tehroy2_BFO = error_theory(N,K,gamma_BFO,processed_error_bits_BFO,alpha0,beta0,theta0,index_max,2);
%     PER_tehroy2_WTA = Error_prob_ireg2(N,N-K,1,deg_column,lambda,rho,processed_error_bits_WTA);            
% 
%     PER_tehroy3_BFO = Error_prob_ireg3(N,N-K,gamma_BFO,deg_column,lambda,rho,processed_error_bits_BFO);    
%     PER_tehroy3_WTA = Error_prob_ireg3(N,N-K,1,deg_column,lambda,rho,processed_error_bits_WTA);    
    
%     try_max = 20;
    
%     PER_tehroy4_BFO = Error_prob_ireg4(N,N-K,gamma_BFO,deg_column,lambda,rho,processed_error_bits_BFO,try_max);    
%     PER_tehroy4_WTA = Error_prob_ireg4(N,N-K,1,deg_column,lambda,rho,processed_error_bits_WTA,try_max);    
    
%     [PER_tehroy5_BFO,BER_tehroy5_BFO5] = Error_prob_ireg5(N,N-K,gamma_BFO,deg_column,lambda,processed_error_bits_BFO,30);
%     [PER_tehroy5_WTA,BER_tehroy5_WTA] = Error_prob_ireg5(N,N-K,1,deg_column,lambda,processed_error_bits_BFO,30);
% 
%     [PER_tehroy6_BFO,BER_tehroy3_BFO] = new_bound(N,N-K,gamma_BFO,deg_column,lambda,processed_error_bits_BFO);
%     [PER_tehroy6_WTA,BER_tehroy3_WTA] = new_bound(N,N-K,1,deg_column,lambda,processed_error_bits_WTA);
    
%     PER_tehroy7_BFO = error_theory(N,K,gamma_BFO,processed_error_bits_BFO,alpha0,beta0,theta0,index_max,7);
%     PER_tehroy7_WTA = Error_prob_ireg7(N,N-K,1,deg_column,lambda,processed_error_bits_WTA,try_max);            
    %----------------------------------------------------------------------

    %---------------------Constructu the Legend----------------------------
    temp_leg = ['n=',num2str(N),', k=',num2str(K)];    
    for i = length(num2str(N)):length(num2str(max(N_in)))
        temp_leg = [temp_leg,' '];
    end
    for i = length(num2str(K)):length(num2str(max(K_in)))
        temp_leg = [temp_leg,' '];        
    end            
    
%     temp_lega = [temp_leg,' - Theory 1'];
%     figure_legend = [figure_legend;temp_lega];
    temp_lega = [temp_leg,' - Theory 2'];
    figure_legend = [figure_legend;temp_lega];
%     temp_lega = [temp_leg,' - Theory 3'];
%     figure_legend = [figure_legend;temp_lega];
%     temp_lega = [temp_leg,' - Theory 4'];
%     figure_legend = [figure_legend;temp_lega];
%     temp_lega = [temp_leg,' - Theory 5'];
%     figure_legend = [figure_legend;temp_lega];
%     temp_lega = [temp_leg,' - Theory 6'];
%     figure_legend = [figure_legend;temp_lega];
%     temp_lega = [temp_leg,' - Theory 7'];
%     figure_legend = [figure_legend;temp_lega];
    temp_lega = [temp_leg,' - Practice'];
    figure_legend = [figure_legend;temp_lega];
    temp_lega = [temp_leg,' - 1st Itr '];
    figure_legend = [figure_legend;temp_lega];
%     temp_lega = [temp_leg,' - BoundNew  '];
%     figure_legend = [figure_legend;temp_lega];
    %----------------------------------------------------------------------

    %----------------------Display The Results on Graph--------------------
    figure(h_error_PER_BFO);
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy_BFO),'b','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy2_BFO),'b-.','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
    hold on
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy3_BFO),'b-s','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy4_BFO),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy5_BFO),'b-.*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy6_BFO),'b-o','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy7_BFO),'b--*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
    plot(sort(processed_error_bits_BFO),sort(processed_PER_BFO),'b--','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     plot(sort(processed_error_bits_BFO_first_itr),sort(processed_PER_BFO_first_itr),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_BFO),sort(PER_tehroy3_BFO),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);

    title('PER for BFO')

    % colormap(summer);
    % colormap(jet);
    set(gca,'FontSize',16);    

%     figure(h_error_BER_BFO);
%     plot(sort(processed_error_bits_BFO),sort(BER_tehroy_BFO),'b','LineWidth',2,'Color',[tanh(ll/length(N)),ll/length(N),1-(ll/length(N))]);
%     hold on
%     plot(sort(processed_error_bits_BFO),sort(processed_BER_BFO),'b--','LineWidth',2,'Color',[tanh(ll/length(N)),ll/length(N),1-(ll/length(N))]);
%     title('BER for BFO')


%     figure(h_error_PER_WTA);
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy_WTA),'b','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy2_WTA),'b-.','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy3_WTA),'b-s','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy4_WTA),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy5_WTA),'b-.*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy6_WTA),'b-o','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy7_WTA),'b--*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(processed_PER_WTA),'b--','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(PER_tehroy3_WTA),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     title('PER for WTA')

%     figure(h_error_BER_WTA);
%     plot(sort(processed_error_bits_WTA),sort(BER_tehroy_WTA),'b','LineWidth',2,'Color',[tanh(ll/length(N)),ll/length(N),1-(ll/length(N))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(processed_BER_WTA),'b--','LineWidth',2,'Color',[tanh(ll/length(N)),ll/length(N),1-(ll/length(N))]);    
%     title('BER for WTA')
    %----------------------------------------------------------------------
    
end
%==========================================================================


% %----------------------Display Total Learn Iterations----------------------
% figure(h_itr);
% bar(hist_x_itr(1,:),hist_out_itr','grouped');
% % colormap(summer);
% colormap(jet);
%     
% title('Total Learning Iterations for Different Constraints','fontsize',16);
% xlabel('Total learning iterations','fontsize',16);
% ylabel('Percentage of converged trials','fontsize',16);
% set(gca,'FontSize',16);
% %--------------------------------------------------------------------------


% figure;
% plot(sort(processed_error_bits_WTA),sort(BER_tehroy5_BFO5),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%     plot(sort(processed_error_bits_WTA),sort(processed_BER_BFO),'b--','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%     hold on
%%
%======================ADD THE LEGEND TO THE FIGURES=======================
figure(h_error_PER_BFO)
legend(figure_legend)
% figure(h_error_PER_WTA)
% legend(figure_legend)
% figure(h_error_BER_BFO)
% legend(figure_legend)
% figure(h_error_BER_WTA)
% legend(figure_legend)
%==========================================================================