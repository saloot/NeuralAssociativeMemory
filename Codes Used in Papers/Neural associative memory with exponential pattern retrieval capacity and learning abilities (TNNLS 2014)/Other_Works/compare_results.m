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

% N_in = 200;
% K_in = 100;
% alpha0= 0.45;
% beta0 = .45;
% theta0 = 0.015;
% gamma_BFO = .95;
% gamma_BFS = 0.95;

N_in = 800;
K_in = 400;
alpha0= 0.95;
beta0 = 1;
theta0 = 0.029;
% gamma_BFO = 1;
% gamma_BFS = 0.999;

gamma_BFO = 0.99;
gamma_BFS = 0.99;

set(0,'defaulttextinterpreter','latex')
no_of_patterns_range = [50,100,200];
% function [processed_error_bits_BFO,processed_PER_WTA,processed_PER_BFO,processed_PER_BFS] = read_Jankowski_results(N_in,K_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS)
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
        error('Undefined WTA input file!')
    
    end
    %----------------------------------------------------------------------

        
    %----------------Process Journal Recall Results------------------------
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
    %----------------------------------------------------------------------                       

    
    %----------------Process Jankwoski's Recall Results------------------------
    final_processed_error_bits_Jankowski_random = [];
    final_processed_error_bits_Jankowski_non_random = [];
    
    final_processed_PER_Jankowski_random = [];
    final_processed_PER_Jankowski_non_random = [];
    
    final_processed_BER_Jankowski_random = [];
    final_processed_BER_Jankowski_non_random = [];
    
    final_processed_count_Jankowski_random = [];
    final_processed_count_Jankowski_non_random = [];
    
    counter = 0;
    for random_flag = 0:1
        for ii = 1:length(no_of_patterns_range)
            no_of_patterns = no_of_patterns_range(ii);
            counter = counter + 1;
            processed_error_bits_Jankowski_random = [];
            processed_error_bits_Jankowski_non_random = [];
    
            processed_PER_Jankowski_random = [];
            processed_PER_Jankowski_non_random = [];
    
            processed_BER_Jankowski_random = [];
            processed_BER_Jankowski_non_random = [];
    
            processed_count_Jankowski_random = [];
            processed_count_Jankowski_non_random = [];
    
            if (random_flag)
                fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/Jankowski/N_',num2str(N),...
                             '_Random/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'r');
            else
                fid = fopen(['/scratch/amir/ITW_Journal/Recall_Results/Jankowski/N_',num2str(N),'_K_',num2str(K),...
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
            
            processed_BER_Jankowski_random = processed_BER_Jankowski_random./processed_count_Jankowski_random;                        
            processed_PER_Jankowski_random = processed_PER_Jankowski_random./processed_count_Jankowski_random;
            
            processed_BER_Jankowski_non_random = processed_BER_Jankowski_non_random./processed_count_Jankowski_non_random;                        
            processed_PER_Jankowski_non_random = processed_PER_Jankowski_non_random./processed_count_Jankowski_non_random;
            
            final_processed_error_bits_Jankowski_random = [final_processed_error_bits_Jankowski_random;processed_error_bits_Jankowski_random];
            final_processed_error_bits_Jankowski_non_random = [final_processed_error_bits_Jankowski_non_random;processed_error_bits_Jankowski_non_random];
    
            final_processed_PER_Jankowski_random = [final_processed_PER_Jankowski_random;processed_PER_Jankowski_random];
            final_processed_PER_Jankowski_non_random = [final_processed_PER_Jankowski_non_random;processed_PER_Jankowski_non_random];
    
            final_processed_BER_Jankowski_random = [final_processed_BER_Jankowski_random;processed_BER_Jankowski_random];
            final_processed_BER_Jankowski_non_random = [final_processed_BER_Jankowski_non_random;processed_BER_Jankowski_non_random];
    
            final_processed_count_Jankowski_random = [final_processed_count_Jankowski_random;processed_count_Jankowski_random];
            final_processed_count_Jankowski_non_random = [final_processed_count_Jankowski_non_random;processed_count_Jankowski_non_random];
            
        end
    end
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


    %---------------------Constructu the Legend----------------------------
    temp_leg = ['$n=',num2str(N),', k=',num2str(K),'$'];    
    for i = length(num2str(N)):length(num2str(max(N_in)))
        temp_leg = [temp_leg,' '];
    end
    for i = length(num2str(K)):length(num2str(max(K_in)))
        temp_leg = [temp_leg,' '];        
    end            
    
    temp_lega = [temp_leg,'- Our Method  '];
    
    figure_legend = [figure_legend;temp_lega];
    
    for ii = 1:length(no_of_patterns_range)
        no_of_patterns = no_of_patterns_range(ii);
        
        temp_lega = [temp_leg,'- [3], $C=', num2str(no_of_patterns),'$'];
        for i = length(num2str(no_of_patterns)):length(num2str(max(no_of_patterns_range)))-1
            temp_lega = [temp_lega,' '];        
        end
        figure_legend = [figure_legend;temp_lega];
        
        
        
        temp_lega = ['n=',num2str(N),', Random - [3], $C=', num2str(no_of_patterns),'$  '];        
        for i = length(num2str(no_of_patterns)):length(num2str(max(no_of_patterns_range)))-1
            temp_lega = [temp_lega,' '];        
        end
        figure_legend = [figure_legend;temp_lega];
    end
    
    %----------------------------------------------------------------------

    %----------------------Display The Results on Graph--------------------    
    figure(h_error_PER_BFO);      
    plot(sort(processed_error_bits_BFO),sort(processed_BER_BFO),'b-','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);    
    hold on
    for ii = 1:length(no_of_patterns_range)        
        plot(sort(final_processed_error_bits_Jankowski_non_random(ii,:)),sort(final_processed_PER_Jankowski_non_random(ii,:)),'b-.','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);    
        plot(sort(final_processed_error_bits_Jankowski_random(ii,:)),sort(final_processed_PER_Jankowski_random(ii,:)),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);    
    end
    

    set(gca,'FontSize',16);    
    %----------------------------------------------------------------------
    
end
%==========================================================================



%======================ADD THE LEGEND TO THE FIGURES=======================
figure(h_error_PER_BFO)
l = legend(figure_legend);
set(l,'Interpreter','latex')           
xlhand = get(gca,'xlabel');            
ylhand = get(gca,'ylabel');            
set(xlhand,'string','$\epsilon$','fontsize',30)            
set(ylhand,'string','Final BER','fontsize',30)
%==========================================================================