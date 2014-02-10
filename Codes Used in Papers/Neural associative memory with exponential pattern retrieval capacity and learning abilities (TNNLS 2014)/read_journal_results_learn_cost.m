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
% alpha0= 0.95;
% beta0 = .85;
% theta0 = 0.15;


% N_in = 800;
% K_in = 400;
% alpha0= 0.75;
% beta0 = 1;
% theta0 = 0.036;
% gamma_BFO = 1;
% gamma_BFS = 0.999;

N_in = 400;
K_in = 200;
alpha0= 0.95;
beta0 = 1;
theta0 = 0.021;


gamma_BFO = 0.99;
gamma_BFS = 0.99;

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

total_cost_avg = [];                        % This is the total learning cost per iteration over the whole ensemble
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
    
   
    
    %------------------Display Total Learn Iterations----------------------
        
    
    %-----Read Total Learn Iterations From the File------    
    fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(1),'.txt'], 'r');
    learn_itr = fscanf(fid, '%d',[1,inf]);        
    no_local_constraints = length(learn_itr);                
    fclose(fid);    
    %----------------------------------------------------
    
    
    %-------------Construct the Output-------------------
    [m,xout] = hist(learn_itr,hist_x_axis_itr);
    hist_x_itr = [hist_x_itr;xout];
    hist_out_itr = [hist_out_itr;m/no_local_constraints];        
    %----------------------------------------------------  
    
    %----------------------------------------------------------------------
    
    
    %--------------------Display the Sparsity Degree-----------------------                        
    
    [deg_row,deg_column,lambda,rho] = deg_dribution(N,K,alpha0,beta0,theta0,index_max);    
            
    %----------------------------------------------------------------------
    
    %-------------------Display Sample Learning Costs---------------------- 
    
    %-----------Select the approprate figure-------------
    if (ll ==1)
        h_cost = figure;        
    else        
        figure(h_cost);        
    end    
    %----------------------------------------------------    
    
    for ind = 1:index_max-1
        %--------Read Iteration Costs From the File----------            
        fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Learn_cost_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(ind),'.txt'], 'r');    
        learn_costs = fscanf(fid, '%f');               
        fclose(fid);    
        %----------------------------------------------------
    
        %----------------Process the costs----------------    
        total_cost = [];
        marker = 0;
        no_consts = 0;
        for i = 1:length(learn_costs)-1
            if ( (learn_costs(i+1) > 0) && (learn_costs(i) < 0))
                temp = exp(learn_costs(marker+1:i))';
                marker = i;
                no_consts = no_consts + 1;
                l1 = length(temp);
                l2 = length(total_cost);
                if (l1 < l2)
                    temp = [temp,zeros(1,l2-l1)];
                    total_cost = total_cost + temp;
                else
                    total_cost = [total_cost,zeros(1,l1-l2)];
                    total_cost = total_cost + temp;
                end
            end
        end
        l1 = length(total_cost_avg);                
        l2 = length(total_cost);                
        if (l1 < l2)                    
            total_cost_avg = [total_cost_avg,zeros(1,l2-l1)];
            total_cost_avg = total_cost + total_cost_avg;                
        else            
            total_cost = [total_cost,zeros(1,l1-l2)];            
            total_cost_avg = total_cost + total_cost_avg;                
        end
    end    
    total_cost_avg = total_cost_avg/index_max;
        
    %----------------------------------------------------    
    
    %------------------Plot the Results------------------    
    plot((total_cost_avg),'b--','LineWidth',2,'Color',[tanh(ll/length(N)),ll/length(N),1-(ll/length(N))]);
    hold on    
        
    title('MSE as a Function of Learning Iteration','fontsize',16);
    xlabel('Learning iterations','fontsize',16);
    ylabel('MSE','fontsize',16);
    set(gca,'FontSize',16);
    %----------------------------------------------------
    
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

%----------------------Display Sparsity Percentage-------------------------
figure(h_sparse);
bar(deg_row/(N),rho','grouped');
% colormap(summer);
colormap(jet);
    
title('Total Sparsity Percentage for Different Variable Nodes','fontsize',16);
xlabel('Total Sparsity Percentage','fontsize',16);
ylabel('Percentage of trials with the specified degree of sparsity','fontsize',16);
set(gca,'FontSize',16);
%--------------------------------------------------------------------------

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