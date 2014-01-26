%==========================================================================
%*************FUNCTION: read_faulty_results_noiseless**********************
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
% recall_algorithm_option: The recall algorithm identifier (0 for winner-take-all, 1 for the original bit flipping and 2 for the simplified bit flipping
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% const_neur_noise: The maximum amount of noise a constraint neuron will "suffer" from
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% processed_pattern_noise: The list of processed number of initial erroneous nodes
% processed_PER: The list of processed Pattern Error Rates
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a faulty neural associative
% memory and reads the result of recall phase when there are no external
% noise in the network.
% The results will then be plotted and compared with theoretical values. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function [processed_pattern_noise,processed_recall_itr,fraction_unsuccess_count] = read_faulty_results_noiseless(N_in,K_in,L_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,pattern_neur_noise_range,const_neur_noise_range)
warning off

%%
%==============================INITIALIZATION==============================
addpath(genpath('/home1/amir/cluster/Common_Library')); 
figure_legend = [];                         % This variable will be field with the legend for the graphs
hist_x_axis_itr = [0:20:1000];             % This is bin places for total learning iterations histogram
PER_processed = zeros(length(pattern_neur_noise_range),length(const_neur_noise_range));
BER_processed = zeros(length(pattern_neur_noise_range),length(const_neur_noise_range));
set(0,'defaulttextinterpreter','latex')
%==========================================================================

%%
%============================PROCESS THE RESULTS===========================
for ll = 1:length(N_in)
    N = N_in(ll);
    K = K_in(ll);
    L = L_in(ll);
    load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           
                       
                 
    for kji = 1:length(pattern_neur_noise_range)
        const_neur_noise = const_neur_noise_range(kji);        
       
                        
        %---------------------------Read Recall Results------------------------
        if (recall_algorithm_option == 0)                            
            fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'const_noise_',num2str(const_neur_noise),'_noiseless.txt'], 'a+');  
            if (fid > -1)                    
                results = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);            
                fclose(fid);                
            else                
                error('Undefined WTA input file!')
            end            
        elseif (recall_algorithm_option == 2)
            fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'const_noise_',num2str(const_neur_noise),'_noiseless.txt'], 'a+');  
            if (fid > -1)                    
                results = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);                    
                fclose(fid);                
            else                
                error('Undefined BFS input file!')
            end            
        elseif (recall_algorithm_option == 1)                              
            fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'const_noise_',num2str(const_neur_noise),'_noiseless.txt'], 'a+');  
            if (fid > -1)                                                                               
                results = fscanf(fid, '%s %f %s %f %s %f',[10,inf]);                    
                fclose(fid);                
            else                
                error('Undefined BFO recall itr input file!')
            end            
        else            
            error('Unknown recall algorithm');            
        end        
        %----------------------------------------------------------------------
               
        %---------------------Process the Results------------------------------            
        unprocessed_pattern_noise = results(2,:);            
        unprocessed_PER = results(6,:);            
        unprocessed_BER = results(10,:);
                
        processed_pattern_noise = [0];   
        processed_PER = [0];            
        processed_BER = [0];            
        processed_count = [1];            
        for i = 1:length(unprocessed_pattern_noise)                
            processed_flag = 0;                
            for j = 1:length(processed_pattern_noise)                    
                if (unprocessed_pattern_noise(i) == processed_pattern_noise(j))                        
                    processed_flag = 1;                        
                    break;                    
                end                
            end
                        
            if (processed_flag == 0)
                processed_pattern_noise = [processed_pattern_noise,unprocessed_pattern_noise(i)];                    
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
            
        %------------------Store the Processed Results----------------------------                       
        if (norm(sort(processed_pattern_noise(1:5))-pattern_neur_noise_range(1:5))<1e-6)
            BER_processed(:,kji) = sort_read_results(processed_pattern_noise(1:5),processed_BER(1:5));
            PER_processed(:,kji) = sort_read_results(processed_pattern_noise(1:5),processed_PER(1:5));
        else
            error('Something is wront with the ranges!')
        end
        %----------------------------------------------------------------------                       

    end
end
%==========================================================================

%==========================PLOT THE RESULTS================================
surf(pattern_neur_noise_range,const_neur_noise_range,PER_processed');
xlabel('$\upsilon$');
ylabel('$\nu$');
zlabel('PER');
set(gca,'FontSize',24)
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
zlhand = get(gca,'zlabel');
set(xlhand,'string','$\upsilon$','fontsize',30)
set(ylhand,'string','$\nu$','fontsize',30)
set(zlhand,'string','PER','fontsize',30)

figure
figure_legend = [];
pp = pattern_neur_noise_range;
pp(~pp) = inf;
    
for i = 1:length(pattern_neur_noise_range)
    plot(const_neur_noise_range,PER_processed(i,:),'b-','LineWidth',6,'Color',[tanh(i/length(pattern_neur_noise_range)),i/length(pattern_neur_noise_range),1-(i/length(pattern_neur_noise_range))]);
    hold on
    
    %---------------------Constructu the Legend----------------------------     
    temp_leg = ['$\upsilon = ',num2str(pattern_neur_noise_range(i)),'$'];    
    for iik = length(num2str(pattern_neur_noise_range(i))):length(num2str(min(pp)))
        temp_leg = [temp_leg,' '];
    end
    figure_legend = [figure_legend;temp_leg];
    %----------------------------------------------------------------------
end
xlabel('$\nu$');    
ylabel('PER');
set(gca,'FontSize',24)
title('PER as a function of constraint neurons noise for various pattern neurons noise parameter')
l = legend(figure_legend);
set(l,'Interpreter','latex')
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','$\nu$','fontsize',30)
set(ylhand,'string','PER','fontsize',30)

figure
figure_legend = [];
pp = const_neur_noise_range;
pp(~pp) = inf;
    
for i = 1:length(const_neur_noise_range)
    plot(pattern_neur_noise_range,PER_processed(:,i)','b-','LineWidth',6,'Color',[tanh(i/length(const_neur_noise_range)),i/length(const_neur_noise_range),1-(i/length(const_neur_noise_range))]);
    hold on
    
    %---------------------Constructu the Legend----------------------------     
    temp_leg = ['$\nu = ',num2str(const_neur_noise_range(i)),'$'];    
    for iik = length(num2str(const_neur_noise_range(i))):length(num2str(min(pp)))
        temp_leg = [temp_leg,' '];
    end
    figure_legend = [figure_legend;temp_leg];
    %----------------------------------------------------------------------
end
xlabel('$\upsilon$');    
ylabel('PER');
set(gca,'FontSize',24)
title('PER as a function of pattern neurons noise for various constyraint neurons noise parameter')
l = legend(figure_legend);
set(l,'Interpreter','latex')
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','$\upsilon$','fontsize',30)
set(ylhand,'string','PER','fontsize',30)
%==========================================================================


%======================SAVE THE PROCESSED RESULTS==========================
% mkdir(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)],'Read_Results');        % Create a specific folder for the current N and K
% save(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Read_Results/final_results_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_alpha0_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '.mat'],'PER_tehroy_single','processed_pattern_noise',...
%     'processed_PER','processed_BER','lambda','rho','deg_column','deg_row','processed_count','m');
%==========================================================================