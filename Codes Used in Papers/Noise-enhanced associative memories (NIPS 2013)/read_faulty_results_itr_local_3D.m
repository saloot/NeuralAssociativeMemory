%==========================================================================
%***************FUNCTION: read_faulty_results_itr_3D***********************
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
% memory and reads the result of recall phase from the appropriate files. 
% The results will then be plotted and compared with theoretical values in
% a 3D manner.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================



% Fig 1-3: 11
% Fig 4-6: 6
function [processed_error_bits,processed_recall_itr_local,fraction_unsuccess_count_local] = read_faulty_results_itr_local_3D(N_in,K_in,L_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,pattern_neur_noise_range,const_neur_noise_range)
warning off


%%
%==============================INITIALIZATION==============================

%--------------------------Set Default Values------------------------------
if ~exist('err_bits','var')
    err_bits = 11;
end

if ~exist('N_in','var')
    N_in = 40;
end

if ~exist('K_in','var')
    K_in = 20;
end

if ~exist('L_in','var')
    L_in = 50;
end

if ~exist('alpha0','var')
    alpha0 = 0.95;
end

if ~exist('beta0','var')
    beta0 = 0.75;
end

if ~exist('theta0','var')
    theta0 = 0.05;
end

if ~exist('gamma_BFO','var')
    gamma_BFO = 0.95;
end

if ~exist('gamma_BFS','var')
    gamma_BFS = 0.95;
end

if ~exist('recall_algorithm_option','var')
    recall_algorithm_option = 1;
end

if ~exist('pattern_neur_noise_range','var')
    pattern_neur_noise_range = [0:.1:.5];
end

if ~exist('const_neur_noise_range','var')
    const_neur_noise_range = [0:.1:.5];
end
%--------------------------------------------------------------------------


%---------------------------Simulation Parameters--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library')); 
figure_legend = [];                         % This variable will be field with the legend for the graphs
hist_x_axis_itr = [0:20:1000];             % This is bin places for total learning iterations histogram
PER_processed = zeros(length(pattern_neur_noise_range),length(const_neur_noise_range));
BER_processed = zeros(length(pattern_neur_noise_range),length(const_neur_noise_range));
set(0,'defaulttextinterpreter','latex')
%--------------------------------------------------------------------------
    
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
        for ihi = 1:length(pattern_neur_noise_range)
            pattern_neur_noise = pattern_neur_noise_range(ihi);
        
                        
%         %---------------------------Read Recall Results------------------------
%         if (recall_algorithm_option == 0)                
%             fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/recall_itr_local_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'pattern_noise_',num2str(pattern_neur_noise),'const_noise_',num2str(const_neur_noise),'.txt'], 'r'); 
%             if (fid > -1)                    
%                 results_itr = fscanf(fid, '%s %d %s %f %s %f',[1000,inf]);
%                 fclose(fid);                
%             else                
%                 error('Undefined WTA input file!')
%             end            
%         elseif (recall_algorithm_option == 2)
%             fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/recall_itr_local_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'pattern_noise_',num2str(pattern_neur_noise),'const_noise_',num2str(const_neur_noise),'.txt'], 'r'); 
%             if (fid > -1)                    
%                 results_itr = fscanf(fid, '%s %d %s %f %s %f',[1000,inf]);
%                 fclose(fid);                
%             else                
%                 error('Undefined BFS input file!')
%             end            
%         elseif (recall_algorithm_option == 1)                              
%             fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/recall_itr_local_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'pattern_noise_',num2str(pattern_neur_noise),'const_noise_',num2str(const_neur_noise),'.txt'], 'r'); 
%             if (fid > -1)                                                                               
%                 results_itr = fscanf(fid, '%s',[1000,inf]);
%                 fclose(fid);                
%             else                
%                 error('Undefined BFO recall itr input file!')
%             end            
%         else            
%             error('Unknown recall algorithm');            
%         end        
%         %----------------------------------------------------------------------
               
        %---------------------Process the Results------------------------------
        [n_rows,n_columns] = size(results_itr);
%         index_values = [];
%         unprocessed_error_bits = [];
%         unprocessed_recall_itr= [];
%         irr = 1;
%         ill = 1;
%         countr = 0;
        
        while ill <=n_columns
            while irr <= n_rows 
                countr = countr + 1;
                if (countr > 1e8)
                    111;
                    counter = 0;
                end
                if (results_itr(irr,ill) == 'e')
                    index_values = [index_values,irr,ill];                    
                    irr = irr + 1;
                    if irr > n_rows
                        ill = ill + 1;
                        irr = mod(irr,n_rows);
                        if (ill > n_columns)
                            break;
                        end
                    end
                     
                    error_str = [results_itr(irr,ill)];
                    irr = irr + 1;
                    if irr > n_rows
                        ill = ill + 1;
                        irr = mod(irr,n_rows);
                        if (ill > n_columns)
                            break;
                        end
                    end
                    error_str = [error_str,results_itr(irr,ill)];
                    unprocessed_error_bits = [unprocessed_error_bits, str2num(error_str)];
                    
                    irr = irr + 4;
                    if irr > n_rows
                        ill = ill + 1;
                        irr = mod(irr,n_rows);
                        if (ill > n_columns)
                            break;
                        end
                    end
                    itr_str = [results_itr(irr,ill)];
                    for ikkl = 1:4
                        irr = irr + 1;
                        if irr > n_rows
                            ill = ill + 1;
                            irr = mod(irr,n_rows);
                            if (ill > n_columns)
                                break;
                            end
                        end
                        if ( (double(results_itr(irr,ill)) >= 48 ) && ( double(results_itr(irr,ill)) <= 57) )
                            itr_str = [itr_str,results_itr(irr,ill)];
                        else
                            break;
                        end
                    end                        
                    
                    unprocessed_recall_itr = [unprocessed_recall_itr,str2num(itr_str)];
                else
                    irr = irr + 1;
                    if irr > n_rows
                        ill = ill + 1;
                        irr = mod(irr,n_rows);
                        if (ill > n_columns)
                            break;
                        end
                    end
                end
            end 
            if (mod(irr,10) == 1)
                111;
            end
            
            if (ill > n_columns)            
                break;                        
            end
            
        end
            
            
%             processed_error_bits = [];
            processed_recall_itr_local = 0;            
            processed_unsuccess_count = 0;
            processed_count = 0;
            for i = 1:length(unprocessed_error_bits)
                if (unprocessed_error_bits(i) == err_bits)    
                    processed_error_bits = err_bits;                                                                                                        
                    
                    processed_count = processed_count + 1;
                    if (unprocessed_recall_itr(i)>=1000)
                        processed_unsuccess_count = processed_unsuccess_count + 1;
                        processed_recall_itr_local = processed_recall_itr_local + 40;
                    else
                        processed_recall_itr_local = processed_recall_itr_local + unprocessed_recall_itr(i);
                    end
                end
            end
            
%             processed_recall_itr = processed_recall_itr - 1000*processed_unsuccess_count;            
            processed_recall_itr_local = processed_recall_itr_local./(processed_count);%-processed_unsuccess_count);
            %----------------------------------------------------------------------                                              
            
        %------------------Store the Processed Results----------------------------                               
        BER_processed(ihi,kji) = processed_recall_itr_local;        
        if (isnan(processed_recall_itr_local))
            111
        end
        %----------------------------------------------------------------------                       

        end
    end
end
%==========================================================================

%==========================PLOT THE RESULTS================================
figure
surf(pattern_neur_noise_range,const_neur_noise_range,BER_processed');
xlabel('$\upsilon$');
ylabel('$\nu$');
zlabel('Final BER');
set(gca,'FontSize',24)
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
zlhand = get(gca,'zlabel');
set(xlhand,'string','$\upsilon$','fontsize',30)
set(ylhand,'string','$\nu$','fontsize',30)
set(zlhand,'string','Final BER','fontsize',30)

figure
figure_legend = [];
pp = pattern_neur_noise_range;
pp(~pp) = inf;
    
for i = 1:length(pattern_neur_noise_range)
    plot(const_neur_noise_range,BER_processed(i,:),'b-','LineWidth',6,'Color',[tanh(i/length(pattern_neur_noise_range)),i/length(pattern_neur_noise_range),1-(i/length(pattern_neur_noise_range))]);
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
ylabel('Final BER');
set(gca,'FontSize',24)
title('BER as a function of constraint neurons noise for various pattern neurons noise parameter')
l = legend(figure_legend);
set(l,'Interpreter','latex')
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','$\nu$','fontsize',30)
set(ylhand,'string','Final BER','fontsize',30)

figure
figure_legend = [];
pp = const_neur_noise_range;
pp(~pp) = inf;
    
for i = 1:length(const_neur_noise_range)
    plot(pattern_neur_noise_range,BER_processed(:,i)','b-','LineWidth',6,'Color',[tanh(i/length(const_neur_noise_range)),i/length(const_neur_noise_range),1-(i/length(const_neur_noise_range))]);
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
ylabel('Final BER');
set(gca,'FontSize',24)
title('BER as a function of pattern neurons noise for various constraint neurons noise parameter')
l = legend(figure_legend);
set(l,'Interpreter','latex')
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');
set(xlhand,'string','$\upsilon$','fontsize',30)
set(ylhand,'string','Final BER','fontsize',30)
%==========================================================================

%==========================================================================

%======================ADD THE LEGEND TO THE FIGURES=======================

% figure(h_error_PER)
% legend(figure_legend)


111;
% figure(h_error_BER_BFO)
% legend(figure_legend)
% figure(h_error_BER_WTA)
% legend(figure_legend)
%==========================================================================

%======================SAVE THE PROCESSED RESULTS==========================
% mkdir(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)],'Read_Results');        % Create a specific folder for the current N and K
% save(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Read_Results/final_results_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_alpha0_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '.mat'],'PER_tehroy_single','processed_pattern_noise',...
%     'processed_PER','processed_BER','lambda','rho','deg_column','deg_row','processed_count','m');
%==========================================================================