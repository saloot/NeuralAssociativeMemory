%==========================================================================
%******************FUNCTION: read_faulty_results***************************
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
% processed_error_bits: The list of processed number of initial erroneous nodes
% processed_PER: The list of processed Pattern Error Rates
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a faulty neural associative
% memory and reads the result of recall phase from the appropriate files. 
% The results will then be plotted and compared with theoretical values. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


% function [processed_error_bits,processed_PER] = read_faulty_results(N_in,K_in,L_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,pattern_neur_noise_range,const_neur_noise_range)
warning off




%%
%==============================INITIALIZATION==============================

%--------------------------Set Default Values------------------------------
if ~exist('err_bits','var')
    err_bits = 50;
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
    pattern_neur_noise_range = [0,.1,.3,.4,.6];
end

if ~exist('const_neur_noise_range','var')
    const_neur_noise_range = [0,.05,.2,.3];
end
%--------------------------------------------------------------------------


%---------------------------Simulation Parameters--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library')); 
addpath(genpath('./Verify_Theory')); 
hist_x_axis_itr = [0:20:1000];             % This is bin places for total learning iterations histogram
set(0,'defaulttextinterpreter','latex')
no_of_try = 1;
%--------------------------------------------------------------------------

%==========================================================================

%%
%============================PROCESS THE RESULTS===========================
for ll = 1:length(N_in)
    N = N_in(ll);
    K = K_in(ll);
    L = L_in(ll);
    load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           
    
    %--------------Calculate the Total Number of Constraints---------------
    m = zeros(1,length(hist_x_axis_itr));
    for l = 1:L
        for i = 1:index_max                
            fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(i),'.txt'], 'r');        
            if (fid >-1)
                learn_itr = fscanf(fid, '%d',[1,inf]);               
                fclose(fid);
            else
                learn_itr = 0;
            end
            
            [m_temp,~] = hist(learn_itr,hist_x_axis_itr);
            m = m+m_temp;
        end
    end  
    m = m/(l*i);
    %----------------------------------------------------------------------                       
                 
    
    for kji = 1:length(const_neur_noise_range)
        figure_legend = [];                         % This variable will be field with the legend for the graphs
        const_neur_noise = const_neur_noise_range(kji);        
%         h_error_PER = figure;
        eval(['h_error_BER_',num2str(kji),' = figure;']);
        for ihi = 1:length(pattern_neur_noise_range)
            pattern_neur_noise = pattern_neur_noise_range(ihi);
            
    
            %---------------------------Read Recall Results------------------------
            if (recall_algorithm_option == 0)
                fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'pattern_noise_',num2str(pattern_neur_noise),'const_noise_',num2str(const_neur_noise),'.txt'], 'r'); 
                if (fid > -1)
                    results = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
                    fclose(fid);
                else
                    error('Undefined WTA input file!')    
                end
            elseif (recall_algorithm_option == 2)
                fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'pattern_noise_',num2str(pattern_neur_noise),'const_noise_',num2str(const_neur_noise),'.txt'], 'r'); 
                if (fid > -1)
                    results = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
                    fclose(fid);
                else
                    error('Undefined BFS input file!')    
                end
            elseif (recall_algorithm_option == 1)
               fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_old_but_valid/clustered_neural_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'pattern_noise_',num2str(pattern_neur_noise),'const_noise_',num2str(const_neur_noise),'.txt'], 'r'); 
                if (fid > -1)
                    results = fscanf(fid, '%s %d %s %f %s %f',[10,inf]);
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
    
            %------------------Read Cluster Degree Distribution--------------------
            %     [deg_row,deg_column,lambda,rho] = read_cluster_degree(N,K,L,alpha0,beta0,theta0);    
            %     [deg_row_within,deg_column_within,lambda_within,rho_within,d_max,d_min,d_ave,m,pe_max,pe_max2,pe_ave,pe_ave2,pe_2_max,pe_2_ave] = read_cluster_degree_local(N,K,L,alpha0,beta0,theta0);
            %----------------------------------------------------------------------
    
            %---------------------------Compare with Theory------------------------
%             Pe_theory_BER = zeros(length(pattern_neur_noise_range),length(const_neur_noise_range),length(processed_error_bits));         
%             Pe_theory_BER= calculate_Pe_theoretical(N,K,L,alpha0,beta0,theta0,processed_error_bits,pattern_neur_noise,const_neur_noise,gamma_BFO,index);                            
%             [Pe_theory_BER,~] = calculate_Pe_theoretical_multiple_trial(N,K,L,alpha0,beta0,theta0,processed_error_bits,pattern_neur_noise,const_neur_noise,gamma_BFO,index,4);                            
            [Pe_theory_BER,~] = calculate_Pe_theoretical(N,K,L,alpha0,beta0,theta0,processed_error_bits,pattern_neur_noise,const_neur_noise,gamma_BFO,index,no_of_try);                                                            
%             [Pe_theory_BER_1,~] = calculate_Pe_theoretical_multiple_trial(N,K,L,alpha0,beta0,theta0,processed_error_bits,pattern_neur_noise,const_neur_noise,gamma_BFO,index,1);
            if (sum(Pe_theory_BER<0)>0)
                111;
            end
            %----------------------------------------------------------------------

            %---------------------Constructu the Legend---------------------------- 
            pp = pattern_neur_noise_range;
            pp(~pp) = inf;
            
            qq = const_neur_noise_range;
            qq(~qq) = inf;
            
            temp_leg = ['$\upsilon = ',num2str(pattern_neur_noise),', \nu = ',num2str(const_neur_noise),'$'];    
            for i = length(num2str(pattern_neur_noise)):length(num2str(min(pp)))
                temp_leg = [temp_leg,' '];
            end
            for i = length(num2str(const_neur_noise)):length(num2str(min(qq)))
                temp_leg = [temp_leg,' '];        
            end            
    
%             temp_lega = [temp_leg,' - Theory  '];
%             figure_legend = [figure_legend;temp_lega];
   
            temp_lega = [temp_leg,'-Thr'];
            figure_legend = [figure_legend;temp_lega];
%             temp_lega = [temp_leg,' - Theory Singl'];
%             figure_legend = [figure_legend;temp_lega];
            temp_lega = [temp_leg,'-Sim'];
            figure_legend = [figure_legend;temp_lega];
            %----------------------------------------------------------------------

%             %----------------------Display The Results on Graph--------------------
%             figure(h_error_PER);
%             plot(sort(processed_error_bits),sort(PER_tehroy_ave),'b','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%             hold on
%             plot(sort(processed_error_bits),sort(PER_tehroy_worst),'b-.','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%             hold on
%             plot(sort(processed_error_bits),sort(PER_tehroy_ave2),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%             hold on
%             plot(sort(processed_error_bits),sort(PER_tehroy_worst2),'b-s','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%             hold on    
%             plot(sort(processed_error_bits),sort(PER_tehroy_single),'b-*','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%             hold on
%             plot(sort(processed_error_bits),sort(PER_tehroy_multi),'b-.','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%             hold on
%             plot(sort(processed_error_bits),sort(processed_PER),'b--','LineWidth',2,'Color',[tanh(ll/length(N_in)),ll/length(N_in),1-(ll/length(N_in))]);
%             hold on

            %----------------------------------------------------------------------
            figure(eval(['h_error_BER_',num2str(kji)]));
            
            plot(sort(processed_error_bits)/N_tot,(sort(Pe_theory_BER)),'b-','LineWidth',6,'Color',[tanh(ihi/length(pattern_neur_noise_range)),ihi/length(pattern_neur_noise_range),1-(ihi/length(pattern_neur_noise_range))]);
            hold on
%             plot(sort(processed_error_bits)/N_tot,(sort(Pe_theory_BER_1)),'b-*','LineWidth',6,'Color',[tanh(ihi/length(pattern_neur_noise_range)),ihi/length(pattern_neur_noise_range),1-(ihi/length(pattern_neur_noise_range))]);
%             hold on
            plot(sort(processed_error_bits)/N_tot,(sort(processed_BER)),'b-.','LineWidth',6,'Color',[tanh(ihi/length(pattern_neur_noise_range)),ihi/length(pattern_neur_noise_range),1-(ihi/length(pattern_neur_noise_range))]);
            hold on
            
            title('BER')    
            
            eval(['figure(h_error_BER_',num2str(kji),');']);
            legend(figure_legend)
            set(gca,'FontSize',24)
        
            l = legend(figure_legend);
            set(l,'Interpreter','latex')
            xlhand = get(gca,'xlabel');
            ylhand = get(gca,'ylabel');
            set(xlhand,'string','$\epsilon$','fontsize',30)
            set(ylhand,'string','Final BER','fontsize',30)
        end
        
        
    end
end
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
% save(['/scratch/amir/Clustered_Neural/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Read_Results/final_results_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_alpha0_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS), '.mat'],'PER_tehroy_single','processed_error_bits',...
%     'processed_PER','processed_BER','lambda','rho','deg_column','deg_row','processed_count','m');
%==========================================================================