%==========================================================================
%******************FUNCTION: read_cluster_results**************************
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
% try_max: The maximum number of iterations in outside-cluster recall algorithm
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% processed_error_bits: The list of processed number of initial erroneous nodes
% processed_PER: The list of processed Pattern Error Rates
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a clustered neural associative
% memory and reads the result of recall phase from the appropriate files. 
% The results will then be plotted and compared with theoretical values. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


% function [processed_error_bits,processed_PER] = read_cluster_results_v2(N_in,K_in,L_in,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,try_max)
% warning off

%%
%==============================INITIALIZATION==============================
addpath(genpath('/home1/amir/cluster/Common_Library')); 
figure_legend = [];                         % This variable will be field with the legend for the graphs
hist_x_axis_itr = [0:20:1000];             % This is bin places for total learning iterations histogram
h_error_PER = figure;
h_error_BER = figure;
%==========================================================================

%%
%============================PROCESS THE RESULTS===========================
simulation_set = 14;
db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/CIFAR_10_Gray_Mixed_Whitened_DB_Q_15_Binary_Projected.mat';       
db_name_in = 'CIFAR_10_Gray_Mixed_Whitened_DB_Q_15_Binary_Projected';
load(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set)]);

current_params = simulation_parameters;

eval(['cluster_size=',num2str(current_params(1)),';']);    
eval(['no_clusters=',num2str(current_params(2)),';']);
eval(['const_learn=',num2str(current_params(3)),';']);    
eval(['alpha0=',num2str(current_params(4)),';']);
eval(['beta0=',num2str(current_params(5)),';']);
eval(['theta0=',num2str(current_params(6)),';']);
eval(['Q=',num2str(current_params(7)),';']);
eval(['learn_itr_max=',num2str(current_params(8)),';']);


N_tot = 4096;
try_max = 400;
Q = 1;



%--------------------Create the Sub-folder If Necessary--------------------
slash_flag = 0;
for i = length(db_file_in):-1:1
    if (strcmp(db_file_in(i),'/'))
        if (slash_flag == 0)
            break;
        else
            slash_flag = slash_flag+1;
        end
    end
end


destination_folder = [db_file_in(1:i-1),'/Recall_Results/Clustered_Version'];
%--------------------------------------------------------------------------

                        
%-----------------------------Read Recall Results--------------------------
fid = fopen([destination_folder,'/n_',num2str(cluster_size),'_consts_',num2str(const_learn),...
    '_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_Q_',num2str(Q),...
    '_itr_',num2str(learn_itr_max),'_',num2str(try_max),'.txt'], 'r');     
if (fid > -1)            
    results = fscanf(fid, '%s %d %s %f %s %f %s %f',[15,inf]);            
    fclose(fid);       
else    
    error('Undefined WTA input file!')   
end
%--------------------------------------------------------------------------        

%-----------------------Process the Results--------------------------------    
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
       
% %------------------Read Cluster Degree Distribution--------------------
[deg_row,deg_column,lambda,rho] = read_cluster_degree_v2(cluster_size,const_learn,no_clusters,N_tot,alpha0,beta0,theta0,Q,learn_itr_max,source_depository);
[deg_row_within,deg_column_within,lambda_within,rho_within,d_max,d_min,d_ave,m,pe_max,pe_max2,pe_ave,pe_ave2] = read_cluster_degree_local_v2(cluster_size,const_learn,no_clusters,alpha0,beta0,theta0,Q,learn_itr_max,source_depository);   
% %----------------------------------------------------------------------
    
%---------------------------Compare with Theory------------------------
% pe_max = 1-(1-pe_max)^(cluster_size-1);
[PER_tehroy_worst,PER_tehroy_ave,BER_tehroy_worst,BER_tehroy_ave] = P_e_clustered_theory_main_v2(cluster_size,N_tot,lambda,rho,processed_error_bits,try_max,pe_max,pe_ave);
% [PER_tehroy_worst,PER_tehroy_ave,BER_tehroy_worst,BER_tehroy_ave] = P_e_clustered_theory_main_v2(cluster_size,N_tot,lambda,rho,processed_error_bits,try_max,pe_max,pe_ave);
%----------------------------------------------------------------------

%-----------------------Constructu the Legend------------------------------
temp_leg = ['$n=',num2str(cluster_size),', m=',num2str(const_learn),', t_{\max}=',num2str(try_max),', T_{\max}=',num2str(learn_itr_max),'$-Thr'];    
figure_legend = [figure_legend;temp_leg];    
temp_leg = ['$n=',num2str(cluster_size),', m=',num2str(const_learn),', t_{\max}=',num2str(try_max),', T_{\max}=',num2str(learn_itr_max),'$-Sim'];    
figure_legend = [figure_legend;temp_leg];    
%--------------------------------------------------------------------------
    
%----------------------Display The Results on Graph--------------------   
figure(h_error_PER);
% plot(sort(processed_error_bits),sort(PER_tehroy_ave),'b','LineWidth',2);    
plot(sort(processed_error_bits),sort(PER_tehroy_worst),'b','LineWidth',2);    
hold on

plot(sort(processed_error_bits),sort(processed_PER),'b--','LineWidth',2);
hold on
title('PER')    
%----------------------------------------------------------------------

figure(h_error_BER);

% plot(sort(processed_error_bits),sort(BER_tehroy_ave),'b','LineWidth',2);
plot(sort(processed_error_bits),sort(BER_tehroy_worst),'b','LineWidth',2);
hold on

plot(sort(processed_error_bits),sort(processed_BER),'b--','LineWidth',2);
hold on
plot(sort(processed_error_bits),sort(processed_error_bits)/N_tot);
title('BER')    
%==========================================================================

%======================ADD THE LEGEND TO THE FIGURES=======================
figure(h_error_PER)
legend(figure_legend)
figure(h_error_BER)
legend(figure_legend)
%==========================================================================

