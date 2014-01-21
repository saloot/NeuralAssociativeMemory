%==========================================================================
%*********************FUNCTION: read_cluster_degree************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% cluster_size: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% no_clusters: The number of clusters if we have multiple levels (no_clusters = 1 for single level)
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% deg_row: The input cluster degrees from a cluster point of view in a clustered neural associative memory
% deg_column: The input pattern nodes degrees from a cluster point of view in a clustered neural associative memory
% lambda: The input pattern nodes degree distribution from a cluster point of view in a clustered neural associative memory
% rho: The input cluster degree distribution from a cluster point of view in a clustered neural associative memory
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function reads the result of the learning phase for a clustered
% neural associative memory and determines the degree distribution of the
% network from a clustered point of view, i.e. the case that all the
% constraint nodes in a cluster are contracted to one node.
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================


function [deg_row,deg_column,lambda,rho] = read_cluster_degree_v2(cluster_size,no_constraints,no_clusters,N_tot,alpha0,beta0,theta0,Q,learn_itr_max,source_depository)

%==============================INITIALIZATION==============================
shift_horiz  = round((N_tot - cluster_size)/(no_clusters-1));
deg_row = [0:cluster_size];
deg_column = [0:ceil(cluster_size/shift_horiz)];
rho = zeros(1,cluster_size+1);
lambda = zeros(1,1+ceil(cluster_size/shift_horiz));
%==========================================================================


%====================CALCULATE DEGREE DISTRIBUTIONS========================

%------------------------Read Weight Matrices------------------------------
s_row = zeros(1,no_clusters);
s_column = zeros(1,N_tot);
    
for l = 1:no_clusters
    file_name = [source_depository,'/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(no_constraints),...
        '_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];
    fid = fopen(file_name,'r');
    if (fid == -1)
        file_name = [source_depository,'/Partial_Convergence/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(no_constraints),...
            '_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];
        fid = fopen(file_name,'r');
    end
    if (fid == -1)
        continue;
    else
        fclose(fid);
    end
    
    load(file_name)
    W = W_total(:,2:end);
    
    for iik = 1:cluster_size            
        if (norm(W(:,iik))>0)               
            s_column((l-1)*shift_horiz + iik) = s_column((l-1)*shift_horiz + iik) + 1;                
            s_row(l) = s_row(l) + 1;
           
        end
    end
end
    
    [rho_temp,deg_row] = hist(s_row,deg_row);
    rho = rho + rho_temp;
    [lambda_temp,deg_column] = hist(s_column,deg_column);       
    lambda = lambda+lambda_temp;
    

%--------------------------------------------------------------------------

%--------------Normalize Row Degree Distribution of Cluster----------------
rho = rho/sum(rho);
%--------------------------------------------------------------------------


%-----------------Normalize Column Degree Distribution---------------------
lambda = lambda/sum(lambda);
%--------------------------------------------------------------------------

%==========================================================================

% %---------------------Display Row Degree Distribution----------------------
% subplot(2,1,1)
% bar(deg_row,rho);
% xlim([0 cluster_size]);
% xlabel('Degree','fontsize',30);
% ylabel('Row degree distribution','fontsize',30);
% set(gca,'FontSize',24);
% %--------------------------------------------------------------------------
% 
% %--------------------Display Column Degree Distribution--------------------
% subplot(2,1,2)
% bar(deg_column,lambda);
% xlim([0 ceil(cluster_size/shift_horiz)]);
% xlabel('Degree','fontsize',30);
% ylabel('Column degree distribution','fontsize',30);
% set(gca,'FontSize',24);
% %--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
