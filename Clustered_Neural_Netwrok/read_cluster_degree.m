%==========================================================================
%*********************FUNCTION: read_cluster_degree************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
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


function [deg_row,deg_column,lambda,rho] = read_cluster_degree(N,K,L,alpha0,beta0,theta0)

%==============================INITIALIZATION==============================
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           

deg_column = [1:L];
deg_row = [1:N*L/deg_cluster_min];
lambda = zeros(1,L);
rho = zeros(1,N*L/deg_cluster_min);
N_tot = N*L/deg_cluster_min;
%==========================================================================


%====================CALCULATE DEGREE DISTRIBUTIONS========================

%------------------------Read Weight Matrices------------------------------
for ind = 1:index_max
    s_row = zeros(1,L);
    s_column = zeros(1,N_tot);
    for l = 1:L
        
        load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);
            n = length(index_l);
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(ind),'.txt'], 'r');
        if (fid > -1)
        W = fscanf(fid, '%f',[n,inf]);
        W = W';
        [m,n] = size(W);
        fclose(fid);
        
        for iik = 1:n
            if (norm(W(:,iik))>0)
                s_column(index_l(iik)) = s_column(index_l(iik)) + 1;
                s_row(l) = s_row(l) + 1;
            end
        end
                
         end
    end
    
    [rho_temp,deg_row] = hist(s_row,deg_row);
    rho = rho + rho_temp;
    [lambda_temp,deg_column] = hist(s_column,deg_column);       
    lambda = lambda+lambda_temp;
    
end
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
% xlim([0 N_tot]);
% xlabel('Degree','fontsize',30);
% ylabel('Row degree distribution','fontsize',30);
% set(gca,'FontSize',24);
% %--------------------------------------------------------------------------
% 
% %--------------------Display Column Degree Distribution--------------------
% subplot(2,1,2)
% bar(deg_column,lambda);
% xlim([0 L]);
% xlabel('Degree','fontsize',30);
% ylabel('Column degree distribution','fontsize',30);
% set(gca,'FontSize',24);
% %--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
