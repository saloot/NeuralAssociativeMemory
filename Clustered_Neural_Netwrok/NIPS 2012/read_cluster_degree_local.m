%==========================================================================
%*******************FUNCTION: read_cluster_degree_local********************
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
% deg_row: The constraints degrees within clusters in a clustered neural associative memory
% deg_column: The pattern nodes degrees within clusters in a clustered neural associative memory
% lambda: The pattern nodes degree distribution within clusters in a clustered neural associative memory
% rho: The cluster degree distribution within clusters in a clustered neural associative memory
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function reads the result of the learning phase for a clustered
% neural associative memory and determines the degree distribution of the
% patterns and constraint neurons within clusters on average, i.e. by
% averaging over clusters and simulated instances. The degrees are also
% normalized to the interval [0,1] with respect to the maximum possible
% degree in the cluster for both pattern and constraint neurons. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function [deg_row,deg_column,lambda,rho] = read_cluster_degree_local(N,K,L,alpha0,beta0,theta0)

%==============================INITIALIZATION==============================
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           

rho = [0:.05:1];
deg_row = [0:.05:1];
lambda = [0:.05:1];
deg_column = [0:.05:1];
N_tot = N*L/deg_cluster_min;
N_const_tot = zeros(1,index_max);
%==========================================================================

%===================COMPUTE TOTAL NUMBER OF CONSTRAINTS====================
for ind = 1:index_max
    
    for l = 1:L
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(ind),'.txt'], 'r');
        
        if (fid >-1)
            fid_flag = 1;
            W = fscanf(fid, '%f',[N,inf]);
            W = W';
            [m,~] = size(W);
            N_const_tot(ind) = N_const_tot(ind) + m;
            fclose(fid);
        end
        
    end
end
%==========================================================================

%=======================CALCULATE DEGREE DISTRIBUTIONS=====================

%----------------------Read the Weight Matrices----------------------------
for ind = 1:index_max
    
    for l = 1:L
        
        
        load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);
        n = length(index_l);
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt'], 'r');
        if (fid >-1)            
            W = fscanf(fid, '%f',[n,inf]);       
            W = W';
            [m,n] = size(W);
            fclose(fid);
   
            %-------Calculate Constraint Neurons Degree Distribution-------
            s_row = sum(abs(sign(W')));
            [rho_temp,deg_row] = hist(s_row/n,deg_row);        
            rho = rho + rho_temp;       
            %--------------------------------------------------------------
        
            %--------Calculate Pattern Neurons Degree Distribution---------
            s_column = sum(sign(abs(W)));                         
            [lambda_temp,deg_column] = hist(s_column/m,deg_column);               
            lambda = lambda+lambda_temp;        
            %--------------------------------------------------------------
        end
    end    
end
%--------------------------------------------------------------------------
% sum(N_const_tot)/length(N_const_tot)

%--------------Normalize Row Degree Distribution of Cluster----------------
rho = rho/sum(rho);
%--------------------------------------------------------------------------


%------------------Normalize Column Degree Distribution--------------------
lambda = lambda/sum(lambda);
%--------------------------------------------------------------------------

%==========================================================================

%=============================SAVE RESULTS=================================
mkdir([pwd,'/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L)],'Read_Results');        % Create a specific folder for the current N and K
save([pwd,'/Recall_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Read_Results/sparsity_results_local_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_alpha0_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.mat'],...
    'lambda','rho','deg_column','deg_row');
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
