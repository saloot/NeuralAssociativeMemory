%==========================================================================
%*******************FUNCTION: read_cluster_degree_local********************
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


function [deg_row,deg_column,lambda,rho,d_max,d_min,d_ave,m,pe_max,pe_max2,pe_ave,pe_ave2] = read_cluster_degree_local_v2(cluster_size,no_constraints,no_clusters,alpha0,beta0,theta0,Q,learn_itr_max,source_depository)

%==============================INITIALIZATION==============================
m = zeros(1,no_clusters);
d_max = zeros(1,no_clusters);
d_min = zeros(1,no_clusters);
d_ave = zeros(1,no_clusters);
m_max = 0;
counter = 0;
%==========================================================================

%===================COMPUTE TOTAL NUMBER OF CONSTRAINTS====================    
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
    load(file_name);                
    W = W_total(:,2:end);
    [a,cluster_size] = size(W);                        
    
    if (a>m_max)                
        m_max = a;    
    end
    
    m(l) = a;    
end
%==========================================================================

%=======================CALCULATE DEGREE DISTRIBUTIONS=====================

%----------------------Read the Weight Matrices----------------------------
lambda = zeros(1,m_max+1);
deg_column = [0:m_max];
rho = zeros(1,cluster_size+1);
deg_row = [0:cluster_size];
pe_max = 0;
pe_max2 = 0;
pe_ave = 0;
pe_ave2 = 0;

p_err = zeros(1,no_clusters);

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
    load(file_name);                
    W = W_total(:,2:end);        
                               
    %-------Calculate Constraint Neurons Degree Distribution-------
    s_row = sum(abs(sign(W')));            
    [rho_temp,deg_row] = hist(s_row,deg_row);                    
    rho_temp = rho_temp(1:cluster_size+1);            
    rho(1:length(rho_temp)) = rho(1:length(rho_temp)) + rho_temp/m(l);                   
    %--------------------------------------------------------------
                    
    %--------Calculate Pattern Neurons Degree Distribution---------
    s_column = sum(sign(abs(W)));                                     
    [lambda_temp,deg_column] = hist(s_column,deg_column);                           
    lambda_temp = lambda_temp(1:m(l)+1);            
    lambda(1:length(lambda_temp)) = lambda(1:length(lambda_temp))+lambda_temp/cluster_size;                    
    %--------------------------------------------------------------
           
            
    dd = sum(abs(sign(W)));
    d_max(l) = max(dd);            
    mind = max(dd);            
    for i = 1:length(dd)    
        if (dd(i)>0)&&(dd(i)<mind)                    
            mind = dd(i);                
        end        
    end
    
    d_min(l) = mind;            
    d_ave(l) = sum(dd)/length(dd);

    pe = ( d_ave(l)/m(l) )^d_min(l);
               
    if (pe>pe_max2)           
        pe_max2 = pe;            
    end
        
    temp_val = ( d_max(l)/m(l) )^d_min(l);
            
    pe_ave = pe_ave + temp_val;            
    pe_ave2 = pe_ave2 + pe;
            
    if (temp_val > pe_max)
        pe_max = temp_val;            
    end
        
    p_err(l) = temp_val;
                      
    counter = counter + 1;
end
pe_ave  = pe_ave/counter;
pe_ave2 = pe_ave2/counter;
%--------------------------------------------------------------------------


%--------------Normalize Row Degree Distribution of Cluster----------------
rho = rho/(counter);
%--------------------------------------------------------------------------


%------------------Normalize Column Degree Distribution--------------------
lambda = lambda/(counter);
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
% xlim([0 m_max]);
% xlabel('Degree','fontsize',30);
% ylabel('Column degree distribution','fontsize',30);
% set(gca,'FontSize',24);
% %--------------------------------------------------------------------------


%--------------------------------------------------------------------------
