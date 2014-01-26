%==========================================================================
%*********************FUNCTION: find_network_characteristics******************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% index: Index of the simulation ensemble
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% min_weights: The minimum absolute value of any connection in EACH cluster of the netwotk
% mean_weights: The average absolute value of all connections in EACH cluster of the netwotk
% d_max: The maximum degree of pattern neurons among ALL clusters
% d_min: The minimum degree of pattern neurons among ALL clusters
% d_ave_max: The maximum average degree of pattern neurons accross all clusters
% Pe_clusters: The upper bound on the probability of not correcting a single external error for each cluster
% no_consts: The number of constraints in each cluster
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function [min_weights,mean_weights,d_max,d_min,d_ave_max,Pe_clusters,no_consts] = find_network_characteristics(N,K,L,alpha0,beta0,theta0,index)

%=============================INITIALIZATION===============================
min_W = 10000;
min_weights = [];
mean_weights =[];
no_consts = [];
Pe_clusters = [];
d_max = 0;
d_min = 10000;
d_ave_max = 0;
Pe_max = 0;
%==========================================================================


%================================MAIN LOOP=================================
for cluster = 1:L                       
    
    %---------------Read the Already Learned Constraints-------------------
    load(['/scratch/amir/Fault_Tolerant_Decoding/Initialization_Files/N_',num2str(N),...
        '_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),...
        '_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(cluster),'.mat']);        
    
    n = length(index_l);                % Load the number of pattern nodes in the corresponding cluster
    
    fid = fopen(['/scratch/amir/Fault_Tolerant_Decoding/Learn_Results/N_',num2str(N),...
        '_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',...
        num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(cluster),'_index_',num2str(index),'.txt'], 'r');                
    
    if (fid > -1)                
        W = fscanf(fid, '%f',[n,inf]);                    
        W = W';                    
        fclose(fid);                                                                
        [m,~] = size(W);                
    else        
        m = 0;
    end
    
    no_consts = [no_consts,m];
    %----------------------------------------------------------------------
    
    
    %--------------------Find Minimum Weight Value-------------------------
    ww = W;        
    ww(~ww) = inf;                                
    if (min(min(abs(ww)))<min_W)
        min_W = min(min(abs(ww)));
        min_W_index = cluster;
    end
    min_weights = [min_weights,min(min(abs(ww)))];
    mean_weights = [mean_weights,mean(min(abs(ww)))];
    %----------------------------------------------------------------------
    
    %------------------Find Max, Min and Average Degree--------------------
    h = sum(abs(W)>0);
    if (max(h)>d_max)
        d_max = max(h);
        max_deg_index = cluster;
    end
    
    hh = h;        
    hh(~hh) = inf;                                
    if (min(hh)<d_min)
        d_min = min(hh);
        min_deg_index = cluster;
    end
    
    d_ave = sum(h)/n;           
    if (d_ave>d_ave_max)
        d_ave_max = d_ave;
        max_deg_ave_index = cluster;
    end
    %----------------------------------------------------------------------
    
    
    %------------------------Find Maximum Pe-------------------------------
    Pe = (d_ave/m)^min(hh);
    Pe_clusters = [Pe_clusters,Pe];
    if (Pe>Pe_max)
        Pe_max = Pe;
        Pe_max_index = cluster;
    end    
    %----------------------------------------------------------------------
        
end
%==========================================================================

