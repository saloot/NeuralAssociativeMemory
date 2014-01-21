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


function [deg_row,deg_column,lambda,rho,d_max,d_min,d_ave,m,pe_max,pe_max2,pe_ave,pe_ave2,pe_2_max,pe_2_ave] = read_cluster_degree_local(N,K,L,alpha0,beta0,theta0)

%==============================INITIALIZATION==============================
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           

N_tot = N*L/deg_cluster_min;
N_const_tot = zeros(1,index_max);
m = zeros(index,L);
d_max = zeros(index,L);
d_min = zeros(index,L);
d_ave = zeros(index,L);
n_max = 0;
m_max = 0;
counter = 0;
%==========================================================================

%===================COMPUTE TOTAL NUMBER OF CONSTRAINTS====================
for ind = 1:index_max
    
    for l = 1:L
        load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);
        n = length(index_l);
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt'], 'r');
        if (fid >-1)
            fid_flag = 1;
            W = fscanf(fid, '%f',[n,inf]);
            W = W';
            [a,n] = size(W);
            if (n>n_max)
                n_max = n;
            end
            if (a>m_max)
                m_max = a;
            end
            m(ind,l) = a;
            N_const_tot(ind) = N_const_tot(ind) + m(ind,l);
            fclose(fid);
        end
        
    end
end
%==========================================================================

%=======================CALCULATE DEGREE DISTRIBUTIONS=====================

%----------------------Read the Weight Matrices----------------------------
lambda = zeros(1,m_max+1);
deg_column = [0:m_max];
rho = zeros(1,n_max+1);
deg_row = [0:n_max];
pe_max = 0;
pe_max2 = 0;
pe_ave = 0;
pe_ave2 = 0;
pe_1_temp_max = 0;
pe_2_temp_max = 0;
pe_2_max = 0;
pe_2_ave = 0;
pe_ave_over_clusters = 0;
p_err = zeros(index_max,L);
for ind = 1:index_max
    
    for l = 1:L
        
        
        load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'_cluster_',num2str(l),'.mat']);
        n = length(index_l);
        fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_',num2str(l),'_index_',num2str(index),'.txt'], 'r');
        if (fid >-1)            
            W = fscanf(fid, '%f',[n,inf]);       
            W = W';
            [a,n] = size(W);
            fclose(fid);
               
            %-------Calculate Constraint Neurons Degree Distribution-------
            s_row = sum(abs(sign(W')));
            [rho_temp,deg_row] = hist(s_row,deg_row);        
            rho_temp = rho_temp(1:n+1);
            rho(1:length(rho_temp)) = rho(1:length(rho_temp)) + rho_temp/a;       
            %--------------------------------------------------------------
        
            %--------Calculate Pattern Neurons Degree Distribution---------
            s_column = sum(sign(abs(W)));                         
            [lambda_temp,deg_column] = hist(s_column,deg_column);               
            lambda_temp = lambda_temp(1:a+1);
            lambda(1:length(lambda_temp)) = lambda(1:length(lambda_temp))+lambda_temp/n;        
            %--------------------------------------------------------------
           

            dd = sum(abs(sign(W)));
            d_max(ind,l) = max(dd);
            mind = max(dd);
            for i = 1:length(dd)
                if (dd(i)>0)&&(dd(i)<mind)
                    mind = dd(i);
                end
            end
            d_min(ind,l) = mind;
            d_ave(ind,l) = sum(dd)/length(dd);

%             pe = lambda_poly(1-rho_poly(1-1/n,rho_temp/a),lambda_temp/n);
            pe = lambda_poly(d_ave(ind,l)/m(ind,l),lambda_temp/n);
%             pe = ( d_ave(ind,l)/m(ind,l) )^d_min(ind,l);
           
            if (pe>pe_max2)
                pe_max2 = pe;
            end
            
            temp_val = ( d_max(ind,l)/m(ind,l) )^d_min(ind,l);
%             temp_val = lambda_poly(d_max(ind,l)/m(ind,l),lambda_temp/n);
            pe_ave = pe_ave + temp_val;
            pe_ave2 = pe_ave2 + pe;
            if (temp_val > pe_max)
                pe_max = temp_val;
            end            
            
             p_err(ind,l) = temp_val;
            %-----------------Correcting two error bits--------------------
%             pe_1 = (1-(1-d_ave(ind,l)/m(ind,l))^2)^d_min(ind,l);
            pe_1 = (min(m(ind,l),2*d_max(ind,l))/m(ind,l))^d_min(ind,l);
%             if (pe_1>pe_1_temp_max)
%                 pe_1_temp_max = pe_1;
%             end
            pe_2_tot = 0;
            for i = 1:length(lambda_temp)
                if (lambda_temp(i) > 0)
                    dp = deg_column(i);
                    pe_2 = 0;
                    for j = ceil(dp/2):dp
                        pe_2 = pe_2 + nchoosek(dp,j) * ((d_ave(ind,l)/m(ind,l))^j) * (1-(d_ave(ind,l)/m(ind,l)))^(dp-j);
                    end
                    pe_2_tot = pe_2_tot + pe_2*lambda_temp(i)/n;
                end
            end
            pe_temp = (1-2/n)*pe_1+(2/n)*pe_2_tot;
%             pe_2 = 0;
%             for i = ceil(d_min(ind,l)/2):d_min(ind,l)
%                 pe_2 = pe_2 + nchoosek(d_min(ind,l),i) * ((d_ave(ind,l)/m(ind,l))^i) * (1-(d_ave(ind,l)/m(ind,l)))^(d_min(ind,l)-i);
%             end
%             
% %             if (pe_2>pe_2_temp_max)
% %                 pe_2_temp_max = pe_2;
% %             end
% %             pe_2_max = (1-2/n)*pe_1_temp_max+(2/n)*pe_2_temp_max;
%             pe_temp = (1-2/n)*pe_1+(2/n)*pe_2;

            if (pe_temp>pe_2_max)
                pe_2_max = pe_temp;
            end
            pe_2_ave = pe_2_ave+ pe_temp;
            %--------------------------------------------------------------
            
            
            counter = counter + 1;
        end
    end    
end
pe_ave  = pe_ave/counter;
pe_ave2 = pe_ave2/counter;
pe_2_ave = pe_2_ave /counter;
%--------------------------------------------------------------------------
% sum(N_const_tot)/length(N_const_tot)

%--------------Normalize Row Degree Distribution of Cluster----------------
% rho = rho/sum(rho);
rho = rho/(counter);
%--------------------------------------------------------------------------


%------------------Normalize Column Degree Distribution--------------------
% lambda = lambda/sum(lambda);
lambda = lambda/(counter);
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
