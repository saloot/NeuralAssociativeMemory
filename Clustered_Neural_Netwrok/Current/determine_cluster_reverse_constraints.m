% function determine_cluster_reverse_constraints

%=============================INITIALIZATION===============================

addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_1/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_1_Train_Layer_2.mat';
db_name_in = 'final_database_vectorized';


%-------------------------Simulation Parameters----------------------------
try_itr = 0;                                                % The number of times we have tried to generate a cluster
try_itr_max = 20000;                                        % The maximum number of times we are going to try to generate "no_of_clusters" clusters
%--------------------------------------------------------------------------

%----------------------------Network Parameters----------------------------
N = 576;                                                    % The toal number of pattern neurons in the network

no_of_clusters = 250;                                        % The total number of clusters we would like to have
cluster_count = 0;                                          % The actual number of clusters we currently have

mean_cluster_size = 150;                                     % average number of CONSTRAINT nodes in a cluster
max_cluster_dimension = 240;                                % The maximum number of PATTERN nodes in a cluster

min_deg_const = 3;                                          % This is the minimum degree that the pattern neurons in each cluster should have

dataset_zero_thr = .075;
dimension_ok_count = 0;
identical_column_ount = 0;
alpha0 = 0.95;
beta0= 0.75;   
theta0 = 0.008;    
cluster_index_global = zeros(no_of_clusters,N+1);
% target_class = 1;
%--------------------------------------------------------------------------
  

%--------------------Create the Sub-folder If Necessary--------------------
% slash_flag = 0;
% for i = length(weight_matrix_file):-1:1
%     if (strcmp(weight_matrix_file(i),'/'))
%         if (slash_flag == 0)
%             break;
%         else
%             slash_flag = slash_flag+1;
%         end
%     end
% end
% 
% 
% destination_folder = [weight_matrix_file(1:i-1),'/Cluster_Indices'];
% if (~exist(destination_folder,'dir'))
%     mkdir(destination_folder);
% end
%--------------------------------------------------------------------------

%---------------------------Load the Training Dataset----------------------
load(db_file_in);
eval(['dataset_learn = ',db_name_in,';']);   
[dataset_size,pattern_length] = size(dataset_learn);
% 
dataset_learn = dataset_learn .*(abs(dataset_learn)>dataset_zero_thr);
dataset_learn  = ones(dataset_size,pattern_length).*((dataset_learn>dataset_zero_thr)-(dataset_learn<-dataset_zero_thr));
% % a = sqrt(sum(dataset_learn'.*dataset_learn'));
% % dataset_learn = dataset_learn./(a'*ones(1,1024));
% % dataset_learn = dataset_learn.*(abs(dataset_learn)>.03);
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
dataset_learn = [ones(dataset_size,1),dataset_learn];
%--------------------------------------------------------------------------

%-----------------------Load Input Weight Matrix---------------------------
% fid = fopen(weight_matrix_file, 'r');               
% if (fid > -1)                            
%     W = fscanf(fid, '%f',[N,inf]);                                
%     W = W';                                
%     fclose(fid);                                                                            
% else    
%     error('Invalid input matrix');
% end
 
W = W_global;
[no_of_constraints,~] = size(W);                      
threshold_zero = .01;
for i = 1:no_of_constraints
    W(i,:) = soft_threshold(W(i,:),threshold_zero);
end

max_y_threshold = 0.05;
display(['-------Convergence!!!--------']);                
display(['Converged count = ',num2str(sum(mean(abs(W*dataset_learn'))<max_y_threshold))]);        
display(['Min degree of W = ',num2str(min(sum(abs(W)>0)))]);
display('-----------------------------------');        
display(' ');   
%--------------------------------------------------------------------------

deg_orig = sum(abs(sign(W)));
[deg_sorted_orig,deg_index_orig] = sort(deg_orig);
deg_cluster_min = 4;
pattern_node_index = 0;
try_per_node_max = 40;
d_min = 6;                  % The desired minimum degree of nodes in a cluster
no_covered_nodes = [];
d_min_in_cluster = 4;
%==========================================================================


%==========================TRY TO GENERATE CLUSTERS========================
while ( (cluster_count < no_of_clusters) && (try_itr <= try_itr_max) )
    
    %------------------Within Loop Initializations-------------------------
    try_itr = try_itr + 1;
    success_flag = 1;                                       % Represents the success of producing new clusters
    pattern_node_index = pattern_node_index + 1;
    deg_cluster_node = 0;
    if (pattern_node_index>N+1)
        pattern_node_index = 1;
    end
    v = deg_index_orig(pattern_node_index);            
    if (v == 1)                      
        continue;        
    end
    %----------------------------------------------------------------------
    
    
    %---------------Generate Random Constraint Nodes Indices---------------
    try_per_node = 0;
    
    while  ( (try_per_node <try_per_node_max) && (deg_cluster_node <deg_cluster_min)) 
        cluster_dimension = 0;
        
        try_per_node = try_per_node + 1;        
        W_cluster = [];
        cluster_index = [];
        
        
        
        cluster_size = 0;
        pattern_node_index_within = 0;
        while (cluster_dimension <max_cluster_dimension )                        
            
            pattern_node_index_within = pattern_node_index_within+1;
            const_index = randi(no_of_constraints+1,no_of_constraints,1);
            deg_v = sum(abs(sign(W(:,v))));
        
            [consts,ind_consts] = sort(const_index.*abs(sign(W(:,v))));
            consts = ind_consts(no_of_constraints-d_min+1:no_of_constraints);
            for i = 1:length(consts)
                W_cluster = [W_cluster;W(consts(i),:)];        
                cluster_index = [cluster_index,consts(i)];  
                cluster_size = cluster_size + 1;
            end            
            
            %-------------Determine the Domain of the Cluster--------------
            if (cluster_size > 1)
                domain = sign(sum(abs(sign(W_cluster))));
            else
                domain = sign(W_cluster);
            end
            %--------------------------------------------------------------
                        
            %--------------Pick the Node with Minimum Degree---------------
            cluster_dimension = sum(sum(abs(sign(W_cluster)))>0);
            
            deg = sum(abs(sign(W))).*sign(sum(abs(sign(W_cluster))));            
            deg(~deg) = inf;
            [deg_sorted,deg_index] = sort(deg);
            v = deg_index(pattern_node_index_within);    
            if (v == 1)                
                v = deg_index(pattern_node_index_within+1);
            end
            %--------------------------------------------------------------
        end
        
        
        %----------------Check If We Have a Valid Cluster------------------
        no_identical_columns = 0;
        no_ok_columns = 0;
        varphi = .75;        
        success_flag = 1;
        for j = 1:cluster_dimension
            if ( (abs(sign(W_cluster(:,j)')*sign(W_cluster(:,j)))/sum(abs((W_cluster(:,j)))) > varphi) && (sum(abs(sign(W_cluster(:,j))))>=d_min_in_cluster-1 ) )
                no_ok_columns = no_ok_columns+1;                            
            end
            for k = j+1:cluster_dimension            
                if (norm(sign(W_cluster(:,j)))>0)
                    if ( (abs(sign(W_cluster(:,j)')*W_cluster(:,k))/sum(abs(W_cluster(:,k))) > varphi) && (sum(abs(sign(W_cluster(:,k))))>=d_min_in_cluster-1 ) )
%                         success_flag = 0;
                    
                        no_identical_columns = no_identical_columns+1;
                    end
                
                end
            end
        end
                
        if ( (no_identical_columns > cluster_dimension/10) || (no_ok_columns < cluster_dimension/15))
            success_flag = 0;
        end
           
%         no_identical_columns
%         no_ok_columns
        if (success_flag)                
            deg_cluster_node = deg_cluster_node + 1;
            cluster_count = cluster_count + 1;
            display(['Total clusters found so far = ',num2str(cluster_count)])
            display(['Number of times tried = ',num2str(try_itr)])
            display(['Number of pattern nodes covered so far = ',num2str(sum(sign(sum(cluster_index_global))))])
            no_covered_nodes = [no_covered_nodes,sum(sign(sum(cluster_index_global)))];
            for j = 1:N+1
                if (  sum(abs(sign(W_cluster(:,j))))>=d_min_in_cluster-1)
                    cluster_index_global(cluster_count,j) = 1;
                end
            end
            
            %----------------Store The Indices for the Cluster-----------------
%             save([destination_folder,'/cluster_indices_N_',num2str(N),'_L_',num2str(no_of_clusters)...
%                 ,'_cluster_index_',num2str(cluster_count),'.mat'],'W_cluster','cluster_index');
            %------------------------------------------------------------------
        end

        
    end
    %------------------------------------------------------------------
        
    
end
        

%==========================================================================
