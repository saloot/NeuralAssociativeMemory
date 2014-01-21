% function determine_cluster_reverse_constraints_v1

%=============================INITIALIZATION===============================

%-------------------------Simulation Parameters----------------------------
try_itr = 0;                                                % The number of times we have tried to generate a cluster
try_itr_max = 20000;                                        % The maximum number of times we are going to try to generate "no_of_clusters" clusters
%--------------------------------------------------------------------------

%----------------------------Network Parameters----------------------------
N = 576;                                                    % The toal number of pattern neurons in the network

no_of_clusters = 250;                                        % The total number of clusters we would like to have
cluster_count = 0;                                          % The actual number of clusters we currently have

mean_cluster_size = 80;                                     % average number of CONSTRAINT nodes in a cluster
max_cluster_dimension = 250;                                % The maximum number of PATTERN nodes in a cluster

min_deg_const = 3;                                          % This is the minimum degree that the pattern neurons in each cluster should have

dataset_zero_thr = .075;
dimension_ok_count = 0;
identical_column_ount = 0;
alpha0 = 0.95;
beta0= 0.75;   
theta0 = 0.008;    
cluster_index_global = zeros(no_of_clusters,N+1);
mebership_index_global = zeros(1,N+1);
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

destination_folder = '/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_1/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Learn_Results/Zero_One/N_576/Cluster_Indices';
%--------------------------------------------------------------------------

%---------------------------Load the Training Dataset----------------------
% load(db_file_in);
% eval(['dataset_learn = ',db_name_in,';']);   
% [dataset_size,pattern_length] = size(dataset_learn);
% 
% dataset_learn = dataset_learn .*(abs(dataset_learn)>dataset_zero_thr);
% dataset_learn  = ones(dataset_size,pattern_length).*((dataset_learn>dataset_zero_thr)-(dataset_learn<-dataset_zero_thr));
% % a = sqrt(sum(dataset_learn'.*dataset_learn'));
% % dataset_learn = dataset_learn./(a'*ones(1,1024));
% % dataset_learn = dataset_learn.*(abs(dataset_learn)>.03);
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
% dataset_learn = [ones(dataset_size,1),dataset_learn];
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
%  
% % W = W_global;
% [no_of_constraints,~] = size(W);                      
% threshold_zero = .005;
% for i = 1:no_of_constraints
%     W(i,:) = soft_threshold(W(i,:),threshold_zero);
% end
% 
% max_y_threshold = 0.05;
% display(['-------Convergence!!!--------']);                
% display(['Converged count = ',num2str(sum(mean(abs(W*dataset_learn'))<max_y_threshold))]);        
% display(['Min degree of W = ',num2str(min(sum(abs(W)>0)))]);
% display('-----------------------------------');        
% display(' ');   
%--------------------------------------------------------------------------

%==========================================================================


%==========================TRY TO GENERATE CLUSTERS========================
while ( (cluster_count < no_of_clusters) && (try_itr <= try_itr_max))
    
    %------------------Within Loop Initializations-------------------------
    try_itr = try_itr + 1;
    success_flag = 1;                                       % Represents the success of producing new clusters
    %----------------------------------------------------------------------
    
    
    %---------------Generate Random Constraint Nodes Indices---------------
    cluster_size = mean_cluster_size + round(sqrt(mean_cluster_size)*randn);
    temp = randperm(no_of_constraints);
    cluster_index = temp(1:cluster_size);
    
    W_cluster = [];
    for i = 1:cluster_size
        W_cluster = [W_cluster;W(cluster_index(i),:)];        
    end
    %----------------------------------------------------------------------
    
    
    
    %----------------Check for Cluster Dimsion Constraint------------------
    cluster_dimension = sum(sum(abs(sign(W_cluster)))>=1);
    if (cluster_dimension > max_cluster_dimension)
        success_flag = 0;
    end

    if (~success_flag)    
        continue
    end
    
    dimension_ok_count = dimension_ok_count + 1;    
    %----------------------------------------------------------------------x
    
    %----------------Check for Minimum Degree Constraint-------------------        
%     W_cluster = W;
%     cluster_dimension = 577;
    no_identical_columns = 0;
    no_ok_columns = 0;
    varphi = .95;
    for j = 1:cluster_dimension
        if ( (abs(sign(W_cluster(:,j)')*sign(W_cluster(:,j)))/sum(abs((W_cluster(:,j)))) > varphi) && (sum(abs(sign(W_cluster(:,j))))>=2 ) )
            no_ok_columns = no_ok_columns+1;                            
        end
        for k = j+1:cluster_dimension            
            if (norm(sign(W_cluster(:,j)))>0)
                if ( (abs(sign(W_cluster(:,j)')*W_cluster(:,k))/sum(abs(W_cluster(:,k))) > varphi) && (sum(abs(sign(W_cluster(:,k))))>=2 ) )
%                     success_flag = 0;
                    
                    no_identical_columns = no_identical_columns+1;
                end
                
            end
        end

        
%         if (~success_flag)
%             break
%         end
    end
    
%     no_ok_columns
%     no_identical_columns
    if ( (no_identical_columns > 2) || (no_ok_columns < cluster_dimension/20))
        success_flag = 0;
    end
    if (~success_flag)    
        continue
    end
    
    identical_column_ount = identical_column_ount + 1;
    %----------------------------------------------------------------------
    
    %----------------Check for Minimum Degree Constraint-------------------
%     tt = sum(abs(sign(W_cluster)));
%     tt(~tt) = inf;
%     d_min = min(tt);
%     if (d_min < min_deg_const)
%         success_flag = 0;
%         continue;
%     end
    %----------------------------------------------------------------------
    
    %----------------Keep Track of Successful Attempts---------------------
    if (success_flag)                
        cluster_count = cluster_count + 1;
        display(' ')
        display(['Total clusters found so far = ',num2str(cluster_count)])
        display(['Number of times tried = ',num2str(try_itr)])
        display(['Number of pattern nodes covered so far = ',num2str(sum(sign(sum(cluster_index_global))))])
        display(['Min number of clusters per each pattern node = ',num2str(min(mebership_index_global))])
        
        s = sum(cluster_index_global);
        
        for j = 1:N+1
            if (  sum(abs(sign(W_cluster(:,j))))>=2 )
                cluster_index_global(cluster_count,j) = 1;
                mebership_index_global(j) = mebership_index_global(j) + 1;
            end
        end
        
        
        %----------------Store The Indices for the Cluster-----------------
        save([destination_folder,'/cluster_indices_N_',num2str(N),'_L_',num2str(no_of_clusters)...
            ,'_cluster_index_',num2str(cluster_count),'.mat'],'W_cluster','cluster_index');
        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
    
    
    
    
        
end
%==========================================================================
