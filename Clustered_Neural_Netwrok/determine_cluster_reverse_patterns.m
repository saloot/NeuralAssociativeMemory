% function determine_cluster_reverse_patterns

%=============================INITIALIZATION===============================

%-------------------------Simulation Parameters----------------------------
try_itr = 0;                                                % The number of times we have tried to generate a cluster
try_itr_max = 2000;                                         % The maximum number of times we are going to try to generate "no_of_clusters" clusters
%--------------------------------------------------------------------------

%----------------------------Network Parameters----------------------------
N = 576;                                                    % The toal number of pattern neurons in the network

no_of_clusters = 50;                                        % The total number of clusters we would like to have
cluster_count = 0;                                          % The actual number of clusters we currently have

mean_cluster_size = 70;                                     % average number of PATTERN nodes in a cluster
max_cluster_dimension = 70;                                % The maximum number of CONSTRAINT nodes in a cluster


% alpha0 = 0.95;
% beta0= 0.75;   
% theta0 = 0.008;    
% target_class = 1;
%--------------------------------------------------------------------------
  

%--------------------Create the Sub-folder If Necessary--------------------
slash_flag = 0;
for i = length(weight_matrix_file):-1:1
    if (strcmp(weight_matrix_file(i),'/'))
        if (slash_flag == 0)
            break;
        else
            slash_flag = slash_flag+1;
        end
    end
end


destination_folder = [weight_matrix_file(1:i-1),'/Cluster_Indices'];
if (~exist(destination_folder,'dir'))
    mkdir(destination_folder);
end
%--------------------------------------------------------------------------


%-----------------------Load Input Weight Matrix---------------------------
fid = fopen(weight_matrix_file, 'r');               
if (fid > -1)                            
    W = fscanf(fid, '%f',[N,inf]);                                
    W = W';                                
    fclose(fid);                                                                            
    [no_of_constraints,~] = size(W);                      
else    
    error('Invalid input matrix');
end
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
    temp = randperm(N);
    cluster_index = temp(1:cluster_size);
    
    W_cluster = [];
    for i = 1:cluster_size
        W_cluster = [W_cluster,W(:,cluster_index(i))];
    end
    %----------------------------------------------------------------------
    
    
    
    %------------Check for the Success of Cluster Generation---------------
    cluster_dimension = correct here sum(sign(sum(abs(sign(W_cluster)))));
    if (cluster_dimension > max_cluster_dimension)
        success_flag = 0;
    end

    identical_columns_count = 0;
    for j = 1:cluster_size
        for k = j+1:cluster_size
            if (norm(sign(W_cluster(:,j)))>0)
                if (sign(W_cluster(:,j)')*W_cluster(:,k)/sum(abs(W_cluster(:,k))) > .95)
                    success_flag = 0;
                    break;
                end
            end
        end

        
        if (~success_flag)
            break
        end
    end        
    %----------------------------------------------------------------------
    
    %----------------Keep Track of Successful Attempts---------------------
    if (success_flag)
        cluster_count = cluster_count + 1; 
        %----------------Store The Indices for the Cluster-----------------
        save([destination_folder,'/cluster_indices_N_',num2str(N),'_L_',num2str(no_of_clusters)...
            ,'_cluster_index_',num2str(cluster_count),'.mat'],'W_cluster','cluster_index');
        %------------------------------------------------------------------
    end
    %----------------------------------------------------------------------
    
    
    
    
        
end
%==========================================================================
