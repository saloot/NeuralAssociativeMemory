cluster_size = 55;
cluster_size = 420;
cluster_size = 80;
cluster_size = 90;
cluster_size = 50;

no_clusters = 140;
no_clusters = 40;
no_clusters = 80;
no_clusters = 150;
no_clusters = 100;

const_learn = min(cluster_size/2,30);
const_learn = 2*cluster_size;
const_learn = 160;
const_learn = 60;
const_learn = 60;
const_learn = 40;

theta0 = 0.01;

learn_itr_max = 250;
learn_itr_max = 750;

Q = 1;
Q = 16;
Q = 8;

db_file_in = ['/scratch/amir/Databases/Caltech-101/Caltech101_Silhouettes/caltech101_silhouettes_28.mat'];
db_name_in = ['X'];

db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/Sparse Filtered/Layer_2/Pooled_3_3/Unnormalized/Pooled_Sparse_Filtered_CIFAR_10_Gray_Mixed_DB_Gray_Class_1_Train_Layer_2.mat';
db_name_in = 'final_database_vectorized';


db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Mixed/CIFAR_10_Gray_Mixed_DB.mat'];  
db_name_in = 'CIFAR_10_Gray_Mixed_DB';

load(db_file_in);
eval(['dataset_learn = ',db_name_in,';']);   
[dataset_size,~] = size(dataset_learn);

%------------------Create Subdirectories if Necessary----------------------
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

destination_folder = [db_file_in(1:i-1),'/Learn_Results/Clustered_Version'];
if (~exist(destination_folder,'dir'))        
    mkdir(destination_folder);    
end
%--------------------------------------------------------------------------

counter = 0;
converged_pattern_index = [];
W_global_tot = [];
deg_dist_tot = zeros(1,cluster_size+1);

for l = 1:no_clusters
    file_name = [destination_folder,'/Simulation_Results_cluster_size_',...
        num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_0.95_beta_1_theta_',num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];
    fid = fopen(file_name,'r');

    if (fid == -1)                    
        file_name = [destination_folder,'/Partial_Convergence/Simulation_Results_cluster_size_'...
            ,num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_0.95_beta_1_theta_',num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];
        fid = fopen(file_name,'r');
    end
    
    if (fid > -1)
        counter = counter + 1;
        load(file_name)
        fclose(fid);
        [mi,~] = size(W_global);
        W_global = soft_threshold_matrix(W_global,2*max_y_threshold+.001);
        fraction_converged = sum(abs(W_global*[ones(dataset_size,1),dataset_learn]')<max_y_threshold)/mi;
        if (mi > 1)
            converged_pattern_index = [converged_pattern_index;fraction_converged];
            deg_dist = sum(abs(sign(W_total)));
            deg_dist_tot = deg_dist_tot + deg_dist;
        end
        W_global_tot = [W_global_tot;W_global];
        if (sum(deg_dist == 0)>0)
            111;
        end
    end
end

deg_dist_tot = deg_dist_tot/counter;
sum(sum(converged_pattern_index>.65)>2)

index_converged_patterns = [];
converged_pats = (sum(converged_pattern_index>.85)>5);
[dataset_size,~] = size(dataset_learn);
for i = 1:dataset_size
    if (converged_pats(i) )
        index_converged_patterns = [index_converged_patterns,i];
    end
end