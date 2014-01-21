%=============================INITIALIZATION===============================

%---------------------------------Initialization---------------------------
req_conv_frac_per_node = .9;
fraction_of_convergence = .9;                                                  % Required fraction of patterns to converge

simulation_set = 4;
load(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set)]);
% load(['./Progressive_Simulation_Results/Big_Simulation_results_set_',num2str(simulation_set)]);

successful_params = [];
%--------------------------------------------------------------------------

%------------------------------Choose Dataset------------------------------
switch dataset_specifier
    case 0                    
        db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Mixed/CIFAR_10_Gray_Mixed_DB.mat'];  
        db_name_in = 'CIFAR_10_Gray_Mixed_DB';
    case 1
        db_file_in = ['/scratch/amir/Databases/CIFAR_10/Color_Mixed/CIFAR_10_Color_Mixed_DB.mat'];  
        db_name_in = 'CIFAR_10_Color_Mixed_DB';
    case 2
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/CIFAR_10_Gray_Mixed_DB_whitened.mat';
        db_name_in = 'CIFAR_10_Gray_Mixed_DB_whitened';
    case 3
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Mixed/Whitened/CIFAR_10_Color_Mixed_DB_whitened.mat';
        db_name_in = 'CIFAR_10_Color_Mixed_DB_whitened';
    case 4
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Mixed/OMP_5_8_50/CIFAR_10_Color_Mixed_OMP_5_8_50.mat';
        db_name_in = 'CIFAR_10_Mixed_OMP_5_8_50';
    case 5
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/Sparse Filtered/Layer_2/Pooled_3_3/Unnormalized/Pooled_Sparse_Filtered_CIFAR_10_Gray_Mixed_DB_Gray_Class_1_Train_Layer_2.mat';
        db_name_in = 'final_database_vectorized';
    case 6
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Class_6/Train_Set/OMP_5/Pooled/pooled_qudrant_class_6_omp_8_50_2500.mat';
        db_name_in = 'pooled_features';        
    case 7
        db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(1),'/Train_Set/OMP_7/CIFAR_10_train_gray_class_',num2str(1),'_OMP_7_dict_50_patch_8_stand.mat'];        
        db_name_in = ['rec_ims'];                    
    case 10
        db_file_in = ['/scratch/amir/Databases/STL_10/Gray_Mixed/STL_10_Gray_Mixed_DB.mat'];  
        db_name_in = 'STL_10_Gray_Mixed_DB';
    case 11
        db_file_in = ['/scratch/amir/Databases/STL_10/Color_Mixed/STL_10_Color_Mixed_DB.mat'];  
        db_name_in = 'STL_10_Color_Mixed_DB';
    case 12
        db_file_in = '/scratch/amir/Databases/STL_10/Gray_Mixed/Whitened/STL_10_Gray_Mixed_DB_whitened.mat';
        db_name_in = 'STL_10_Gray_Mixed_DB_whitened';
    case 13
        db_file_in = '/scratch/amir/Databases/STL_10/Color_Mixed/Whitened/STL_10_Color_Mixed_DB_whitened.mat';
        db_name_in = 'STL_10_Color_Mixed_DB_whitened';
    case 14    
        db_file_in = ['/scratch/amir/Databases/STL_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/STL_10_train_gray_scale_class_',num2str(target_class),'.mat'];            
        db_name_in = ['STL_10_Gray_DB_class_',num2str(target_class)];                        
    otherwise        
        error('Invalid dataset!')   
end
%--------------------------------------------------------------------------

load(db_file_in);    
eval(['dataset_learn = ',db_name_in,';']);          
eval(['clear ',db_name_in,';']);                       
[dataset_size,N] = size(dataset_learn);



%--------------Add a Column to the Dataset for the Threshold---------------
dataset_learn_orig = [ones(dataset_size,1),dataset_learn];
%--------------------------------------------------------------------------

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

    destination_folder = [db_file_in(1:i-1),'/Learn_Results/Zero_One/Clustered_Version'];

%==========================================================================


%===============================MAIN LOOP==================================
parameters = [];
total_cost = [];
average_cost = [];
average_cost2 = [];
good_parameters = [];
params_converged = {};
% params_converged = [];
counter_converged_jobs = 0;
converged_patterns_per_cluster = {};
for ii = 1:size(simulation_parameters,1)
    current_params = simulation_parameters(ii,:);

    eval(['cluster_size=',num2str(current_params(1)),';']);
    eval(['no_clusters=',num2str(current_params(2)),';']);
    eval(['const_learn=',num2str(current_params(3)),';']);    
    eval(['alpha0=',num2str(current_params(4)),';']);
    eval(['beta0=',num2str(current_params(5)),';']);
    eval(['theta0=',num2str(current_params(6)),';']);
    eval(['Q=',num2str(current_params(7)),';']);
    eval(['learn_itr_max=',num2str(current_params(8)),';']);    
       
    %--------------------------Quantize if Necessary---------------------------
    if (Q > 0)    
        dataset_learn = round(Q*tanh(sigma1*dataset_learn_orig-mean1));
        dataset_learn = dataset_learn;    
    else
        dataset_learn = dataset_learn_orig;
    end
    %--------------------------------------------------------------------------

    
    W_global_tot = [];  
    converged_patterns_per_cluster_mat = [];
    for cluster_index = 1:no_clusters         
        
        load(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(cluster_index)]);
        save_flag = save_place(end,2);            
        if (save_flag == 2)           
            destination_folder_save = [destination_folder,'/Partial_Convergence'];   
        else    
            destination_folder_save = destination_folder;
        end
                   
        
    
        file_name = [destination_folder_save,'/Simulation_Results_cluster_size_',num2str(cluster_size),...
            '_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
            num2str(theta0),'_clustere_',num2str(cluster_index),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];
        fid = fopen(file_name,'r');
                                    
        if (fid > -1)                  
            counter_converged_jobs = counter_converged_jobs + 1;                                  
            successful_params = [successful_params;current_params];
            clear 'best_y_thr'                
            load(file_name);     
            W_global_tot = [W_global_tot;W_global];    
            %----------------Original Converge Count-------
            max_y_threshold = 0.05;
            cost = (W_global*dataset_learn');            
            fraction_converged_const_per_node = sum(abs(cost)<max_y_threshold)/const_learn;             
            deg_dist_orig = sum(abs(sign(W_total))>0);
            converged_count_orig = sum(fraction_converged_const_per_node>req_conv_frac_per_node);            
            %-------------------------------------
            
            [best_weight_thr,good_weight_thr,converged_count_max,good_deg_dist,best_y_thr,good_y_thr] = find_best_thresholds_clustered(W_global,W_total,dataset_learn);                        
            params_converged{counter_converged_jobs,1} = converged_count_max;
            params_converged{counter_converged_jobs,2} = best_weight_thr;
            params_converged{counter_converged_jobs,3} = good_weight_thr;
            params_converged{counter_converged_jobs,4} = good_deg_dist;
            params_converged{counter_converged_jobs,5} = converged_count_orig;
            params_converged{counter_converged_jobs,6} = deg_dist_orig;
            params_converged{counter_converged_jobs,7} = best_y_thr;
            params_converged{counter_converged_jobs,8} = good_y_thr;
            
            %---Determine Which Patterns Have Converged for Each Cluster---
            W_global = soft_threshold_matrix(W_global,best_weight_thr);
            cost = (W_global*dataset_learn');            
            fraction_converged_const_per_node = sum(abs(cost)<(max(best_weight_thr-0.0025,0)))/const_learn; 
            converged_flag = (fraction_converged_const_per_node>req_conv_frac_per_node);
            converged_patterns_per_cluster_mat = [converged_patterns_per_cluster_mat,converged_flag'];
            %--------------------------------------------------------------

%             [perfect_y_thr,perfect_weight_thr,good_y_thr,good_weight_thr,best_converged_count_matrix,good_converged_count_matrix,best_mean_degrees,good_mean_degrees] = find_best_thresholds_clustered(W_global,W_total,dataset_learn,index_pattern_neurons);     
            
%             params_converged{counter_converged_jobs,1} = perfect_y_thr;
%             params_converged{counter_converged_jobs,2} = perfect_weight_thr;
%             params_converged{counter_converged_jobs,3} = good_y_thr;
%             params_converged{counter_converged_jobs,4} = good_weight_thr;                              
%             params_converged{counter_converged_jobs,5} = best_converged_count_matrix;
%             params_converged{counter_converged_jobs,6} = good_converged_count_matrix;
%             params_converged{counter_converged_jobs,7} = best_mean_degrees;
%             params_converged{counter_converged_jobs,8} = good_mean_degrees;
        end        
    end
    if (converged_count_max > 50)
%         ii
%         sum(converged_patterns_per_cluster_mat)
        
    end
    converged_patterns_per_cluster{ii,1} = converged_patterns_per_cluster_mat;
    converged_patterns_per_cluster{ii,2} = sum(converged_patterns_per_cluster_mat);
                
end
%==========================================================================
 converged_patterns_per_cluster

%==========================SORT THE RESULTS================================
best_params = [];
for i = 1:counter_converged_jobs
    converged_count = params_converged{i,1};
                
    best_weight_thr = params_converged{i,2};        
    good_weight_thr = params_converged{i,3};
    good_deg_dist = params_converged{i,4};                
    converged_count_orig = params_converged{i,5};
    deg_dist_orig = params_converged{i,6};
    current_params = successful_params(i,:);    
    
    if (converged_count > fraction_of_convergence*dataset_size)
        i
        best_params = [best_params;current_params,mean(good_deg_dist),mean(good_weight_thr)];
    end        
    
end


bb = cell2mat(params_converged(:,1));
[val,ind] = sort(bb,'descend');
quantized_good_params = [];
val_quantized = [];
for i = 1:counter_converged_jobs
    if (successful_params(ind(i),7)>0)        
        quantized_good_params = [quantized_good_params;successful_params(ind(i),:)];
        val_quantized = [val_quantized,val(i)];
    end
end
further_sim_params = [];
for i = 1:size(quantized_good_params,1)
    if ( (val(i) < .75 * val(1))) 
        break;
    else
        further_sim_params = [further_sim_params;successful_params(ind(i),:)];
    end
end
%==========================================================================
        