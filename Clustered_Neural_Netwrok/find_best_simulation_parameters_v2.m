%=============================INITIALIZATION===============================

%---------------------------------Initialization---------------------------
req_conv_frac_per_node = .9;
fraction_of_convergence = .9;                                                  % Required fraction of patterns to converge
addpath(genpath('/home1/amir/cluster/Common_Library'));        
simulation_set = 15;
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
    case 15
        db_file_in = ['/scratch/amir/Databases/Caltech-101/Caltech101_Silhouettes/caltech101_silhouettes_28.mat'];
        db_name_in = ['X'];
    case 16
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/CIFAR_10_Gray_Mixed_DB_Q_15_Binary.mat';
        db_name_in = 'CIFAR_10_Gray_Mixed_DB_Q_15_Binary';
    case 17
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/CIFAR_10_Gray_Mixed_DB_Q_15_Binary_Projected.mat';
        db_name_in = 'CIFAR_10_Gray_Mixed_DB_Q_15_Binary_Projected';        
    case 18
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Class_6/Train_Set/OMP_5/Pooled/pooled_qudrant_class_6_omp_8_50_2500_vectorized_Q_15_Binary';
        db_name_in = 'pooled_features_Q_15_Binary';
    case 19
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Mixed_DB_Gray_Class_1_Train_Layer_2_Q_15_Binary.mat';
        db_name_in = 'final_database_vectorized_Q_15_Binary';
    case 20
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/CIFAR_10_Gray_Mixed_Whitened_DB_Q_15_Binary.mat';
        db_name_in = 'CIFAR_10_Gray_Mixed_Whitened_DB_Q_15_Binary';               
    case 20.5
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/CIFAR_10_Gray_Mixed_Whitened_DB_Q_7_Binary.mat';
        db_name_in = 'CIFAR_10_Gray_Mixed_Whitened_DB_Q_7_Binary'; 
    case 25
        db_file_in = '/scratch/amir/Databases/CIFAR_10/Gray_Mixed/Whitened/CIFAR_10_Gray_Mixed_Whitened_DB_Q_15_Binary_Projected.mat';
        db_name_in = 'CIFAR_10_Gray_Mixed_Whitened_DB_Q_15_Binary_Projected';        
    otherwise        
        error('Invalid dataset!')   
end
%--------------------------------------------------------------------------

load(db_file_in);    
eval(['dataset_learn = ',db_name_in,';']);          
eval(['clear ',db_name_in,';']);                       
[dataset_size,N] = size(dataset_learn);
dataset_learn = [ones(dataset_size,1),dataset_learn];

% if (Q>0)
%     %---------------Determine Quantization Parameters------------------
%     a = matrix2vector(dataset_learn);        
%     [d,c] = hist(abs(a),4000);
%     q1 = 0;        
%     for i = 1:length(d)
%         if (sum(d(i:end))/sum(d) <.975)                
%             q1 = c(i-1);                
%             ind1 = i-1;                
%             break;            
%         end        
%     end    
%     mean1 = q1;
%     q1 = q1-c(1);
%                         
%     for i = length(d):-1:1
%         if (sum(d(1:i))/sum(d) <.975)                
%             q2 = c(i+1);                
%             ind2 = i+1;                
%             break;            
%         end        
%     end
%         
%     sigma1 = 1/(q2-q1);                          % Map the value of q1 to 1
%     %------------------------------------------------------------------
% end

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

    destination_folder = [db_file_in(1:i-1),'/Learn_Results/Clustered_Version'];

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
%     if (Q > 0)    
%         dataset_learn = round(Q*tanh(sigma1*dataset_learn_orig-mean1));
%         dataset_learn = dataset_learn;    
%     else
%         dataset_learn = dataset_learn_orig;
%         
%     end
    %--------------------------------------------------------------------------

    if (Q == 0)
        continue;
    end
    W_global_tot = []; 
    converged_matrix = [];
    converged_matrix = zeros(no_clusters,dataset_size);
    converged_constraints_tot = [];
    converged_patterns_per_cluster_mat = [];
    super_graph_connectivity = [];
    for cluster_index = 1:no_clusters         
        
        fid = fopen(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(cluster_index),'.mat'],'r');
        if (fid > -1)
            fclose(fid);
            load(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(cluster_index)]);
        else
%             fclose(fid);
            continue;
        end
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
            fclose(fid);
            
            successful_params = [successful_params;current_params];
            clear 'best_y_thr'                
            load(file_name); 
            temp = zeros(1,N);
            
            cost = W_global*dataset_learn';
            [mi,~] = size(W_global);  
            if (mi > 1)
                
                dataset_recall_temp = dataset_recall(:,1+index_pattern_neurons);
                dataset_recall_temp = [ones(dataset_size,1),dataset_recall_temp];
                cc = W_total*dataset_recall_temp';
                ccc = (abs(cc) > max_y_threshold).* (sign(cc));
                feedback = W_total'*ccc;
                feedback = feedback./((sum(abs(W_total))')*ones(1,dataset_size));
    
                converged_matrix(cluster_index,:) = (sum((abs(feedback)>.25))==0);

                
                temp(index_pattern_neurons) = 1;
                counter_converged_jobs = counter_converged_jobs + 1;                                  
%                 converged_matrix = [converged_matrix;sum(abs(cost)<.045)>=mi * 0.75];
                super_graph_connectivity = [super_graph_connectivity;temp];
            end
            W_global_tot = [W_global_tot;W_global];   
            [m,~] = size(W_global);
            converged_constraints_tot = [converged_constraints_tot,m];
            
%            W_global = soft_threshold_matrix(W_global,.002);
%         W_global = W_global./(sqrt(sum(W_global'.*W_global'))'*ones(1,size(W_global,2)));        
%         W_global_deg_dist = sum(abs(sign(W_global)));
%         
%         %---------------------Remove Zero Columns--------------------------
%         [val,ind] = sort(W_global_deg_dist);
%         for i = 1:length(val)
%             if (val(i) > 0)
%                 break;
%             end
%         end
%         W_global = W_global(:,ind(i:end));
%         W_global_deg_dist = sum(abs(sign(W_global)));          
%         %----------------------------------------------------------------------
% 
%         
%         dataset_temp = dataset_learn(:,ind(i:end));
%         cost = (W_global*dataset_temp');            
%         x_min = min(min(dataset_learn));
%         x_max = max(max(dataset_learn));
% 
%         %------------Find closest Patterns to the Matrix W_global_tot----------
%         [no_learned_patterns,max_y_threshold,dataset_temp_out] = find_closest_projection(W_global,dataset_size,dataset_temp,varphi,x_max,x_min);
%         %----------------------------------------------------------------------

            
        else
%             fclose(fid);
        end
    end
    
    
    %-------DETERMINE THE MIN AND MAX ACCEPTABLE WEIGHT THRESHOLDS---------
    max_weight_threshold = 0;
    min_weight_threshold = 0;
    threshold_flag = 0;
    W_global_tot_orig = W_global_tot;
    max_no_correctable_patterns = 0;
        
    dd = converged_matrix'*super_graph_connectivity./(ones(dataset_size,1)*sum(super_graph_connectivity));
    sum(sum(dd'>.25)>.75*4096)
    111;
%     for weight_thr = 0:.001:.25        
%         W_global = soft_threshold_matrix(W_global_tot_orig,weight_thr);        
%         W_global = W_global./(sqrt(sum(W_global'.*W_global'))'*ones(1,size(W_global_tot,2)));        
%         W_global_deg_dist = sum(abs(sign(W_global)));
%         
%         %---------------------Remove Zero Columns--------------------------
%         [val,ind] = sort(W_global_deg_dist);
%         for i = 1:length(val)
%             if (val(i) > 6)
%                 break;
%             else
%                 W_global(:,ind(i)) = W_global(:,ind(i)) + .05*random_vector(size(W_global,1),6-val(i))';
%             end
%         end
% %         W_global = W_global(:,ind(i:end));        
%         W_global = W_global./(sqrt(sum(W_global'.*W_global'))'*ones(1,size(W_global_tot,2)));        
%         W_global_deg_dist = sum(abs(sign(W_global)));
%         %----------------------------------------------------------------------
% 
%         
% %         dataset_temp = dataset_learn(:,ind(i:end));
%         dataset_temp = dataset_learn;
%         cost = (W_global*dataset_temp');            
%         x_min = min(min(dataset_learn));
%         x_max = max(max(dataset_learn));
% 
%         %------------Find closest Patterns to the Matrix W_global_tot----------
%         [no_learned_patterns,max_y_threshold,dataset_temp_out] = find_closest_projection(W_global,dataset_size,dataset_temp,varphi,x_max,x_min);
%         %----------------------------------------------------------------------
% 
%         dataset_out = dataset_learn;
%         dataset_out(:,ind(i:end)) = dataset_temp_out;
%         
%         if (no_learned_patterns > max_no_correctable_patterns)
%             max_no_correctable_patterns = max_no_correctable_patterns;
%             best_y = max_y_threshold;
%             best_weight = weight_thr;
%         end
%     
%     end
    111;
    
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
        