function [best_y,best_weight] = find_best_thresholds_clustered(W_global_orig,W_total_orig,dataset_in)

     
%%
%=============================INITIALIZATION===============================

%--------------------------Simulation Parameters---------------------------
fraction_of_convergence = .9;
req_conv_frac_per_node = 0.9;
converged_count_max = 0;
best_weight_thr = 0;
best_y_thr =0;
best_weight = 0;
best_y = 0;
[dataset_size,~] = size(dataset_in);
[const_learn,n] = size(W_total_orig);
good_y_thr = [];
good_weight_thr = [];
good_deg_dist = [];
ww = W_global_orig;        
ww(~ww) = inf;       
n_l2 = size(W_global_orig,2);
min_wight = min(min(abs(ww)));
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
%--------------------------------------------------------------------------



%----------------Conduct Prechecks on the Weight Matrix--------------------
W_global_orig_deg_dist = sum(abs(sign(W_total_orig))>0);
frac_zero_nodes_orig = .00001 + sum(W_global_orig_deg_dist==0)/n;        
if (frac_zero_nodes_orig > .5)
    bad_matrix_flag = 1;    
    display('The original W_global matrix has a lot of zero degree pattern neurons');   
    return
% else    
%     frac_zero_nodes_orig
%     display('The original W_global matrix is good!');   
end
%--------------------------------------------------------------------------

%==========================================================================


%%
%==========DETERMINE THE MIN AND MAX ACCEPTABLE WEIGHT THRESHOLDS==========
max_weight_threshold = 0;
min_weight_threshold = 0;
threshold_flag = 0;
for weight_thr = 0.001:.001:.25
        
    W_total = soft_threshold_matrix(W_total_orig,weight_thr);
    W_total = W_total./(sqrt(sum(W_total'.*W_total'))'*ones(1,size(W_total,2)));
    W_global_deg_dist = sum(abs(sign(W_total)));
    
    frac_low_degree_nodes = sum( (W_global_deg_dist <=1).*(W_global_deg_dist >0))/n;            
    frac_zero_degree_nodes = sum(W_global_deg_dist == 0)/n;            
    frac_high_degree_nodes = sum(W_global_deg_dist>=const_learn-1)/n;
    
    if (threshold_flag == 0)
        if ( (frac_high_degree_nodes < .15) && (frac_low_degree_nodes < .15) && (frac_zero_degree_nodes < .45) )
            min_weight_threshold = weight_thr;
            threshold_flag = 1;
        end
    else
        if ( (frac_high_degree_nodes > .15) || (frac_zero_degree_nodes >= .45) || (frac_low_degree_nodes >= .08) )
            max_weight_threshold = weight_thr;
            break;
        end
    end
end
%==========================================================================


max_no_correctable_patterns = 0;
%%
%========DETERMINE THE BEST UPDATE THRESHOLD FOR CONSTRAINT NODES==========
for weight_thr = 0:.005:max_weight_threshold    
% for max_y_threshold = .001:.005:.2    
    W_global = soft_threshold_matrix(W_global_orig,weight_thr);
    W_global = W_global./(sqrt(sum(W_global'.*W_global'))'*ones(1,n_l2));
    
    W_total = soft_threshold_matrix(W_total_orig,weight_thr);
    W_total = W_total./(sqrt(sum(W_total'.*W_total'))'*ones(1,size(W_total,2)));
    W_total_deg_dist = sum(abs(sign(W_total)));
    W_global_deg_dist = sum(abs(sign(W_global)));
        
    
    %-----------------------Remove Zero Columns----------------------------
    [val,ind] = sort(W_total_deg_dist);
    for i = 1:length(val)
        if (val(i) > 0)
            break;
        end
    end
    W_total = W_total(:,ind(i:end));
    W_total_deg_dist = sum(abs(sign(W_total)));
   
    dataset_in_temp = dataset_in(:,ind(i:end));
    cost = (W_total*dataset_in_temp');            
    varphi = 0.75;
    [no_correctable_patterns,max_y_threshold] = check_correction_one_error_clustered(W_total,dataset_size,cost,varphi);
    
    
    if (no_correctable_patterns > max_no_correctable_patterns)
        max_no_correctable_patterns = no_correctable_patterns;
        best_y = max_y_threshold;
        best_weight = weight_thr;
    end
    
    
    
end
%==========================================================================





% function [perfect_y_thr,perfect_weight_thr,good_y_thr,good_weight_thr,best_converged_count_matrix,good_converged_count_matrix,best_mean_degrees,good_mean_degrees] = find_best_thresholds_clustered(W_global_orig,W_total_orig,dataset_in,index_pattern_neurons)
% 
% 
%      
% %%
% %=============================INITIALIZATION===============================
% 
% %--------------------------Simulation Parameters---------------------------
% fraction_of_convergence = .9;
% req_conv_frac_per_node = 0.9;
% converged_count_max = 0;
% good_weight_thr = 0;
% 
% [dataset_size,~] = size(dataset_in);
% [const_learn,n] = size(W_total_orig);
% 
% good_weight_thr = [];
% good_y_thr = [];
% good_deg_dist = [];
% perfect_y_thr = [];
% perfect_weight_thr = []; 
% good_mean_degrees = [];
% best_mean_degrees = [];
% good_converged_count_matrix = [];
% best_converged_count_matrix = [];
% 
% ww = W_global_orig;        
% ww(~ww) = inf;       
% 
% min_wight = min(min(abs(ww)));
% %--------------------------------------------------------------------------
% 
% %-------------------------Other Initializations----------------------------
% addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
% %--------------------------------------------------------------------------
% 
% 
% 
% %----------------Conduct Prechecks on the Weight Matrix--------------------
% W_global_orig_deg_dist = sum(abs(sign(W_total_orig)));
% frac_zero_nodes_orig = .00001 + sum(W_global_orig_deg_dist==0)/n;        
% if (frac_zero_nodes_orig > .5)
%     bad_matrix_flag = 1;    
%     display('The original W_global matrix has a lot of zero degree pattern neurons');   
%     return
% % else    
% %     frac_zero_nodes_orig
% %     display('The original W_global matrix is good!');   
% end
% %--------------------------------------------------------------------------
% 
% %==========================================================================
% 
% 
% %%
% %==========DETERMINE THE MIN AND MAX ACCEPTABLE WEIGHT THRESHOLDS==========
% max_weight_threshold = 0;
% min_weight_threshold = 0;
% threshold_flag = 0;
% for weight_thr = min_wight:.005:.5
%         
%     W_total = soft_threshold_matrix(W_total_orig,weight_thr);
%     W_global_deg_dist = sum(abs(sign(W_total)));
%     
%     frac_low_degree_nodes = sum(W_global_deg_dist <=1)/n;            
%     frac_high_degree_nodes = sum(W_global_deg_dist>=const_learn)/n;
%     
%     if (threshold_flag == 0)
%         if ( (frac_high_degree_nodes < .1) && (frac_low_degree_nodes < .35) )
%             min_weight_threshold = weight_thr;
%             threshold_flag = 1;
%         end
%     else
%         if ( (frac_high_degree_nodes > .1) || (frac_low_degree_nodes >= .35) )
%             max_weight_threshold = weight_thr;
%             break;
%         end
%     end
% end
% %==========================================================================
% 
% 
% success_flag = 0;
% perfect_y_thr = [];
% good_weight_thr = [];
% perfect_weight_thr = [];
% deg_zero_flags = [];
%         
% 
% %%
% %========DETERMINE THE BEST UPDATE THRESHOLD FOR CONSTRAINT NODES==========
% for weight_thr = min_weight_threshold:.005:max_weight_threshold    
%     W_global = soft_threshold_matrix(W_global_orig,weight_thr);
%     W_total = soft_threshold_matrix(W_total_orig,weight_thr);
%     W_global_deg_dist = sum(abs(sign(W_total)));
%     
%     for max_y_threshold = .001:.005:.2                  
%         cost = (W_global*dataset_in');            
%         fraction_converged_const_per_node = sum(abs(cost)<max_y_threshold)/const_learn; 
%         converged_count = sum(fraction_converged_const_per_node>req_conv_frac_per_node);        
%         
%         if (converged_count > fraction_of_convergence*dataset_size)           
%             success_flag = 1;
%             break;                    
%         end
%     end
%     
%     if (success_flag) 
%         if (max_y_threshold < weight_thr-.001)
%             perfect_y_thr = [perfect_y_thr,max_y_threshold];
%             perfect_weight_thr = [perfect_weight_thr,weight_thr];
%             best_mean_degrees = [best_mean_degrees,mean(W_tilde_deg_dist)];
%             best_converged_count_matrix = [best_converged_count_matrix,converged_count];
%         else
%             [one_error_flag,deg_zero_nodes] = check_correction_one_error_clustered(W_global,max_y_threshold,dataset_in,index_pattern_neurons);
%             if (one_error_flag)
%                 good_y_thr = [good_y_thr,max_y_threshold];
%                 good_weight_thr = [good_weight_thr,weight_thr];                
%                 good_mean_degrees = [good_mean_degrees,mean(W_tilde_deg_dist)];
%                 good_converged_count_matrix = [good_converged_count_matrix,converged_count];
%             end
%         end
%     end
%     
% %     if (converged_count > converged_count_max)
% %         converged_count_max = converged_count;
% %         good_weight_thr = weight_thr;
% %     end
% %     
% %     if (converged_count > fraction_of_convergence*dataset_size)                       
% %         good_weight_thr = [good_weight_thr,weight_thr];
% %         good_deg_dist = [good_deg_dist;W_global_deg_dist];
% %     end
% end
% %==========================================================================
% 
