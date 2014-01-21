% function [max_no_correctable_patterns,max_y_threshold,dataset_out] = find_closest_projection(W_global_tot,dataset_size,dataset_in_orig,varphi_in,x_max,x_min)

%=============================INITIALIZATION===============================

%---------------------------Load the Input Dataset-------------------------
load(db_file_in);
eval(['dataset_in = ',db_name_in,';']);
% % eval(['clear', db_name_in]);
% % [dataset_size,N] = size(dataset_in_orig);
% % dataset_out = dataset_in_orig;
% dataset_in = [ones(dataset_size,1),dataset_in];
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
itr_max = 20;
[const_learn,n] = size(W_global_tot);

req_conv_frac_per_node = .9;
fraction_of_convergence = .95;
max_no_correctable_patterns = 0;
% 
best_y = 0;
converged_count_orig = 0;
max_no_correctable_patterns = 0;
converged_count_per_itr = [];    
cost_per_itr = [];
itr_per_itr = [];
varphi_inv = 0.25;
psi = 0.005;
try_itr_max = 10;


%--------------------------------------------------------------------------

%==========================================================================

for mu = 1:dataset_size
    pattern = dataset_in(mu,:);
    
            
    success_flag = 0;
    try_itr = 0;
    while ( (success_flag == 0) && (try_itr < try_itr_max) )
        try_itr = try_itr+1;        
        consts_so_far = 1;
        success_flag = 1;
        for l = 1:no_clusters
                               
            consts = converged_constraints_tot(l);                    
            if (consts > 1)                    
                index_pattern_neurons = index_pattern_tot(l,:);                
                W_total = W_total_tot(consts_so_far:consts_so_far+consts-1,:);                
                consts_so_far = consts_so_far + consts;            
            else            
                consts_so_far = consts_so_far + consts;
                continue;            
            end
                        
            x_temp = [thr1,pattern(index_pattern_neurons)];                                   
            W_total_deg_dist_0 = sum(abs((W_total)));
            updated_flag = 0;
            W_total_deg_dist = sum(abs((W_total))>0);
            for itr = 1:itr_max
        
                %------------Calculate the Feedback and Decision Parameters------------
                cost = (W_total*x_temp');
                c = (abs(cost) > psi).* (sign(cost));
                feedback = sign(W_total)'*c;            
                normalized_feedback = feedback./(W_total_deg_dist');
                %----------------------------------------------------------------------
    
                %----------------------Update the Dataset------------------------------
                [val,ind] = max(abs(normalized_feedback));
                if (val > varphi_inv)
                    x_temp(ind) =  min(max(x_temp(ind)-sign(normalized_feedback(ind)),x_min),x_max);
                    updated_flag = 1;
                else
                    break;
                end
            
                updated_neurons_orig = (abs(normalized_feedback)>varphi_inv);    
               
            end
            if (val < varphi_inv)
                
                dataset_in(mu,index_pattern_neurons) = x_temp(2:end);
                if (updated_flag)
                    success_flag = 0;
                end
            end
        
            
        end
    end       
    if (success_flag)
        converged_count_orig = converged_count_orig+ 1;        
    end
    converged_count_per_itr = [converged_count_per_itr,converged_count_orig]    
end
