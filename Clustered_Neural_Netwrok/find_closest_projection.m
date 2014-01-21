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
itr_max = 100;
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
varphi0 = 0.75;
psi = 0.005;

W_total_deg_dist_0 = sum(abs((W_global_tot)));
W_total_deg_dist = sum(abs((W_global_tot))>0);
%--------------------------------------------------------------------------

%==========================================================================

for mu = 1:dataset_size
    pattern = [1,dataset_in(mu,:)];
    
    for itr = 1:itr_max
        
        %------------Calculate the Feedback and Decision Parameters------------
        cost = (W_global_tot*pattern');
        c = (abs(cost) > psi).* (sign(cost));
        feedback = W_global_tot'*c;            
        normalized_feedback = feedback./(W_total_deg_dist_0');
        %----------------------------------------------------------------------
    
        %----------------------Update the Dataset------------------------------
        [val,ind] = max(abs(normalized_feedback));
        if (val > varphi)
            pattern(ind) =  min(max(pattern(ind)-sign(normalized_feedback(ind)),x_min),x_max);
        else
            break;
        end
            
        updated_neurons_orig = (abs(normalized_feedback)>varphi);    
   
        %     pattern = min(max(pattern-(updated_neurons_orig.*sign(normalized_feedback))',x_min),x_max);
        %----------------------------------------------------------------------
    
        %--------------Calculate the Number of Converged Patterns--------------                
        
        
        %----------------------------------------------------------------------
    end
    if (val < varphi)
        converged_count_orig = converged_count_orig+ 1;
        dataset_in(mu,:) = pattern(2:end);
        converged_count_per_itr = [converged_count_per_itr,converged_count_orig];    
    end
%     converged_count_orig = sum(sum(updated_neurons_orig)==0);    
    cost_per_itr = [cost_per_itr,mean(mean(abs(cost)))];
    itr_per_itr = [itr_per_itr,itr];
   
end

