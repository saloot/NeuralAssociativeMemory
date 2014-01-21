function [max_no_correctable_patterns,max_y_threshold] = check_correction_one_error_clustered(W_total,dataset_size,cost,varphi)

[const_learn,n] = size(W_total);
W_total_deg_dist = sum(abs((W_total))>0);

req_conv_frac_per_node = .9;
fraction_of_convergence = .95;
max_no_correctable_patterns = 0;
% 
best_y = 0;
max_no_correctable_patterns = 0;
for max_y_threshold = .005:.002:.2                             
    W_total_deg_dist_0 = sum(abs((W_total)));
    feedback = W_total'*(sign(cost-max_y_threshold).*(abs(cost)>max_y_threshold));            
    normalized_feedback = feedback./(W_total_deg_dist_0'*ones(1,dataset_size));
    normalized_feedback = normalized_feedback.*(ones(dataset_size,1)*(W_total_deg_dist>2))';        

    updated_neurons_orig = (abs(normalized_feedback)>varphi);    
    converged_count_orig = sum(sum(updated_neurons_orig)==0);
    
    patterns_correctibility_matrix = zeros(dataset_size,n);
    for i = 1:n   
%         if (W_total_deg_dist(i) > 2)
        cost_temp = cost + W_total(:,i)*ones(1,dataset_size);
        
        feedback = W_total'*(sign(cost_temp-max_y_threshold).*(abs(cost_temp)>max_y_threshold));    
        normalized_feedback = feedback./(W_total_deg_dist_0'*ones(1,dataset_size));
        normalized_feedback = normalized_feedback.*(ones(dataset_size,1)*(W_total_deg_dist>2))';
        updated_neurons = (abs(normalized_feedback)>varphi);            
        converged_count = sum((abs(updated_neurons(i,:))==1).*(sum(abs(updated_neurons))==1));
        
        %-----------------------------WTA---------------------------------
        [val,ind] = max(abs(normalized_feedback));
        patterns_correctibility_matrix(:,i) = (ind==i);
        converged_count = sum(ind == i);
        %-----------------------------------------------------------
%         else
%             patterns_correctibility_matrix(:,i) = ones(dataset_size,1);
%         end
%         no_of_ok_patterns = [no_of_ok_patterns,converged_count];

    end
    correctable_patterns = (sum(patterns_correctibility_matrix')/n > 0.95).*(sum(updated_neurons_orig)==0);
    if (sum(correctable_patterns)> max_no_correctable_patterns)    
        max_no_correctable_patterns = sum(correctable_patterns);
        best_y = max_y_threshold;
    end
end

