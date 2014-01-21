% function [max_no_correctable_patterns,max_y_threshold,dataset_out] = find_closest_projection(W_total,dataset_size,dataset_in_orig,varphi_in,x_max,x_min)

%=============================INITIALIZATION===============================

%---------------------------Load the Input Dataset-------------------------
load(db_file_in);
eval(['dataset_in_orig = ',db_name_in,';']);   
[dataset_size,N] = size(dataset_in_orig);
dataset_out = dataset_in_orig;
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
itr_max = 10;
[const_learn,n] = size(W_total);
W_total_deg_dist = sum(abs((W_total))>0);

req_conv_frac_per_node = .9;
fraction_of_convergence = .95;
max_no_correctable_patterns = 0;
% 
best_y = 0;
max_no_correctable_patterns = 0;
%--------------------------------------------------------------------------

%==========================================================================


%=======================FIND P==================================
ww=W_total;
ww(~ww) = inf;
W_Min = min(min(abs(ww)));
W_Min = sort(min(abs(ww)));
for i = 319:length(W_Min)
    if ( (W_Min(i)==Inf) || (W_Min(i)>.1) )
        break;
    end
end
W_Min = mean(W_Min(1:i-1));
%==========================================================================


for max_y_threshold = max(W_Min-.01,0.001):.005:max(W_Min-.01,0.001)
%     dataset_in = dataset_in_orig;
    max_i = 1200;
%      converged_count_per_itr = [];
%    cost_per_itr = [];
%  converged_count_not_learned_per_itr = [];
    varphi0 = .95;
    phi_itr = 0;
    for i = 2:max_i
%         varphi = min(1-3/(4+sqrt(i)),.885);
        if (mod(i,75) == 2)
            phi_itr = phi_itr + 1;            
            varphi = max(5*varphi0/(4+phi_itr),.5);
        end
        
%             varphi = 0.95/(sqrt(log(i)));
        cost = (W_total*dataset_in');    
        W_total_deg_dist_0 = sum(abs((W_total)));
        feedback = W_total'*(sign(cost-max_y_threshold).*(abs(cost)>max_y_threshold));            
        normalized_feedback = feedback./(W_total_deg_dist_0'*ones(1,dataset_size));
        if ( (i< .25*max_i) || ( (i > .5*max_i) && (i< .5*max_i)) )
            normalized_feedback = normalized_feedback.*(ones(dataset_size,1)*(W_total_deg_dist>15/(log(1+phi_itr))))';        
        end

%         [val,ind] = max(abs(normalized_feedback));
%         dataset_in(:,ind) =  min(max(dataset_in(:,ind)-sign(normalized_feedback(ind,:)'),x_min),x_max);
        
        updated_neurons_not_learned = (abs(normalized_feedback)>varphi).*(abs(normalized_feedback)<.8);    
        updated_neurons_orig = (abs(normalized_feedback)>varphi);    
        
        dataset_in = min(max(dataset_in-(updated_neurons_orig.*sign(normalized_feedback).*(ones(N+1,1)*(sum(updated_neurons_orig)<5*N/sqrt(i))))',x_min),x_max);
        converged_count_orig = sum(sum(updated_neurons_orig)==0);
        converged_count_per_itr = [converged_count_per_itr,converged_count_orig];
        converged_count_not_learned_per_itr = [converged_count_not_learned_per_itr, sum(sum(updated_neurons_not_learned)==0)];
        cost_per_itr = [cost_per_itr,mean(mean(abs(cost)))];
        111;
        i
    end
%     dataset_out = dataset_in;
%     patterns_correctibility_matrix = zeros(dataset_size,n);
%     correctable_patterns = (sum(patterns_correctibility_matrix')/n > 0.95).*(sum(updated_neurons_orig)==0);
%     if (sum(correctable_patterns)> max_no_correctable_patterns)    
    if (converged_count_orig > max_no_correctable_patterns)
        max_no_correctable_patterns = sum(converged_count_orig);
        dataset_out = dataset_in;
    end
end

