m = 50;
p_converege = 0.8;

min_req_converged_const = .8*m;
temp = 0;
for i = ceil(min_req_converged_const):m
    temp = temp + nchoosek(m,i) * ((p_converege)^i) * ((1-p_converege)^(m-i));
end
temp


no_clusters = 12;
p_converege = .5;

min_req_converged_const = 4;%.15*no_clusters;
temp = 0;
for i = ceil(min_req_converged_const):no_clusters
    temp = temp + nchoosek(no_clusters,i) * ((p_converege)^i) * ((1-p_converege)^(no_clusters-i));
end
temp


converged_matrix = zeros(no_clusters,dataset_size);

dataset_recall_temp = dataset_recall(:,1+index_pattern_neurons);
dataset_recall_temp = [ones(dataset_size,1),dataset_recall_temp];
cc = W_total*dataset_recall_temp';
ccc = (abs(cc) > max_y_threshold).* (sign(cc));
feedback = W_total'*ccc;
feedback = feedback./((sum(abs(W_total))')*ones(1,dataset_size));

converged_matrix(l,:) = (sum((abs(feedback)>.25))==0);
