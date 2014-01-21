%==========================================================================
%***************FUNCTION: clustered_neural_recall**************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% no_clusters: The number of clusters if we have multiple levels (no_clusters = 1 for single level)
% no_simulated_instances: The number of patterns considered in the recall phase.
% err_bits: Number of initial noisy nodes.
% max_noise_amp: The maximum amplitude of noise during the recall process
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
% recall_algorithm_option: The recall algorithm identifier (0 for winner-take-all, 1 for the original bit flipping and 2 for the simplified bit flipping
% try_max: The maximum number of iterations in outside-cluster recall algorithm
% index: The index of the simulation setup among various random scenarios
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a clustered neural associative
% memory and performs the recall phase. More specifically, 
% the function first reads the weight matrix found in the learning phase. 
% Then the function performs the recall phase introduced in our NIPS 2012
% paper, i.e. perform the recall algorithm within each cluster one after
% another. If the recall process was successful, the state of neurons is
% maintained, and reverted back to their original version otherwise. This
% process is repeated try_max times, after which a recall error is declared
% if the output of the algorithm is not equal to the original pattern.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function clustered_neural_recall(cluster_size,const_learn,no_clusters,db_file_in,db_name_in,...
    no_simulated_instances,err_bits,max_noise_amp,alpha0,beta0,theta0,...
    varphi,psi,try_max,Q,learn_itr_max,simulation_set)



%=============================INITIALIZATION===============================

%---------------------------Load the Training Dataset----------------------
load(db_file_in);
eval(['dataset_recall = ',db_name_in,';']);   
[dataset_size,N] = size(dataset_recall);
dataset_recall_projected = dataset_recall;
%--------------------------------------------------------------------------

%------------------------Compute the Projected Dataset---------------------
% if (no_of_PC > 0)
%     [COEFF,SCORE,latent] = princomp(dataset_recall);
%     dataset_projected = SCORE(:,1:no_of_PC)*COEFF(:,1:no_of_PC)';
%     dataset_recall_orig = dataset_recall;
%     dataset_recall = dataset_projected;
% end
y_min = min(min(dataset_recall));
y_max = max(max(dataset_recall));
%--------------------------------------------------------------------------

%--------------------------Quantize if Necessary---------------------------
% if (Q > 0)    
%     dataset_recall = round(Q*tanh(sigma1*dataset_recall-mean1));    
% end
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
thr1 = max(max(abs(dataset_recall)));
% dataset_recall = [thr1*ones(dataset_size,1),dataset_recall];
%--------------------------------------------------------------------------


%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 

error_count = 0;
bit_error_count = 0;
error_count_target = 0;
bit_error_count_target = 0;
learning_error = zeros(1,dataset_size);
recall_error = zeros(1,dataset_size);
error_correction_gain = 0;
%--------------------------------------------------------------------------

%--------------------Create the Sub-folder If Necessary--------------------
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


destination_folder_recall = [db_file_in(1:i-1),'/Recall_Results/Clustered_Version'];
destination_folder_learn = [db_file_in(1:i-1),'/Learn_Results/Clustered_Version'];
if (~exist(destination_folder_recall,'dir'))
    mkdir(destination_folder_recall);
end
%--------------------------------------------------------------------------

%==========================================================================



%%
%===============================MAIN LOOP==================================
for net_simul_itr = 1:no_simulated_instances                % Simulate the error correction procedure for the given ensemble.                                
        
    %-----------------------------Generate Noise---------------------------
    nois = zeros(1,N);                                  % This is the noise added to the whole pattern of length N*L_in        
    pp = 1+floor((N-1)*rand(1,err_bits));                
    for h = 1:err_bits        
        nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
    end        
    %----------------------------------------------------------------------                            
                
    %-----------------------Generate the Pattern---------------------------
    mu = 1+floor((dataset_size-1)*rand);            % Pick a pattern index at random                   
%     mu = 9960;
    pattern = dataset_recall(mu,:);    
    projected_pattern = dataset_recall_projected(mu,:);
    learning_error(mu) = sum(abs((pattern-projected_pattern)));    
    %----------------------------------------------------------------------                                
                                                                                                                
    %------------Initialize the Network with a Noisy SubPattern------------
    x = pattern+nois;                                       % Initialize the network with a nosiy version of the subpattern               
    x = min(max(x,y_min),y_max);
    nois = x - pattern;
    %----------------------------------------------------------------------
            
        
    %------------------------Perform the Recall Step-----------------------
    success_flag = 0;    
    try_itr = 0;
    
    W_global_tot = [];
    W_total_tot = [];
    index_pattern_tot = zeros(no_clusters,cluster_size);
    converged_constraints_tot = zeros(1,no_clusters);
    for l = 1:no_clusters                                
        file_name = ['Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(l),'.mat'];
        fid = fopen(file_name,'r');
        if (fid > -1)
            fclose(fid);            
            load(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(l)]);            
        else
            continue;
        end
        save_flag = save_place(end,2);                        
        if (save_flag == 2)                           
            destination_folder_save = [destination_folder_learn,'/Partial_Convergence'];               
        else            
            destination_folder_save = destination_folder_learn;
        end        
        %-----------------------Store the Connectivity Matrix----------------------
        file_name = [destination_folder_save,'/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...            
        num2str(theta0),'_clustere_',num2str(l),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'];            
        fid = fopen(file_name,'r');        
        if (fid > -1)      
            fclose(fid);
            load(file_name)
            W_global_tot = [W_global_tot;W_global];
            W_total_tot = [W_total_tot;W_total];
            [mi,~] = size(W_total);
            index_pattern_tot(l,:) = index_pattern_neurons;
            converged_constraints_tot(l) = mi;
        end
    end
    
    counter = 0;
    while ((try_itr < try_max) && (success_flag == 0))
        success_flag = 1;
        try_itr = try_itr+1;
        consts_so_far = 1;
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
            
            
            x_temp = [thr1,x(index_pattern_neurons)];                                    
            pattern_temp = [thr1,pattern(index_pattern_neurons)]; 
            nois_effective = (x_temp - pattern_temp);
            projected_pattern_temp =[thr1,projected_pattern(index_pattern_neurons)];
            if  (sum((abs(nois_effective)).*sum(abs(sign(W_total)))) > 0)
                111;
            end
            %--------------------------------------------------------------    
            
            %---------------------Iterate Until Convergence----------------            
%             [x_out] = recall_step_clustered_real_value(W,x_temp,.05,y_min,y_max,0);    
            [x_out,update_term] = recall_step_clustered_v2(W_total,x_temp,varphi,psi,y_min,y_max,Q);
            %--------------------------------------------------------------
            
            %----------Verify the Success of the Recall Algorithm----------
%             if (norm(W*x_out')< norm(W*x_temp'))       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                       
            if ( (norm(W_total*x_out') < norm(W_total*x_temp')) )
                if (norm(x_out-pattern_temp) > 0)
                    111;
                end
                x(index_pattern_neurons) = x_out(2:end);                                 
%                 if ( (norm(x_out-x_temp) > .01) && (~update_term) )
%                     11;                
%                 end
%                 if (norm(W_total*x_out') > .05)
%                     success_flag = 0;                
%                 end
                
            end
            
            if (norm(W_total*[1,x(index_pattern_neurons)]') > .001)
                success_flag = 0;
            end
            
            %--------------------------------------------------------------
            
        end
    end
    %----------------------------------------------------------------------
                                                          
                
    %--------------------Calculate Recall Error Count----------------------
    x_target = dataset_recall_projected(mu,:);    
    recall_error(mu) = sum(abs(sign(x_target-x)));
    if (sum(recall_error(mu))>0)
        error_count_target = error_count_target+1;
        bit_error_count_target = bit_error_count_target+recall_error(mu);
    end
    
    if (norm(x-pattern)>.01)       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.
        error_count = error_count+1;
        bit_error_count = bit_error_count+sum(abs(x-pattern)>1e-1);
%         sum(abs(sign(x-pattern)));        
    else
        111;
    end
    %----------------------------------------------------------------------
    
    %--------------------------Display Progress----------------------------
    if (mod(net_simul_itr, 20) == 0)
        error_count/net_simul_itr
        bit_error_count/net_simul_itr/N
    end
     %----------------------------------------------------------------------

     error_correction_gain = error_correction_gain + sum(abs(nois))/bit_error_count;
    
    if (mod(net_simul_itr,100) == 0)
    %------------------Transform Error Count to Error Rate---------------------
    PER = error_count/net_simul_itr;
    BER = bit_error_count/net_simul_itr/N;
    error_correction_gain = error_correction_gain/net_simul_itr;
    %--------------------------------------------------------------------------


    
    %----------------Store the Bit and Pattern Error Rates-----------------
    fid = fopen([destination_folder_recall,'/n_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
           num2str(theta0),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'_',num2str(try_max),'.txt'], 'a+');  
    fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t gain \t %f \t',err_bits,PER,BER,error_correction_gain);
    fprintf(fid,'\n');
    fclose(fid);
    %----------------------------------------------------------------------
    end
    
end    
        
            
    


    
