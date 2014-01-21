%==========================================================================
%********************FUNCTION: clustered_neural_learn**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% no_clusters: The number of clusters if we have multiple levels (no_clusters = 1 for single level)
% const_learn: Number of constraints which must be learned during the learning phase
% cluster_index: The cluster for which the learning should be done
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% index: The index of the simulation setup among various random scenarios
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function implements the idea mentioned in our journal paper to learn 
% a sparse vector W which is orthogonal to a set of given patterns. 

% The code starts by reading patterns at random and adjusting the weights 
% accordingly such that in the end we have Wx = 0, for all data vectors x 
% in the training set. The learning approach which is very similar (almost 
% identical to) the one proposed in the paper "The subspace learning 
% algorithm as a formalism for pattern recognition and neural networks" 
% plus a penalty term to make the results sparse. 

% The code stops once a maximum number of iterations is passed or
% convergence is reached. Once finished, the learnt vectors are added to
% the end of the appropriatefile. 
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================

function clustered_neural_learn(no_clusters,cluster_size,cluster_index,const_learn,alpha0,beta0,theta0,db_file_in,db_name_in,learn_itr_max,Q,simulation_set,shift_horiz)
 

%%
%=============================INITIALIZATION===============================

save_place = [];
cluster_size
cluster_index
no_clusters
%---------------------------Load the Training Dataset----------------------
load(db_file_in);
eval(['dataset_learn = ',db_name_in,';']);   
[dataset_size,N] = size(dataset_learn);
%--------------------------------------------------------------------------


%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path

b=clock;                                                                        % Initialize the seed for random number generation with the clock value.
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*b))); 


max_y_threshold = .005;                                                          % The maximum acceptable projection for a pattern on a weight vector
min_W = .1;                                                                     % The minimum absolute value of the weight vectors
fraction_of_convergence = .995;                                                   % Required fraction of patterns to converge
req_frac_conv_const = .9;                                                       % Required number of constraints to converge for each pattern

index_pattern_neurons = [1+shift_horiz*(cluster_index-1):min(cluster_size+shift_horiz*(cluster_index-1),N)];

W_global = [];
W_total = [];
%--------------------------------------------------------------------------

%------------------------Compute the Projected Dataset---------------------
% if (no_of_PC > 0)
%     [COEFF,SCORE,latent] = princomp(dataset_learn);
%     dataset_projected = SCORE(:,1:no_of_PC)*COEFF(:,1:no_of_PC)';
%     dataset_learn = dataset_projected;
% end
%--------------------------------------------------------------------------

%--------------------------Quantize if Necessary---------------------------
% if (Q > 0)    
%     dataset_learn = dataset_learn-min(min(dataset_learn));
%     dataset_learn = dataset_learn/max(max(abs(dataset_learn)));
%     KK = ceil(log(Q)/log(2));
%     
%     dataset_learn = round(Q*dataset_learn);
%     vectorized_dataset = zeros(dataset_size,KK*N);
%     for i = 1:dataset_size
%         temp = dec2bin(dataset_learn(i,:),KK);
%         current_vector = [];
%         for j = 1:N
%             temp2 = temp(j,:);
%             message = zeros(1,KK);
%             for k = 1:KK    
%                 message(k) = temp2(k) - 48;
%             end
%             current_vector = [current_vector,message];
%         end
%         vectorized_dataset(i,:) = current_vector;
%         111;
%     end   
% %     max_y_threshold = max_y_threshold/Q;
% end
if (Q <= 1)
    sigma1 = 1;
    mean1 = 0;
end
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
% thr1 = max(max(abs(dataset_learn)));
thr1 = 1;
dataset_learn = [thr1*ones(dataset_size,1),dataset_learn];
dataset_learn_temp = dataset_learn(:,1+index_pattern_neurons);
dataset_learn_temp = [thr1*ones(dataset_size,1),dataset_learn_temp];
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


destination_folder = [db_file_in(1:i-1),'/Learn_Results/Clustered_Version'];
if (~exist(destination_folder,'dir'))
    mkdir(destination_folder);
end
%--------------------------------------------------------------------------

%-------------Update the Number of Constraints in the Process--------------
fid = fopen([destination_folder,'/job_being_process_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
    num2str(theta0),'_clustere_',num2str(cluster_index),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.txt'],'r');
if (fid > -1)
    r = fscanf(fid, '%d %d',[2,inf]);
    job_flag = r(1);                                                 
    job_id = r(2);
    fclose(fid);
    fid = fopen([destination_folder,'/job_being_process_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
        num2str(theta0),'_clustere_',num2str(cluster_index),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.txt'],'w');
    fprintf(fid, '%d %d', [1 job_id]);
    fclose(fid);
else
%     error('This cant be!')
end
%--------------------------------------------------------------------------

%==========================================================================
    


%%
%=============================LEARNING PHASE=============================== 
learned_consts = 0;
total_itr = 0;
total_itr_max = 2 * const_learn;

while ( (learned_consts <= const_learn)  && (total_itr < total_itr_max) )

    total_itr = total_itr + 1;
    
    %------------------Randomly Initialize the Weight Matrix---------------
    W = random_vector(cluster_size+1,round(1*(1+rand)*log(cluster_size)));
    % W = W';
%     W = rand(1,cluster_size+1);
    W = W./(sqrt(sum(W'.*W'))'*ones(1,cluster_size+1));    
    %----------------------------------------------------------------------
    
    %---------------------Adjust Simulation Parameters---------------------
    theta = theta0;       
    theta_itr = 1;    
    cost_per_itr_av = []; 
    no_converged_patts_per_itr = [];
    alph_itr = 0;        
    %----------------------------------------------------------------------

    %--------------------------Main Learning Loop--------------------------           
    for learn_itr = 1:learn_itr_max                         % Repeat the steps below several times in order to enhance the quality of the learning phase        
        
        
        if (mod(learn_itr,30) == 1)
            alph_itr = alph_itr + 1;
            alph0 = max(14*alpha0/(14+alph_itr),0.0005);        
        end        
                
        for ii = 1:dataset_size
            
        
            %------------------Pick a Pattern at Random--------------------
            mu = 1+floor((dataset_size-1)*rand); 
            x_temp = dataset_learn(mu,:);
            x = [thr1,x_temp(1+index_pattern_neurons)];                                    
            %--------------------------------------------------------------
        
            %---------------------Update the Weight Vector-----------------
            y = (W*x');                 % Calculate the projection of x on the vector W(j,:).
            
            if (norm(x)>0)
                alph = alph0/norm(x)^2;
            else
                alph = 0;
            end
            if (norm(W) > .001)                                                    
                W = soft_threshold(W - 1.05*alph*y*(x-y*W)+ .85*alph*(1-norm(W)^2)*W,theta);                                               
                W = W';
%                 W = W - alph*y*( x - (y*W/(norm(W))^2))-beta0*(soft_threshold_inv(W,theta))';
            end            
            %--------------------------------------------------------------
        
            
            %------------------Check for Numerical Errors------------------
            if (max(max(isnan(W)))>0)    
                W = rand(1,cluster_size+1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta0/(2+theta_itr);
                display('Error Nan');        
                break;                                    
            end
            
            if (norm(W) < .001)                
                display('Error zero!');                    
                W = rand(1,cluster_size+1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta0/(1+theta_itr);
            end
            %--------------------------------------------------------------
                                                                       
        end 
        
        %------------------Calcuate the Simulation Cost--------------------
        converged_count = sum(abs(W*dataset_learn_temp')<max_y_threshold);
        cost_one_itr = sum(abs(W*dataset_learn_temp'))/dataset_size;
        %------------------------------------------------------------------        
        
        %--------------------Update the Simulation Costs-------------------
        cost_per_itr_av = [cost_per_itr_av,cost_one_itr];
        no_converged_patts_per_itr = [no_converged_patts_per_itr,converged_count];
        ww = W;
        ww(~ww) = inf;       
        %------------------------------------------------------------------        
        

        %--------------------------Display Progress------------------------                
        if (mod(learn_itr,20)==0)            
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);
            display(['Converged count = ',num2str(converged_count)]);
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);
            display(['Average Norm-0 of W = ',num2str(mean(sum(abs(sign(W(2:end))))))]);
            display(['Maximum projection = ',num2str(max(max(abs(W*dataset_learn_temp'))))]);  
            display(['Minum of W = ',num2str(min(min(abs(ww))))]);  
            display(['Max deg of pattern neuron = ',num2str(max(sum(abs(sign(W_global(:,2:end))))))]);  
            display(['Number of connected pattern neurons = ',num2str(sum(sum(abs(sign(W_global(:,2:end))))>0))]);  
            display('-----------------------------------');
            111;
            display(' ');             
        end           
               
        if (mod(learn_itr,15) == 0)
            theta_itr = theta_itr + 1;
            theta = 3*theta0/(3+theta_itr);              
        end
        
        if (mod(learn_itr,85) == 0)
            W = soft_threshold(W,max_y_threshold+.0015);
            W = W';
        end
        %--------------------------------------------------------------
         
            

        %----------------------Check Convergence---------------------------
        ww = W;
        ww(~ww) = inf;

        if ( (converged_count>fraction_of_convergence*dataset_size) && (norm(W)>.1))
            W = soft_threshold(W,1*max_y_threshold+0.0015);
            W = W';
            converged_count = sum(abs(W*dataset_learn_temp')<max_y_threshold)
                            
            if ( (converged_count>fraction_of_convergence*dataset_size) && (norm(W)>.1))
                break;
            end
        end        
        %------------------------------------------------------------------        
        

    end  
    %======================================================================

    
    %============================SAVE THE RESULTS==========================    
    W = soft_threshold(W,1*max_y_threshold+0.0015);            
    W = W';
    converged_count = sum(abs(W*dataset_learn_temp')<max_y_threshold);
    if ( (converged_count>fraction_of_convergence*dataset_size) && (norm(W)>.1))
       save_flag = 1;
    else
        save_flag = 2;
    end
    
    if (save_flag == 2)            
        destination_folder_save = [destination_folder,'/Partial_Convergence'];            
        if (~exist(destination_folder_save,'dir'))                
            mkdir(destination_folder_save)
        end
        
    else        
        destination_folder_save = destination_folder;        
    end
       
    W_temp = W;            
    W2 = zeros(1,N+1);                
    W2(1) = W_temp(1);        
    W2(1+index_pattern_neurons) = W_temp(2:end);
    

    W_temp = [W_global;W2];
    
                   
    success_flag = 1;         
   [mi,~] = size(W_global);  
   for ii = 1:mi
       var1 = abs(soft_threshold(W_global(ii,:),0.05));
       var2 = abs(soft_threshold(W2,0.05));
       if ( (norm(var1-var2) < .1) )
           success_flag = 0;
           break;
       end
   end
   if (converged_count < fraction_of_convergence * dataset_size)
       success_flag = 0;
   end
   
   if (success_flag)
       W_global = [W_global;W2];
       W_total = [W_total;W];                
       learned_consts = learned_consts + 1;
       Learn_itr_cell{learned_consts} = learn_itr;
        
       save_place = [save_place;learned_consts,save_flag];
       save(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(cluster_index)],'save_place');
       %-----------------------Store the Connectivity Matrix----------------------
       save([destination_folder_save,'/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
           num2str(theta0),'_clustere_',num2str(cluster_index),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'],...
           'W_global','W_total','index_pattern_neurons','cost_per_itr_av','Learn_itr_cell','Q','sigma1','mean1');
       %--------------------------------------------------------------------------
        [mi,~] = size(W_global);
        avg_no_of_converged_const = sum(abs(W_global*dataset_learn')<max_y_threshold)/mi;
   end
  
   learned_consts
   total_itr
   %======================================================================
   
end

% w_thr = .05;
% % W_total2 = soft_threshold_matrix(W_total,w_thr + .005);
% fraction_converged = sum(abs(W_total2*dataset_learn_temp')<max_y_threshold)/mi;
% converged_count = sum(fraction_converged > req_frac_conv_const)

fid = fopen([destination_folder,'/job_being_process_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
    num2str(theta0),'_clustere_',num2str(cluster_index),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.txt'],'w');
    fprintf(fid, '%d %d', [2 job_id]);    
fclose(fid);
    



   