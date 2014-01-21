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

% function clustered_neural_learn(no_clusters,cluster_size,cluster_index,const_learn,alpha0,beta0,theta0,db_file_in,db_name_in,learn_itr_max,Q,sigma1,mean1,simulation_set,no_of_PC)
 

%%
%=============================INITIALIZATION===============================

%---------------------------Load the Training Dataset----------------------
load(db_file_in);
eval(['dataset_learn = ',db_name_in,';']);   
[dataset_size,N] = size(dataset_learn);
%--------------------------------------------------------------------------


%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path

b=clock;                                                                        % Initialize the seed for random number generation with the clock value.
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*b))); 


max_y_threshold = .01;                                                          % The maximum acceptable projection for a pattern on a weight vector
min_W = .1;                                                                     % The minimum absolute value of the weight vectors
fraction_of_convergence = .9;                                                   % Required fraction of patterns to converge
req_frac_conv_const = .99;                                                       % Required number of constraints to converge for each pattern
%--------------------------------------------------------------------------

%------------------------Compute the Projected Dataset---------------------
if (no_of_PC > 0)
    [COEFF,SCORE,latent] = princomp(dataset_learn);
    dataset_projected = SCORE(:,1:no_of_PC)*COEFF(:,1:no_of_PC)';
    dataset_learn = dataset_projected;
end
%--------------------------------------------------------------------------

%--------------------------Quantize if Necessary---------------------------
if (Q > 0)    
    dataset_learn = round(Q*tanh(sigma1*dataset_learn-mean1));
    dataset_learn = dataset_learn;
    max_y_threshold = max_y_threshold;
end
%--------------------------------------------------------------------------

%--------------Add a Column to the Dataset for the Threshold---------------
thr1 = max(max(abs(dataset_learn)));
dataset_learn = [thr1*ones(dataset_size,1),dataset_learn];
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
total_itr_max = 10;
while ( (total_itr < total_itr_max) )
    
    %-----------------------In-Loop Initializations------------------------
    temp = randperm(N);
    index_pattern_neurons = temp(1:cluster_size);

    W_global = [];
    W_total = [];
    total_itr = total_itr + 1;
    
    dataset_learn_temp = dataset_learn(:,1+index_pattern_neurons);
    dataset_learn_temp = [thr1*ones(dataset_size,1),dataset_learn_temp];
    %----------------------------------------------------------------------
        
    %------------------Randomly Initialize the Weight Matrix---------------
    W = [];
    for i = 1:const_learn
        W = [W;random_vector(cluster_size+1,2*round(log(cluster_size)))];
    end
%     W = rand(const_learn,cluster_size+1);
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
        learn_itr      
        if (mod(learn_itr,420) == 1)            
            alph_itr = alph_itr + 1;
            alph0 = max(4*alpha0/(3+alph_itr),0.0005);        
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
                W = soft_threshold_matrix(W - alph*(y*x-diag(ones(const_learn,1)-(sum(W'.*W'))')*W),theta);%-beta0*W.*(1-(tanh(100*W.^2)).^2);                    
            end            
            %--------------------------------------------------------------
        
            
            %------------------Check for Numerical Errors------------------
            if (max(max(isnan(W)))>0)    
                W = rand(cluster_size+1,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
                error('Error Nan');        
                break;                                    
            end
            
            if (norm(W) < .001)                
                error('Error zero!');                    
                W = rand(cluster_size+1,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta0/(2+theta_itr);
            end
            %--------------------------------------------------------------
                                                           
        end 
        
        %-------------Calculate Number of Converged Patterns---------------                        
        no_converged_consts = sum(abs(W*dataset_learn_temp')<max_y_threshold);
        converged_count = sum(no_converged_consts > req_frac_conv_const * const_learn);
        cost_one_itr = sum(sum((W*dataset_learn_temp').*(W*dataset_learn_temp')));
        %------------------------------------------------------------------        
        
        %--------------------Update the Simulation Costs-------------------
        cost_per_itr_av = [cost_per_itr_av,cost_one_itr];
        no_converged_patts_per_itr = [no_converged_patts_per_itr,converged_count];        
        ww = W;
        ww(~ww) = inf;
        %------------------------------------------------------------------        
                

        %--------------------------Display Progress------------------------
        if (mod(learn_itr,40)==0)
            
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);
            display(['Converged count = ',num2str(converged_count)]);
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);            
            display(['Maximum projection = ',num2str(max(max(abs(W*dataset_learn_temp'))))]);  
            display(['Minum of W = ',num2str(min(min(abs(ww))))]);  
            display(['Max deg of pattern neuron = ',num2str(max(sum(abs(sign(W(:,2:end))))))]);  
            display(['Number of connected pattern neurons = ',num2str(sum(sum(abs(sign(W(:,2:end))))>0))]);  
            display('-----------------------------------');
            111;
            display(' ');             
        end           
        
        
        if (mod(learn_itr,10) == 0)
            theta_itr = theta_itr + 1;
            theta = 1*theta0/(1+theta_itr);              
        end
        
        if (mod(learn_itr,100) == 0)
            W = soft_threshold_matrix(W,max_y_threshold+.005);
            
        end
        %--------------------------------------------------------------
         
            

        %----------------------Check Convergence---------------------------        
        if ( (converged_count>fraction_of_convergence*dataset_size) && (norm(W)>.1))
            W = soft_threshold_matrix(W,max_y_threshold+0.005);
            no_converged_consts = sum(abs(W*dataset_learn_temp')<max_y_threshold);
            converged_count = sum(no_converged_consts > req_frac_conv_const * const_learn);
            cost_one_itr = sum(sum((W*dataset_learn_temp').*(W*dataset_learn_temp')));

            if ( (converged_count>fraction_of_convergence*dataset_size) && (norm(W)>.1))
                break;
            end
        end        
       %------------------------------------------------------------------        
        

    end  
    %======================================================================
    
    
    %============================SAVE THE RESULTS==========================
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
       
    
    for i = 1:const_learn       
        W_temp = W(i,:);
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
        if (success_flag)
            W_global = [W_global;W2];
            W_total = [W_total;W_temp];                                        
        end
    end
           
    [mi,~] = size(W_global);
    
    if (mi > .75*const_learn)              
        save(['./Progressive_Simulation_Results/Simulation_results_set_',num2str(simulation_set),'_cluster_',num2str(cluster_index)],'save_flag');            
        %-----------------------Store the Connectivity Matrix----------------------            
        save([destination_folder_save,'/Simulation_Results_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...                
            num2str(theta0),'_clustere_',num2str(cluster_index),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.mat'],...        
            'W_global','W_total','index_pattern_neurons','cost_per_itr_av','Learn_itr','Q','sigma1','mean1');        
        %--------------------------------------------------------------------------   
        break;
    end    
    %======================================================================
    
end

fid = fopen([destination_folder,'/job_being_process_cluster_size_',num2str(cluster_size),'_consts_',num2str(const_learn),'_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',...
        num2str(theta0),'_clustere_',num2str(cluster_index),'_Q_',num2str(Q),'_itr_',num2str(learn_itr_max),'.txt'],'w');
    
fprintf(fid, '%d %d', [2 job_id]);    
fclose(fid);
    