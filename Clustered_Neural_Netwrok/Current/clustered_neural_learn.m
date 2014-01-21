%==========================================================================
%********************FUNCTION: clustered_neural_learn**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
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

% function clustered_neural_learn(N,L,const_learn,cluster_index,alpha0,beta0,theta0,db_file_in,db_name_in)
% 
L = 800;

dataset_zero_thr = .075;
target_class = 6;

no_nonzero = 10;
dict_size = 80;
patch_size = 8;

db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/CIFAR_10_train_gray_class_',num2str(target_class),'_OMP_',num2str(no_nonzero),'_dict_',num2str(dict_size),'_patch_',num2str(patch_size),'_stand.mat'];
db_name_in = 'rec_ims';

db_file_in = '/scratch/amir/Databases/CIFAR_10/Color_Class_6/Train_Set/OMP_5/Pooled/pooled_qudrant_class_6_omp_8_50_2500.mat';
db_name_in = 'pooled_features';

db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Mixed/CIFAR_10_Gray_Mixed_DB.mat'];  
db_name_in = 'CIFAR_10_Gray_Mixed_DB';

% db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/Sparse Filtered/Layer_2/Pooled_3_3/Pooled_Sparse_Filtered_CIFAR_10_Gray_Class_',num2str(target_class),'_Train_Layer_2.mat'];
% db_name_in = 'final_database_vectorized';

% db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/CIFAR_10_train_gray_scale_class_',num2str(target_class),'.mat'];
% db_name_in = ['CIFAR_10_Gray_DB_class_',num2str(target_class)];

% db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/Whitened/CIFAR_10_train_gray_scale_class_',num2str(target_class),'_whitened.mat'];
% db_name_in = ['CIFAR_10_Gray_DB_class_',num2str(target_class),'_whitened'];

%%

%=============================INITIALIZATION===============================

%---------------------------Load the Training Dataset----------------------
load(db_file_in);
eval(['dataset_learn = ',db_name_in,';']);   
% dataset_learn = dataset_learn';

dataset_learn = round(5*tanh(.7*dataset_learn));

[dataset_size,pattern_length] = size(dataset_learn);
N = pattern_length;
% dataset_learn = dataset_learn .*(abs(dataset_learn)>dataset_zero_thr);
% dataset_learn  = ones(dataset_size,pattern_length).*((dataset_learn>dataset_zero_thr)-(dataset_learn<-dataset_zero_thr));
% dataset_learn = abs(dataset_learn);
% a = sqrt(sum(dataset_learn'.*dataset_learn'));
% dataset_learn = dataset_learn./(a'*ones(1,1024));
% dataset_learn = dataset_learn.*(abs(dataset_learn)>.03);
%--------------------------------------------------------------------------
 
%--------------Add a Column to the Dataset for the Threshold---------------
dataset_learn = [ones(dataset_size,1),dataset_learn];
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


destination_folder = [db_file_in(1:i-1),'/Learn_Results/Zero_One/Clustered_Version/N_',num2str(N),'L_',num2str(L)];
if (~exist(destination_folder,'dir'))
    mkdir(destination_folder);
end
%--------------------------------------------------------------------------

%==========================================================================

    
temp = randperm(N);
index_pattern_neurons = temp(1:100);
const_learn = 25;
alpha0=.5;
beta0=.75;
theta0=.005;
   
%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
dataset_zero_thr = .05;                                                         % This is the threshold below which we set entries in the dataset to zero

b=clock;                                                                        % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*b))); 

learn_itr_max = 2500;                                                           % The maximum number of learning iterations
max_y_threshold = .08;                                                          % The maximum acceptable projection for a pattern on a weight vector
min_W = .1;                                                                     % The minimum absolute value of the weight vectors
fraction_of_convergence = .9;                                                   % Required fraction of patterns to converge
req_frac_conv_const = .9;                                                       % Required number of constraints to converge for each pattern
n = length(index_pattern_neurons);                                              % n is the total number of pattern neurons in the cluster
%--------------------------------------------------------------------------


%-------------Update the Number of Constraints in the Process--------------
fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'r');
if (fid > -1)           
    W = fscanf(fid, '%f',[N+1,inf]);                
    W = W';            
    fclose(fid);                                    
    [m,~] = size(W);            
else    
    m = 0;    
end
%--------------------------------------------------------------------------

%==========================================================================
    

%%
%=============================LEARNING PHASE===============================
         
%--------------------Randomly Initialize the Weight Matrix-----------------
% W = random_vector(n+1,2*round(log(n)));
% W = W';
W = rand(const_learn,n+1)-.5;
W = W./(sqrt(sum(W'.*W'))'*ones(1,n+1));    
%--------------------------------------------------------------------------
    
%-----------------------Adjust Simulation Parameters-----------------------
theta = theta0;       
theta_itr = 1;    
cost_per_itr_av = [];                
alph_itr = 0;        
%--------------------------------------------------------------------------

    %--------------------------Main Learning Loop--------------------------           
    for learn_itr = 1:learn_itr_max                         % Repeat the steps below several times in order to enhance the quality of the learning phase        
        max_y = 0;
        
        if (mod(learn_itr,512) == 1)            
            alph_itr = alph_itr + 1;
            alph0 = max(alpha0/alph_itr,0.0005);        
        end        
        
        cost_one_itr = 0;        
        cost_short_term =0;
        for ii = 1:dataset_size
            
        
            %------------------Pick a Pattern at Random--------------------
            mu = 1+floor((dataset_size-1)*rand); 
            x_temp = dataset_learn(mu,:);
            x = [1];                        
            for i = 1:length(index_pattern_neurons)
                x = [x,x_temp(1+index_pattern_neurons(i))];                                
            end    
            
            %--------------------------------------------------------------
        
            %---------------------Update the Weight Vector-----------------
            y = tanh(W*x');                 % Calculate the projection of x on the vector W(j,:).
            
            if (norm(x)>0)
                alph = alph0/norm(x)^2;
            else
                alph = 0;
            end
            if (norm(W) > .001)                                
%                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,0.05);W.*(1-(tanh(100*W(j,:).^2)).^2)
                for ij = 1:const_learn 
                    W(ij,:) = soft_threshold(W(ij,:) - alph*y(ij)*( x - (y(ij)*(W(ij,:)/(norm(W(ij,:))^2)))),theta);
%                     W(ij,:) = soft_threshold(W(ij,:) - alph*y(ij)*(x*(2*norm(W(ij,:))^2-1) - y(ij)*W(ij,:)),theta);
%                     W(ij,:) = soft_threshold(W(ij,:) - alph*y(ij)*(x-y(ij)*W(ij,:))+ alph*(1-norm(W(ij,:))^2)*W(ij,:),theta);                           
                end
                % W = soft_threshold_matrix(W - alph*y*( x - (y'*(W./(sqrt(sum(W'.*W'))'*ones(1,n+1))))),theta0);%-beta0*W.*(1-(tanh(100*W.^2)).^2);                    

            end            
            %--------------------------------------------------------------
        
            
            %------------------Check for Numerical Errors------------------
            if (max(max(isnan(W)))>0)    
                W = rand(n+1,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
                error('Error Nan');        
                break;                                    
            end
            
            if (norm(W) < .001)                
                error('Error zero!');                    
                W = rand(n+1,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
            end
            %--------------------------------------------------------------

            %-----------------Update the Simulation Costs------------------
            cost_one_itr = cost_one_itr + mean(abs(W*x'));            
            cost_short_term = cost_short_term+mean(abs(W*x'));
            
            if (norm(x)>0)
                cos_phi = abs(y)/(norm(W)*norm(x));
            else
                cos_phi = 1;
            end
            
            if (max(abs(W*x'))>max_y)
                max_y = max(abs(W*x'));
            end                                                      
            %--------------------------------------------------------------
                                                           
            
            
        end 
                
        %--------------------Update the Simulation Costs-------------------
        cost_per_itr_av = [cost_per_itr_av,cost_short_term];
        cost_short_term = 0;
        avg_no_of_converged_const = sum(abs(W2*dataset_learn')<max_y_threshold)/const_learn;
        %------------------------------------------------------------------                
        
        
        %--------------Map the Weight Vector to Actual Size----------------
        W_temp = W;
        W2 = zeros(const_learn,N+1);        
        W2(1) = W_temp(1);
        for i = 1:length(index_pattern_neurons)                
            W2(:,1+index_pattern_neurons(i)) = W_temp(:,i+1);            
        end       
        
        if (sum(max(abs(W2*dataset_learn'))<max_y_threshold)>fraction_of_convergence*dataset_size)
            W2=soft_threshold_matrix(W2,min_W);         
        end
        
        ww = W2;
        ww(~ww) = inf;
        %------------------------------------------------------------------
        
        %--------------------------Display Progress------------------------
        
        
        if (mod(learn_itr,10)==0)            
            
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);
            display(['Converged count = ',num2str(sum(avg_no_of_converged_const>req_frac_conv_const))]);
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);
            display(['Average Norm-0 of W = ',num2str(mean(sum(abs(sign(W(:,2:end))))))]);
            display(['Maximum projection = ',num2str(max(max(abs(W2*dataset_learn'))))]);  
            display(['Minum of W = ',num2str(min(min(abs(ww))))]);  
            display(['Max deg of pattern neuron = ',num2str(max(sum(abs(sign(W(:,2:end))))))]);  
            display(['Number of connected pattern neurons = ',num2str(sum(sum(abs(sign(W(:,2:end))))>0))]);  
            display('-----------------------------------');
            111;
            display(' ');             
        end           
        
        
        if (mod(learn_itr,120) == 0)
            theta_itr = theta_itr + 1;
            theta = 1*theta1/(1+theta_itr);            
%             if (sum((abs(dataset_learn*W2)<max_y_threshold))<dataset_size/2)
%                 W = rand(n,1);
%                 W = W/norm(W);
%             end
        end
        %--------------------------------------------------------------
         
            

        %----------------------Check Convergence---------------------------
        ww = W;
        ww(~ww) = inf;

        if ( (sum(avg_no_of_converged_const>req_frac_conv_const)>fraction_of_convergence*dataset_size) &&(min(min(abs(ww)))>min_W)&&(norm(W)>.1))
            break;
        end
        %------------------------------------------------------------------        
        

    end  
    %======================================================================


    
%============================SAVE THE RESULTS==============================

%-----------------------Store the Connectivity Matrix----------------------
W_temp = W;            
W2 = zeros(const_learn,N+1);                
W2(1) = W_temp(1);        
for i = 1:length(index_pattern_neurons)                            
    W2(:,1+index_pattern_neurons(i)) = W_temp(:,i+1);                        
end

fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'w');
for i = 1:const_learn        
    fprintf(fid, '%f \t',W2(i,:));                
    fprintf(fid, '\n');            
end
fclose(fid);
%--------------------------------------------------------------------------

%--------------Keep Track of Patterns Learned in This Cluster--------------
converged_flag = (max(abs(W2*dataset_learn')<max_y_threshold));
%--------------------------------------------------------------------------    
    
%---------------------Remove Unconnected Pattern Neurons-------------------    
s = sum(abs(sign(W(:,2:end))));
temp = [];
for j = 1:n
    if (s(j) > 0)
        temp = [temp,index_pattern_neurons(j)];
    end
end
index_pattern_neurons = temp;
%--------------------------------------------------------------------------
    
    save([destination_folder,'/Simulation_Flags_And_Learned_Pattern_Cluster_',num2str(cluster_index),'.mat'],'converged_flag','index_pattern_neurons','W','cost_per_itr_av');
    
    fid = fopen([destination_folder,'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'a');
    fprintf(fid, '%d \t',learn_itr);
    fprintf(fid, '\n');
    fclose(fid);    
        
    
    
%==========================================================================


   