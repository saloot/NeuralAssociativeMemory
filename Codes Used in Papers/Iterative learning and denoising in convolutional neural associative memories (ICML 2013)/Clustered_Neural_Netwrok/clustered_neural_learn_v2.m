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
W_global = [];
N = 1024;
L = 800;
L2 = 1200;
dataset_zero_thr = .075;
target_class = 1;

no_nonzero = 7;
dict_size = 50;
patch_size = 8;

db_file_in = ['/scratch/amir/Databases/CIFAR_10/Gray_Scale_Class_',num2str(target_class),'/Train_Set/CIFAR_10_train_gray_class_',num2str(target_class),'_OMP_',num2str(no_nonzero),'_dict_',num2str(dict_size),'_patch_',num2str(patch_size),'_stand.mat'];
db_name_in = 'rec_ims';

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
dataset_learn = dataset_learn';

[dataset_size,pattern_length] = size(dataset_learn);
N = pattern_length;
% dataset_learn = dataset_learn .*(abs(dataset_learn)>dataset_zero_thr);
% dataset_learn  = ones(dataset_size,pattern_length).*((dataset_learn>dataset_zero_thr)-(dataset_learn<-dataset_zero_thr));
% dataset_learn = abs(dataset_learn);
% a = sqrt(sum(dataset_learn'.*dataset_learn'));
% dataset_learn = dataset_learn./(a'*ones(1,1024));
% dataset_learn = dataset_learn.*(abs(dataset_learn)>.03);
%--------------------------------------------------------------------------
% 
% %------------------------------Quantization--------------------------------
% % partition = [0.0001,.00015,.0002,.00025,.0003,.00035,.0004,.00045,.0005,.0006,.0007,.0008,.0009,.001,.002,.004,.006,.008,.01];
% % partition = [0.0001,.00015,.0002,.00025,.0003,.00035,.0004,.00045,.0005,.0006,.0007,.0008,.0009,.001,.002,.0025,.003,.0035,.004,.0045,.005,.0055,.006,.0065,.007,.0075,.008,.0085,.009,.0095,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1];
% % partition = [0.0001,.00015,.0002,.00025,.0003,.00035,.0004,.00045,.0005,.0006,.0007,.0008,.0009,.001,.002,.0025,.003,.0035,.004,.0045,.005,.0055,.006,.0065,.007,.0075,.008,.0085,.009,.0095,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.12,.14,.16,.18]
% partition = [0:.2:1];
% codebook = zeros(1,length(partition)+1);
% codebook(2:length(partition)+1) = partition;
% db_temp2 = [];
% for i = 1:dataset_size
%     [index,quants] = quantiz(dataset_learn(i,:),partition,codebook);    
%     db_temp2 = [db_temp2;index];        
% end
% %--------------------------------------------------------------------------
% 
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


for cluster_index = 1:L2

    
temp = randperm(N);
index_pattern_neurons = temp(1:50);
alpha0=.95;
beta0=.75;
theta0=.01;

const_learn = 1;
    
% %----------------Load the Saved Initialized Parameters---------------------
% load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_L_',num2str(L2),...
%     '/Clustered_parameters_N_',num2str(N),'_L_',num2str(L2),'_cluster_index_',num2str(cluster_index),'.mat']);           
% %--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
dataset_zero_thr = .05;                                                         % This is the threshold below which we set entries in the dataset to zero

b=clock;                                                                        % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*b))); 

learn_itr_max = 500;                                                           % The maximum number of learning iterations
max_y_threshold = .024;                                                          % The maximum acceptable projection for a pattern on a weight vector
min_W = .025;                                                                   % The minimum absolute value of the weight vectors
fraction_of_convergence = .999;                                                  % Required fraction of patterns to converge
d_max = 22;                                                                     % The maximum allowed degree for pattern neurons
n = length(index_pattern_neurons);                                              % n is the total number of pattern neurons in the cluster
%--------------------------------------------------------------------------




%-------------Update the Number of Constraints in the Process--------------
fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'r'); 
if (fid > -1)
    current_count = max(fscanf(fid, '%d'),0);            
    fclose(fid);    
else    
    current_count = 0;            
end        

fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'r');
if (fid > -1)           
    W = fscanf(fid, '%f',[n+1,inf]);                
    W = W';            
    fclose(fid);                                    
    [m,~] = size(W);            
else    
    m = 0;    
end

const_to_learn = min(const_learn,N-m-current_count);

fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'w'); 
fprintf(fid,'%d',current_count+const_to_learn);
fclose(fid);
%--------------------------------------------------------------------------

%==========================================================================
    

%%
%=============================LEARNING PHASE===============================
const_this_round = 0;
for itr = 1: const_to_learn

    alpha1 = alpha0;
    theta1 = theta0;
    %-----------------------Read Current Learned Matrix--------------------
%     fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'r');
%     if (fid > -1)           
%         W_tot = fscanf(fid, '%f',[n+1,inf]);                
%         W_tot = W_tot';            
%         fclose(fid);                                    
%         [m,~] = size(W_tot);            
%     else    
%         m = 0;    
%     end
    [m,~] = size(W_global);
    if (m >= 2)
       deg_tot = sum(abs(sign(W_global)));
    else
       deg_tot = zeros(1,N+1);
    end                        
    %----------------------------------------------------------------------
       
    %----------------------------------------------------------------------

    %------------------Randomly Initialize the Weight Matrix---------------
    W = random_vector(n+1,2*round(log(n)));
    W = W';
    W = rand(n+1,1);
    W = W/norm(W);
    %----------------------------------------------------------------------
    
    %---------------------Adjust Simulation Parameters--------------------- 
    theta = theta1;   
    theta_itr = 1;
    cost_per_itr_av = [];            
    alph_itr = 0;    
    %----------------------------------------------------------------------

    %--------------------------Main Learning Loop--------------------------           
    for learn_itr = 1:learn_itr_max                         % Repeat the steps below several times in order to enhance the quality of the learning phase        
%         if (mod(learn_itr,round(log(2+const_this_round)*25))==0)
% %             W = random_vector(N,2*round(log(N)));
% %             W = W';
%             W = rand(n+1,1);
%             W = W/norm(W);
%             theta_itr = theta_itr + 1;                
%             theta = theta1/(2+theta_itr);
%             alpha1 = alpha0/theta_itr;
%         end

        max_y = 0;
        
%         alph0 = max(50*alpha1/(50+(learn_itr)),0.0005);
        if (mod(learn_itr,202) == 1)            
            alph_itr = alph_itr + 1;
            alph0 = max(alpha1/alph_itr,0.0005);        
        end        
        
        cost_one_itr = 0;        
        cost_short_term =0;
        for ii = 1:dataset_size
            
            %------------------Pick a Pattern at Random--------------------
            mu = 1+floor((dataset_size-1)*rand); 
            x_temp = dataset_learn(mu,:);
            x = [1];            
            deg = [0];
            for i = 1:length(index_pattern_neurons)
                x = [x,x_temp(1+index_pattern_neurons(i))];
                deg = [deg,deg_tot(1+index_pattern_neurons(i))];
                
            end    
            
            %--------------------------------------------------------------
        
            %---------------------Update the Weight Vector-----------------
            y = tanh(x*W);                 % Calculate the projection of x on the vector W(j,:).
            
            if (norm(x)>0)
                alph = alph0/norm(x)^2;
            else
                alph = 0;
            end
            if (norm(W) > .001)                                
%                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,0.05);W.*(1-(tanh(100*W(j,:).^2)).^2)
%                  W = soft_threshold(W - alph*y*( x' - (y*W/(norm(W))^2)),theta0);%-beta0*W.*(1-(tanh(100*W.^2)).^2);
%                  W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
%                  W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta)-.05*W.*(1-(tanh(100*W.^2)).^2).*(1+tanh(100*(deg-d_max)))';
                W = soft_threshold(W - alph*y*(x'-y*W)+ alph*(1-norm(W)^2)*W-.025*W.*(1-(tanh(100*W.^2)).^2).*((1+tanh(.3*(deg-d_max))')),theta); 
%                  
            end            
            
%             W = soft_threshold(W,theta);                       % Set entries less than 0.0001 to zero.
            %--------------------------------------------------------------
        
            
            %------------------Check for Numerical Errors------------------
            if (sum(isnan(W))>0)    
                W = rand(n+1,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
                display('Error Nan');        
                break;                                    
            end
            
            if (norm(W) < .001)                
                display('Error zero!');                    
                W = rand(n+1,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
            end
            %--------------------------------------------------------------

            %-----------------Update the Simulation Costs------------------
            cost_one_itr = cost_one_itr + abs(x*W);
            cost_short_term = cost_short_term+norm(x*W);
            
            if (norm(x)>0)
                cos_phi = abs(y)/(norm(W)*norm(x));
            else
                cos_phi = 1;
            end
            
            if (abs(x*W)>max_y)
                max_y = abs(x*W);
            end                                                      
            %--------------------------------------------------------------
                                                           
            
            
        end 
                
        %--------------------Update the Simulation Costs-------------------
        cost_per_itr_av = [cost_per_itr_av,cost_short_term];
        cost_short_term = 0;
        %------------------------------------------------------------------                
        
        
        %--------------Map the Weight Vector to Actual Size----------------
        W_temp = W;
        W2 = zeros(N+1,1);        
        W2(1) = W_temp(1);
        for i = 1:length(index_pattern_neurons)                
            W2(1+index_pattern_neurons(i)) = W_temp(i+1);            
        end       
        
        if (sum(abs(dataset_learn*W2)<max_y_threshold)>fraction_of_convergence*dataset_size)
            W2=soft_threshold(W2,min_W);         
        end
        
        ww = W2;
        ww(~ww) = inf;
        %------------------------------------------------------------------
        
        %--------------------------Display Progress------------------------
        if (mod(learn_itr,10)==0)            
            
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);
            display(['Converged count = ',num2str(sum((abs(dataset_learn*W2)<max_y_threshold)))]);
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);
            display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
            display(['Maximum projection = ',num2str(max(abs(dataset_learn*W2)))]);  
            display(['Minum of W = ',num2str(min(abs(ww)))]);  
            display('-----------------------------------');
            111;
            display(' ');             
        end           
        
        
        if (mod(learn_itr,4) == 0)
            theta_itr = theta_itr + 1;
            theta = 4*theta1/(4+theta_itr);            
%             if (sum((abs(dataset_learn*W2)<max_y_threshold))<dataset_size/2)
%                 W = rand(n,1);
%                 W = W/norm(W);
%             end
        end
        %--------------------------------------------------------------
         
            

        %----------------------Check Convergence---------------------------
        ww = W;
        ww(~ww) = inf;

        if ((sum((abs(dataset_learn*W2))<max_y_threshold)>=fraction_of_convergence*dataset_size)&&(min(abs(ww))>min_W)&&(norm(W)>.1))
            break;
        end
        %------------------------------------------------------------------        
        

    end  
    %======================================================================


    %==========================SAVE THE RESULTS============================    
    if ((sum((abs(dataset_learn*W2))<max_y_threshold)>=fraction_of_convergence*dataset_size)&&(min(abs(ww))>min_W)&&(norm(W)>.1))                                
        const_this_round = const_this_round+1;
        display(['-------Convergence!!!--------']);        
        display(['Converged count = ',num2str(sum((abs(dataset_learn*W2)<max_y_threshold)))]);        
        display(['Cost in last iteration = ',num2str(cost_one_itr)]);
        display(['Norm-2 of W = ',num2str(norm(W))]);
        display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
        display(['Maximum projection = ',num2str(max(abs(dataset_learn*W2)))]);  
        display(['Minum of W = ',num2str(min(abs(ww)))]);  
        display('-----------------------------------');
        display(' ');   
        
        fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'a');
        fprintf(fid, '%f \t',W');
        fprintf(fid, '\n');
        fclose(fid);
                
            
%         W_temp = [W_global;W2'];
        
        
%         if (rank(sign(W_temp))>rank(sign(W_global)))
         success_flag = 1;
         [mi,~] = size(W_global);
         for i = 1:mi
             if ((norm(sign(W2')-sign(W_global(mi,:)))<.1) || (norm(sign(W2')+sign(W_global(mi,:)))<.1))
                 success_flag = 0;
                 break;
             end
         end
         if (success_flag)
            W_global = [W_global;W2'];
            fid = fopen([destination_folder,'/Weigh_matrix_global_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'a');
            fprintf(fid, '%f \t',W2');
            fprintf(fid, '\n');
            fclose(fid);
        end
        
        fid = fopen([destination_folder,'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'a');
        fprintf(fid, '%d \t',learn_itr);
        fprintf(fid, '\n');
        fclose(fid);    
        
        fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'r'); 
        if (fid > -1)
            current_count = fscanf(fid, '%d');
            fclose(fid);
        else
            current_count = 0;
            display('error!');
        end
    
        fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'clustered_index_',num2str(cluster_index),'.txt'], 'w'); 
        fprintf(fid,'%d',max(current_count-1,0));
        fclose(fid);
    end
    %======================================================================
    
end
sum(sign(sum(abs(sign(W_global)))))
rank(sign(W_global))
size(W_global)
end
   