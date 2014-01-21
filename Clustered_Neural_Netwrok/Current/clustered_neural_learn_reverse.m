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

function clustered_neural_learn_reverse(N,const_learn,alpha0,beta0,theta0,db_file_in,db_name_in)

%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
% load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(index),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                         % Include the library of common functions in the search path
dataset_zero_thr = .05;                                                         % This is the threshold below which we set entries in the dataset to zero

b=clock;                                                                        % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*b))); 

learn_itr_max = 4000;                                                           % The maximum number of learning iterations
max_y_threshold = .05;                                                          % The maximum acceptable projection for a pattern on a weight vector
min_W = .005;                                                                   % The minimum absolute value of the weight vectors
fraction_of_convergence = .98;                                                  % Required fraction of patterns to converge
d_max = 15;                                                                     % The maximum allowed degree for pattern neurons
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


destination_folder = [db_file_in(1:i-1),'/Learn_Results/Zero_One/N_',num2str(N)];
if (~exist(destination_folder,'dir'))
    mkdir(destination_folder);
end
%--------------------------------------------------------------------------


%---------------------------Load the Training Dataset----------------------
load(db_file_in);
eval(['dataset_learn = ',db_name_in,';']);   
[dataset_size,pattern_length] = size(dataset_learn);

dataset_learn = dataset_learn .*(abs(dataset_learn)>dataset_zero_thr);
dataset_learn  = ones(dataset_size,pattern_length).*((dataset_learn>dataset_zero_thr)-(dataset_learn<-dataset_zero_thr));
% a = sqrt(sum(dataset_learn'.*dataset_learn'));
% dataset_learn = dataset_learn./(a'*ones(1,1024));
% dataset_learn = dataset_learn.*(abs(dataset_learn)>.03);
%--------------------------------------------------------------------------

%-------------Update the Number of Constraints in the Process--------------
fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'r'); 
if (fid > -1)
    current_count = max(fscanf(fid, '%d'),0);            
    fclose(fid);    
else    
    current_count = 0;            
end        

fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'r');
if (fid > -1)           
    W = fscanf(fid, '%f',[N,inf]);                
    W = W';            
    fclose(fid);                                    
    [m,~] = size(W);            
else    
    m = 0;    
end

const_to_learn = min(const_learn,N-m-current_count);

fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'w'); 
fprintf(fid,'%d',current_count+const_to_learn);
fclose(fid);
%--------------------------------------------------------------------------

%==========================================================================
    

%%
%=============================LEARNING PHASE===============================
for itr = 1: const_to_learn

    alpha1 = alpha0;
    theta1 = theta0;
    %-----------------------Read Current Learned Matrix--------------------
    fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'r');
    if (fid > -1)           
        W_tot = fscanf(fid, '%f',[N,inf]);                
        W_tot = W_tot';            
        fclose(fid);                                    
        [m,~] = size(W_tot);            
    else    
        m = 0;    
    end
    
    if (m >= 2)
       deg = sum(abs(sign(W_tot)));
    else
       deg = zeros(1,N);
    end                        
    %----------------------------------------------------------------------

    %------------------Randomly Initialize the Weight Matrix---------------
    W = random_vector(N,2*round(log(N)));
    W = W';
    W = rand(N,1);
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
        if (mod(learn_itr,round(log(1+m)*35))==0)
%             W = random_vector(N,2*round(log(N)));
%             W = W';
            W = rand(N,1);
            W = W/norm(W);
            theta_itr = theta_itr + 1;                
            theta = theta1/(2+theta_itr);
            alpha1 = alpha0/theta_itr;
        end

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
            x = dataset_learn(mu,:);
            %--------------------------------------------------------------
        
            %---------------------Update the Weight Vector-----------------
            y = x*W;                 % Calculate the projection of x on the vector W(j,:).
            
            if (norm(x)>0)
                alph = alph0/norm(x)^2;
            else
                alph = 0;
            end
            if (norm(W) > .001)                                
%                 W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,0.05);W.*(1-(tanh(100*W(j,:).^2)).^2)
%                  W = soft_threshold(W - alph*y*( x' - (y*W/(norm(W))^2)),theta0);%-beta0*W.*(1-(tanh(100*W.^2)).^2);
%                  W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
%                  W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta)-.05*W.*(1-(tanh(100*W.^2)).^2).*(1+tanh(.5*(deg-d_max)))';
                W = soft_threshold(W - alph*y*(x'-y*W)+ alph*(1-norm(W)^2)*W-.02*W.*(1-(tanh(100*W.^2)).^2).*((1+tanh(.05*(deg-d_max))')),theta); 
%                  
            end            
            
            W = soft_threshold(W,theta);                       % Set entries less than 0.0001 to zero.
            %--------------------------------------------------------------
        
            
            %------------------Check for Numerical Errors------------------
            if (sum(isnan(W))>0)    
                W = rand(N,1);
                W = W/norm(W);
                theta_itr = theta_itr + 1;                
                theta = 2*theta1/(2+theta_itr);
                display('Error Nan');        
%                 break;                                    
            end
            
            if (norm(W) < .001)                
                display('Error zero!'); 
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
        
              
        
        if (sum(abs(dataset_learn*W)<max_y_threshold)>fraction_of_convergence*dataset_size)
            W=soft_threshold(W,min_W);         
        end
        
        ww = W;
        ww(~ww) = inf;
        
        %--------------------------Display Progress------------------------
        if (mod(learn_itr,2)==0)
            theta_itr = theta_itr + 1;
            theta = 3*theta1/(3+theta_itr);            
            
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);
            display(['Converged count = ',num2str(sum((abs(dataset_learn*W)<max_y_threshold)))]);
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);
            display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
            display(['Maximum projection = ',num2str(max(abs(dataset_learn*W)))]);  
            display(['Minum of W = ',num2str(min(abs(ww)))]);  
            display('-----------------------------------');
            111;
            display(' ');             
        end           
        %------------------------------------------------------------------
         
            
        %----------------------Check Convergence---------------------------
        ww = W;
        ww(~ww) = inf;

        if ((sum((abs(dataset_learn*W))<max_y_threshold)>=fraction_of_convergence*dataset_size)&&(min(abs(ww))>min_W)&&(norm(W)>.1))
            break;
        end
        %------------------------------------------------------------------
        
        

    end  
    %======================================================================


    %==========================SAVE THE RESULTS============================    
    if ((sum((abs(dataset_learn*W))<max_y_threshold)>=fraction_of_convergence*dataset_size)&&(min(abs(ww))>min_W)&&(norm(W)>.1))                                
                        
        display(['-------Convergence!!!--------']);        
        display(['Converged count = ',num2str(sum((abs(dataset_learn*W)<max_y_threshold)))]);        
        display(['Cost in last iteration = ',num2str(cost_one_itr)]);
        display(['Norm-2 of W = ',num2str(norm(W))]);
        display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
        display(['Maximum projection = ',num2str(max(abs(dataset_learn*W)))]);  
        display(['Minum of W = ',num2str(min(abs(ww)))]);  
        display('-----------------------------------');
        display(' ');   
        
        fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'a');
        fprintf(fid, '%f \t',W');
        fprintf(fid, '\N');
        fclose(fid);
    

        fid = fopen([destination_folder,'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'a');
        fprintf(fid, '%d \t',learn_itr);
        fprintf(fid, '\N');
        fclose(fid);

%         fid = fopen(['/scratch/amir/Clustered_Neural/Learn_Results/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/Learn_cost_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'a');
%         fprintf(fid, '%f \t',log(cost_per_itr_av));
%         fprintf(fid, '\N');
%         fclose(fid);
    
        fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'r'); 
        if (fid > -1)
            current_count = fscanf(fid, '%d');
            fclose(fid);
        else
            current_count = 0;
            display('error!');
        end
    
        fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'.txt'], 'w'); 
        fprintf(fid,'%d',max(current_count-1,0));
        fclose(fid);
    end
    %======================================================================
    
end
111;

   