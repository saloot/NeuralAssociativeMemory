%==========================================================================
%********************FUNCTION: clustered_neural_learn**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N_horiz_in: The number of columns in the 2D patterns
% N_vert_in: The number of rows in the 2D patterns
% L_horiz_in: The number of clusters within a neural plane.
% L_vert_in: The number of neural planes
% const_learn: Number of constraints which must be learned during the learning phase
% cluster_index_horiz_in_in: The horizontal index of the corresponding cluster
% cluster_index_vert_in_in: The vertical index of the corresponding cluster
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% db_name_in: The name of the variable used as the database in the corresponding .mat file
% db_file_in: The FULL name of the .mat file containing the database. Example: /home/amir/database.mat'.
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

function spatially_neural_learn(N_horiz,N_vert,L_horiz,L_vert,const_learn,cluster_index_horiz_in_in,cluster_index_vert_in_in,alpha0,beta0,theta0,db_name_in,db_file_in)

%%
%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Spatially_Coupled/Initialization_Files/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_parameters_N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'.mat']);           
learn_itr_max = 1000;
%--------------------------------------------------------------------------

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                 % Include the library of common functions in the search path
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
max_y_threshold = .01;                                                  % The threshold below which we count the projections as zero
n = N_loc_horiz*N_loc_vert;                                             % The total number of neurons in each cluster
if (~exist('N_const','var'))
    N_const =N_horiz;                                                   % The maximum possible number of constraints which a network should learn
end

a=clock;                                                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------


%-------------------------Load the Database--------------------------------
load(db_file_in);
eval(['dataset_learn=',db_name_in,';']);
[pattern_learn_number,~] = size(dataset_learn);
%--------------------------------------------------------------------------


%--------------------Create the Sub-folder If Necessary--------------------
slash_flag = 0;
for i = length(db_file_in):-1:1
    if (strcmp(db_file_in(i),'/'))
        if (slash_flag == 2)
            break;
        else
            slash_flag = slash_flag+1;
        end
    end
end

destination_folder = [db_file_in(1:i-1),'/Learn_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert)];

if (~exist(destination_folder,'dir'))
    mkdir(destination_folder);
end
%--------------------------------------------------------------------------

%==========================================================================


%%
%=========CALCULATE THE ACTUAL NUMBER OF CONSTRAINTS TO BE LEARNED=========

%------------Determine the Number of Already Learned Constraints-----------
fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'r');        

if (fid > -1)    
    W = fscanf(fid, '%f',[n,inf]);            
    W = W';            
    fclose(fid);                                                        
    [m,~] = size(W);            
else    
    m = 0;    
end
%--------------------------------------------------------------------------

%--------------Take Into Account the Constraints Being Learned-------------
fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'r'); 

if (fid > -1)                
    current_count = max(fscanf(fid, '%d'),0);               
    fclose(fid);    
else    
    current_count = 0;
end
%--------------------------------------------------------------------------

%---------Determine the Actual Number of Constraints to be Learned---------
no_const_learn = min(const_learn,N_const-m-current_count);
%--------------------------------------------------------------------------        

%---------------Update the Number of Constraints Being Learned-------------
fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'w');
fprintf(fid,'%d',current_count+no_const_learn);
fclose(fid);
%--------------------------------------------------------------------------        

%==========================================================================

%%
%=============================LEARNING PHASE===============================
for itr = 1:no_const_learn                 % Do untill all MC's are learned                             
   
    %------------------Randomly Initialize the Weight Matrix---------------
    W = random_vector(n,min(6*round(log(n)),n));
    W = W';
%     W = rand(n,1);
    W = W/norm(W);

    %----------------------------------------------------------------------
    
    %---------------------Adjust Simulation Parameters---------------------
    theta = theta0;   
    theta_itr = 1;
    cost_per_itr_av = [];                    
    %----------------------------------------------------------------------

    
    %--------------------------Main Learning Loop--------------------------           
    for learn_itr = 1:learn_itr_max                         % Repeat the steps below several times in order to enhance the quality of the learning phase        

        if (mod(learn_itr,100)==0)            
            W = rand(n,1);
            W = W/norm(W);            
            theta_itr = theta_itr + 1;                
            theta = theta0/(2+theta_itr);        
        end
        
        %---------------Adjust In-loop Simulation Parameters---------------
        cost_short_term =0;        
        cost_one_itr = 0;        
        max_y = 0;
        converged_count = 0;                        
        alph0 = max(50*alpha0/(50+alph_itr),0.0005);                
        %------------------------------------------------------------------
        
        for ii = 1:pattern_learn_number
            
            %------------------Pick a Pattern at Random--------------------
            mu = 1+floor((pattern_learn_number-1)*rand);                                            
            x = dataset_learn(mu,:);
            IMG = vector2matrix(x,N_vert);

            if (cluster_index_horiz_in>1)
                if (cluster_index_vert_in>1)
                    ind_horiz_start = 1+(cluster_index_horiz_in-1)*N_loc_horiz/deg_horiz;
                    ind_horiz_finish = ind_horiz_start+N_loc_horiz-1;
                    ind_vert_start = 1+(cluster_index_vert_in-1)*N_loc_vert/deg_vert;
                    ind_vert_finish = ind_vert_start+N_loc_vert-1;
                else
                    ind_horiz_start = 1+(cluster_index_horiz_in-1)*N_loc_horiz/deg_horiz;
                    ind_horiz_finish = ind_horiz_start+N_loc_horiz-1;
                    ind_vert_start = 1;
                    ind_vert_finish = N_loc_vert;
                end
            else
                if (cluster_index_vert_in>1)
                    ind_horiz_start = 1;
                    ind_horiz_finish = N_loc_horiz;
                    ind_vert_start = 1+(cluster_index_vert_in-1)*N_loc_vert/deg_vert;
                    ind_vert_finish = ind_vert_start+N_loc_vert-1;
                else
                    ind_horiz_start = 1;
                    ind_horiz_finish = N_loc_horiz;
                    ind_vert_start = 1;
                    ind_vert_finish = N_loc_vert;
                end
            end
            
            A = IMG(ind_horiz_start:ind_horiz_finish,ind_vert_start:ind_vert_finish);
            
            x = matrix2vector(A);                                               
            %--------------------------------------------------------------
            
            %-------------------Adjust the Learning Rate-------------------
            if (norm(x)>0.0001)
                alph = alph0/norm(x)^2;            
            else
                alph = 0;
            end            
            %--------------------------------------------------------------
            
            %---------------------Update the Weight Vector-----------------
            y = x*W;                                     % Calculate the projection of x on the vector W(j,:).
            
            if (norm(W) > .01)                
                % W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);%W.*(1-(tanh(100*W(j,:).^2)).^2)
                W = soft_threshold(W - alph*y*( x' - (y*W/(norm(W))^2)),theta);%-beta0*W.*(1-(tanh(100*W.^2)).^2);                
                % W = W - alph*y*( x' - (y*W/(norm(W))^2))-beta0*soft_threshold_inv(W,theta);
                %  W = soft_threshold(W - alph*y*(x'*(2*norm(W)^2-1) - y*W),theta);
                % W = soft_threshold(W - alph*y*(x'-y*W)+ alph*(1-norm(W)^2)*W,theta);                           
            end                                    
            % W = soft_threshold(W,.001);
            %--------------------------------------------------------------
                       
            %------------------Check for Numerical Errors------------------           
            if (sum(isnan(W))>0)    
                display('Error Nan'); 
                break;                                    
            end
            
            if (norm(W) < .001)                
                display('Error zero!');                
                break;                                    
            end           
            %--------------------------------------------------------------

            %-----------------Update the Simulation Costs------------------
            cost_one_itr = cost_one_itr + abs(x*W);
            cost_short_term = cost_short_term+norm(x*W);

            if (abs(x*W)> max_y)
                max_y = abs(x*W);
            end
            
            if (norm(x*W) < max_y_threshold)
                converged_count = converged_count + 1;
            end
            %--------------------------------------------------------------                                                                                      
        end
                       
        
        %-----------------------Display Progress-----------------------       
        if (mod(learn_itr,5)==0)
            display(['-------Learn iteration: ',num2str(learn_itr),'--------']);
            display(['Converged count = ',num2str(converged_count)]);
            display(['Cost in last iteration = ',num2str(cost_one_itr)]);
            display(['Norm-2 of W = ',num2str(norm(W))]);
            display(['Norm-0 of W = ',num2str(sum(abs(sign(W))))]);
            display(['Maximum projection = ',num2str(max(abs(max_y)))]);              
            display('-----------------------------------');
            display(' ');            
        end
        %------------------------------------------------------------------                
        
        
        %--------------------Update the Simulation Costs-------------------
        cost_per_itr_av = [cost_per_itr_av,cost_short_term];
        cost_short_term = 0;
        %------------------------------------------------------------------                
        
        %---------------------------Update theta---------------------------
        if (mod(learn_itr,100)==0)
            theta_itr = theta_itr + 1;
            theta =theta0/(theta_itr);
        end                            
        %------------------------------------------------------------------
        
        %----------------------Check Convergence---------------------------
        if (max_y<max_y_threshold)
            if (surely_flag == 1)                    
                break;                
            else                
                surely_flag = 1;                
            end            
        end      
        %------------------------------------------------------------------ 
    

    end  
    %======================================================================
    

    
    %==========================SAVE THE RESULTS============================    
    if (max_y<max_y_threshold)
        fid = fopen([destination_folder,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'a');
        fprintf(fid, '%f \t',W');
        fprintf(fid, '\n');
        fclose(fid);
   
        fid = fopen([destination_folder,'/Learn_itr_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'a');
        fprintf(fid, '%d \t',learn_itr);
        fprintf(fid, '\n');
        fclose(fid);

        fid = fopen([destination_folder,'/Learn_cost_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'a');
        fprintf(fid, '%f \t',log(cost_per_itr_av));
        fprintf(fid, '\n');
        fclose(fid);
    
        fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'r'); 
        if (fid > -1)
            current_count = fscanf(fid, '%d');
            fclose(fid);
        else
            current_count = 0;
            display('error!');
        end
    
        fid = fopen([destination_folder,'/learn_count_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(cluster_index_horiz_in_in),'_cluster_vert_',num2str(cluster_index_vert_in_in),'.txt'], 'w'); 
        fprintf(fid,'%d',max(current_count-1,0));
        fclose(fid);        
    end
    %======================================================================
    
end   