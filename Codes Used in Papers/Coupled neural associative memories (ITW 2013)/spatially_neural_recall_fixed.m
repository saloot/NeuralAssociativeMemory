%==========================================================================
%************FUNCTION: spatially_neural_recall_fixed***********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N_horiz: The number of columns in the 2D patterns
% N_vert: The number of rows in the 2D patterns
% L_horiz: The number of clusters within a neural plane.
% L_vert: The number of neural planes
% no_simulated_instances: The number of patterns considered in the recall phase.
% p_err_init: Initial 'bit' error probability. It can be a vector.
% max_noise_amp: The maximum amplitude of noise during the recall process
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
% recall_algorithm_option: The recall algorithm identifier (0 for winner-take-all, 1 for the original bit flipping and 2 for the simplified bit flipping
% try_max: The maximum number of iterations for the global recall algorithm
% fixed_index: The 'length' of the rectangular window which is used to fix the boundaries during the recall phase. Put 0 for the uncostrained system
% db_name: The name of the variable used as the database in the corresponding .mat file
% db_file: The FULL name of the .mat file containing the database. Example: /home/amir/database.mat'.
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a spatially-coupled neural associative
% memory and performs the recall phase. More specifically, 
% the function first reads the weight matrix found in the learning phase. 
% Then the function performs the recall phase introduced in our Allerton 2012
% paper, i.e. perform the recall algorithm within each cluster one after
% another and then in each plane, one after another. If the recall process 
% was successful, the state of neurons is maintained, and reverted back to 
% their original version otherwise. This process is repeated try_max times, 
% after which a recall error is declared if the output of the algorithm is 
% not equal to the original pattern.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function spatially_neural_recall_fixed(N_horiz,N_vert,L_horiz,L_vert,no_simulated_instances,p_err_init,max_noise_amp,alpha0,beta0,theta0,gamma_BFO,gamma_BFS,recall_algorithm_option,try_max,fixed_index,db_name_in,db_file_in)

%=============================INITIALIZATION===============================

%-------------------------Simulation Parameters----------------------------
fixed_set_in = ones(N_horiz,N_vert);
if (fixed_index > 0)
    fixed_set_in(1:fixed_index,1:fixed_index) = 0;
    fixed_set_in(1:fixed_index,N_vert-fixed_index+1:N_vert) = 0;
    fixed_set_in(N_horiz-fixed_index+1:N_horiz,1:fixed_index) = 0;
    fixed_set_in(N_horiz-fixed_index+1:N_horiz,N_vert-fixed_index+1:N_vert) = 0;
elseif (fixed_index < 0)
    error ('Invalid fixed set size!');
end

% N_loc_horiz = 6;
% N_loc_vert = 6;
% window_shift_horiz = 2;
% window_shift_vert = 2;

n = N_loc_horiz*N_loc_vert;                                              % The total number of neurons in each cluster
%--------------------------------------------------------------------------

%----------------Load the Saved Initialized Parameters---------------------
load(['/scratch/amir/Spatially_Coupled/Initialization_Files/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_parameters_N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'.mat']);           
%--------------------------------------------------------------------------

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                 % Include the library of common functions in the search path
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

destination_folder = [db_file_in(1:i-1),'/Recall_Results/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert)];

if (~exist(destination_folder,'dir'))
    mkdir(destination_folder);
end
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 

error_count = 0;                                % The number of pattern errors captured during the simulations
bit_error_count = 0;                            % The number of bit errors captured during the simulations
y_min = -100;                                   % The minimum value a pattern neuron can assume
y_max = 100;                                    % The maximum value a pattern neuron can assume
%--------------------------------------------------------------------------

%-------------------------Load the Database--------------------------------
load(db_file_in);
eval(['dataset_recall = ',db_name_in,';']);                                
[pattern_learn_number,n_tot] = size(dataset_recall);
%--------------------------------------------------------------------------

%==========================================================================


%%
%===============================MAIN LOOP==================================
for net_simul_itr = 1:no_simulated_instances                % Simulate the error correction procedure for the given ensemble.                                
        
    %-----------------------------Generate Noise---------------------------    
    pp = rand(N_horiz,N_vert);                
    nois = (pp<p_err_init).*(((-1*ones(N_horiz,N_vert)).^randi(2,N_horiz,N_vert)));     % This is the noise added to the whole pattern of length N*L_in        
    nois = nois.*fixed_set_in;                                                          % Make sure there is no noise in the fixed boundaries
    %----------------------------------------------------------------------                            
                
    %-----------------------Pick a Random Pattern--------------------------    
    mu = 1+floor((pattern_learn_number-1)*rand);                                        % Pick a pattern index at random               
    pattern_tot = dataset_recall(mu,:);
    x_tot = pattern_tot + matrix2vector(nois);                            
    pattern_tot = vector2matrix(pattern_tot,N_vert);
    x_tot = vector2matrix(x_tot,N_vert);
    %----------------------------------------------------------------------        
        
    %------------------Adjust Simulation Parameters------------------------
    success_flag = 0;    
    try_itr = 0;
    %----------------------------------------------------------------------        
       
    %------------------------Perform the Recall Step-----------------------    
    while ((try_itr < try_max) && (success_flag == 0))
    
        %-----------------Update the Simulation Parameters-----------------
        success_flag = 1;
        unsuccess_count = 0;
        try_itr = try_itr+1;
        %------------------------------------------------------------------
        for l_v = 1:L_vert
            for l_h = 1:L_horiz                            
                n = N_loc_horiz*N_loc_vert;
                
                %---------------------Read the Weight Matrix-------------------                
                fid = fopen([source_folder_learn,'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_cluster_horiz_',num2str(l_h),'_cluster_vert_',num2str(l_v),'.txt'], 'r');                
                if (fid == -1)
                   continue
                end
                W = fscanf(fid, '%f',[n,inf]);                                        
                fclose(fid);
                W = W';      
                %----------------------------------------------------------                
                
                %---Determine the Coordination of the Rectangular Patch----
                if (l_h>1)
                    if (l_v>1)
                        ind_horiz_start = 1+(l_h-1)*window_shift_horiz;
                        ind_horiz_finish = ind_horiz_start+N_loc_horiz-1;
                        ind_vert_start = 1+(l_v-1)*window_shift_vert;
                        ind_vert_finish = ind_vert_start+N_loc_vert-1;
                    else
                        ind_horiz_start = 1+(l_h-1)*window_shift_horiz;
                        ind_horiz_finish = ind_horiz_start+N_loc_horiz-1;
                        ind_vert_start = 1;
                        ind_vert_finish = N_loc_vert;
                    end
                else
                    if (l_v>1)
                        ind_horiz_start = 1;
                        ind_horiz_finish = N_loc_horiz;
                        ind_vert_start = 1+(l_v-1)*window_shift_vert;
                        ind_vert_finish = ind_vert_start+N_loc_vert-1;
                    else
                        ind_horiz_start = 1;
                        ind_horiz_finish = N_loc_horiz;
                        ind_vert_start = 1;
                        ind_vert_finish = N_loc_vert;
                    end
                end
                %----------------------------------------------------------
                            
                %-----------------Generate the Sub-patterns----------------
                x_temp = matrix2vector(x_tot(ind_horiz_start:ind_horiz_finish,ind_vert_start:ind_vert_finish));
                fixed_set_temp = matrix2vector(fixed_set_in(ind_horiz_start:ind_horiz_finish,ind_vert_start:ind_vert_finish));                
                pattern_temp = matrix2vector(pattern_tot(ind_horiz_start:ind_horiz_finish,ind_vert_start:ind_vert_finish));                
                %----------------------------------------------------------
                                                               
                %-------------------Iterate Until Convergence--------------
                [x_out] = recall_step_spatial_fixed(W,x_temp,fixed_set_temp,n,gamma_BFO,gamma_BFS,max_noise_amp,y_min,y_max,recall_algorithm_option);                              
                %----------------------------------------------------------                            
                
                %--------Verify the Success of the Recall Algorithm--------
                if (norm(x_out-pattern_temp)<.05)                                   % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                                                                                                         
                    x_tot(ind_horiz_start:ind_horiz_finish,ind_vert_start:ind_vert_finish) = vector2matrix(x_out,N_loc_vert);                 
                else                         
                    unsuccess_count = unsuccess_count+1;
                    success_flag = 0;                                    
                end
                %----------------------------------------------------------
            
            end
            
        end
        
    end
    %----------------------------------------------------------------------
                                                          
                
    %--------------------Calculate Recall Error Count----------------------
    if (norm(x_tot-pattern_tot)>.05)                                            % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.        
        error_count = error_count+1;
        bit_error_count = bit_error_count+sum(sum(abs(sign(x_tot-pattern_tot))));        
    end
    %----------------------------------------------------------------------
    
    %--------------------------Display Progress----------------------------
    if (mod(net_simul_itr, 1) == 0)
        display(['-------------Iteration: ',num2str(net_simul_itr),'-------------']);
        display(['Error rate so far: ',num2str(error_count/net_simul_itr)]);        
        display(' ');
    end
    %----------------------------------------------------------------------
    
end    
%==========================================================================        
            

%%
%========================STORE THE RESULTS================================
    
%------------------Transform Error Count to Error Rate---------------------
PER = error_count/net_simul_itr;
BER = bit_error_count/net_simul_itr/n_tot;
%--------------------------------------------------------------------------
 
%------------------Store the Bit and Pattern Error Rates-------------------
if (fixed_index >0)
    if (recall_algorithm_option == 0)
        fid = fopen([destination_folder,'/spatially_neural_results_fixed_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
        fprintf(fid, 'e \t %f \t per \t %f \t ber \t %f \t',p_err_init,PER,BER);
        fprintf(fid,'\n');
        fclose(fid);
    elseif (recall_algorithm_option == 1)
        fid = fopen([destination_folder,'/spatially_neural_results_fixed_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
        fprintf(fid, 'e \t %f \t per \t %f \t ber \t %f \t',p_err_init,PER,BER);
        fprintf(fid,'\n');
        fclose(fid);
    elseif (recall_algorithm_option == 2)
        fid = fopen([destination_folder,'/spatially_neural_results_fixed_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
        fprintf(fid, 'e \t %f \t per \t %f \t ber \t %f \t',p_err_init,PER,BER);
        fprintf(fid,'\n');
        fclose(fid);
    else
        error('Unknown recall algorithm');
    end
else
    if (recall_algorithm_option == 0)
        fid = fopen([destination_folder,'/spatially_neural_results_WTA_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
        fprintf(fid, 'e \t %f \t per \t %f \t ber \t %f \t',p_err_init,PER,BER);
        fprintf(fid,'\n');
        fclose(fid);
    elseif (recall_algorithm_option == 1)
        fid = fopen([destination_folder,'/spatially_neural_results_BFO_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
        fprintf(fid, 'e \t %f \t per \t %f \t ber \t %f \t',p_err_init,PER,BER);
        fprintf(fid,'\n');
        fclose(fid);
    elseif (recall_algorithm_option == 2)
        fid = fopen([destination_folder,'/spatially_neural_results_BFS_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_gamma_BFO_',num2str(gamma_BFO),'_gamma_BFS_',num2str(gamma_BFS),'.txt'], 'a+');  
        fprintf(fid, 'e \t %f \t per \t %f \t ber \t %f \t',p_err_init,PER,BER);
        fprintf(fid,'\n');
        fclose(fid);
    else
        error('Unknown recall algorithm');
    end
end
%--------------------------------------------------------------------------

%==========================================================================

    
