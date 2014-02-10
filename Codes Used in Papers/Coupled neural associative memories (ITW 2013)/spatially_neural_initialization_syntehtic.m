%==========================================================================
%****************CODE: spatially_neural_initialization*********************
%==========================================================================

%----------------------------CODE DESCRIPTION------------------------------
% This piece of code initializes the necessary parameters for the clustered
% neural associative memory and stores them on appropriate files so that
% other functions can have access to them. The advantage of storing them on
% hard disk is that we can run a parallel version of the code which makes
% it much more faster.

%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

% function spatially_neural_initialization(N_horiz,N_v)

%=============================INITIALIZATION===============================
if (~exist('initialization_done','var'))    % If already not initialized by the GUI...
    
    %-------------------------Network Parameters---------------------------    
    N_horiz = 64;                           % The 'width' of the 2D patterns
    N_vert = 64;                            % The 'height' of the 2D patterns    
    L_vert = 29;                            % The number of the neural 'planes'
    L_horiz = 29;                           % The number of clusters within each neural plane
    deg_horiz = 4;                          % The number of clusters each pattern neuron should be connected to
    deg_vert = 4;                           % The number of planes each pattern neuron should be connected to
            
    N_loc_horiz = 8;                        % The 'width' of the rectangular window used to divide the 2D patterns into patches
    N_loc_vert = 8;                         % The 'height' of the rectangular window used to divide the 2D patterns into patches
    
    N_const = N_loc_horiz*N_loc_vert/2;     % The maximum possible number of constraints which a network should learn
    %----------------------------------------------------------------------
               
    %-------------------------Simulation Parameters------------------------    
    learn_itr_max = 2000;                   % This is the number of times that the learning phase is repeated for the patterns in the training set       
    pattern_learn_number = 10000;           % This is the number of patterns used in the learning process.    
    dataset_synthetic = [];                 % The dataset containing the generated patterns
    dataset_address = ['/scratch/amir/Spatially_Coupled/Initialization_Files/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/syntehtic_db.mat'];
    %----------------------------------------------------------------------
    
    %--------------------------Neural Parameters---------------------------
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.       
    %----------------------------------------------------------------------       
    
    %----------------Generator Matrix Parameters---------------------------           
    deg_row_G_local = 8;                                        % This is the average number of non-zero elements in each row of the block corresponding to one neural 'plane'
    deg_column_G_local = 4;                                     % This is the average number of non-zero elements in each column of the block corresponding to one neural 'plane'
    delta_vert = 2;                                             % This is the amount of shift within each block to construct the generating block of each plane
    delta_N = N_loc_horiz/deg_horiz;                            % The amount of horizontal shift in each iteration
    K_loc = N_loc_horiz/2;                                      % This is the number of rows in each block for constructing the generator matrix for each 'plane' 
    delta_K = floor(.8*(N_loc_horiz-K_loc)/(deg_horiz-1));      % The amount of vertical shift in each iteration
    %----------------------------------------------------------------------    
    
    %----------------Check the Validty of Input Parameters-----------------
    if (abs(N_horiz - N_loc_horiz*(L_horiz-1+deg_horiz)/deg_horiz )>0)
        error('Invalid horizontal clustering!');
    elseif (abs(N_vert - N_loc_vert*(L_vert-1+deg_vert)/deg_vert )>0)
        error('Invalid vertical clustering!');
    end
    
    if (delta_K <=0)
        error('Invalid input parameters!');
    end
    %----------------------------------------------------------------------            
    
end

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                 % Include the library of common functions in the search path
%--------------------------------------------------------------------------

%==========================================================================

%%
%================CREATE THE GENERATOR MATRIX & PATTERNS====================

%---------------------Create the Generator Matrix--------------------------
[G_tot,K_tot] = generator_matrix_vert(N_horiz,N_vert,N_loc_horiz,N_loc_vert,K_loc,deg_column_G_local,deg_row_G_local,delta_N,delta_K,delta_vert);
%--------------------------------------------------------------------------

%------------------------Generate the Patterns List------------------------
KK=round(log(pattern_learn_number)/log(2));
mu_list = zeros(KK+1,pattern_learn_number);
for i = 1:pattern_learn_number
    randindex = randperm(K_tot);
    mu_list(1,i) = 1+floor((pattern_learn_number-1)*rand);
    mu_list(2:KK+1,i) = randindex(1:KK);         
end
%--------------------------------------------------------------------------

%------------Generate the Patterns and Store Them In a Database------------
for ii = 1:pattern_learn_number                            
    mu = 1+floor((pattern_learn_number-1)*rand);   
    temp = dec2bin(mu_list(1,mu),KK);      
    randindex = mu_list(2:KK+1,mu);    
    message = zeros(1,K_tot);                                                   % Generate the message from the index                                    
    for j = 1:KK                
        message(randindex(j)) = (z_max-z_min)*(temp(j) - 48)+z_min;                                         
    end        
    x = message*G_tot;
    dataset_synthetic = [dataset_synthetic;x];            
end
%--------------------------------------------------------------------------

%==========================================================================

%%
%============================STORE THE SETUP===============================
if (~exist(['/scratch/amir/Spatially_Coupled/Initialization_Files/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert)],'dir'))
    mkdir(['/scratch/amir/Spatially_Coupled/Initialization_Files'],['N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert)]);        % Create a specific folder for the current N and K
end

save(dataset_address,'dataset_synthetic');

clear dataset_synthetic

save(['/scratch/amir/Spatially_Coupled/Initialization_Files/N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/spatially_neural_parameters_synthetic_N_horiz_',num2str(N_horiz),'_N_vert_',num2str(N_vert),'_L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'.mat']);
%==========================================================================





