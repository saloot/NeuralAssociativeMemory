%==========================================================================
%********************************READ ME***********************************
%==========================================================================

%--------------------------------Summary-----------------------------------
% This piece of code initializes the necessary parameters for the constraint
% enforcing neural network and stores them on appropriate files so that
% other functions can have access to them. The advantage of storing them on
% hard disk is that we can run a parallel version of the code which makes
% it much more faster.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

%%
%=============================INITIALIZATION===============================

if (~exist('initialization_done_by_main','var'))    % If already not initialized by the GUI..
    %-------------------------Network Parameters---------------------------
    N = 1600;                               % N is the number of neurons in network.
    K = 800;                                % K is the number of message bits.
    N_const = N-K;                          % N_cost represents the number of constraints.        
    deg_column_G = 15;
    deg_row_G = N*deg_column_G/K;    
    %----------------------------------------------------------------------

    %-------------------------Simulation Parameters------------------------
    learn_itr_max = 10000;                              % This is the number of times that the learning phase is repeated for the patterns in the training set       
    index_max = 50;                                     % This is the maximum number of random scenarios generated for simulation
    pattern_learn_number = min(100000,2^K);             % This is the number of patterns used in the learning process.
    KK = round(log(pattern_learn_number)/log(2));
    %----------------------------------------------------------------------

    
    %--------------------------Neural Parameters---------------------------
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.                             
    %----------------------------------------------------------------------

    %--------------------------Other Initializations-----------------------    
    addpath(genpath('../Common_Library'));                                  % Include the library of common functions in the search path

    a=clock;                                                                % Initialize the seed for random number generation with the clock value.
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*a))); 
    %----------------------------------------------------------------------
    
    
    %---Determine if the Dataset is Synthetic or Contains Real Patterns----
    real_dataset_flag = 0;                   % If 0, patterns will be generated synthetically. If 1, the user should specify the path to the dataset
        
    if real_dataset_flag
        db_file_name = ['CIFAR_10_train_gray_scale_class_4.mat'];
        db_file_folder = ['./Database/CIFAR_10']    % Make sure that there is no "/" or "\" at the end of the folder name
        db_name_in = 'CIFAR_10_Gray_DB_class_4';    % The name of the matrix that contains the patterns as its rows
    end
    %----------------------------------------------------------------------
    
end
initialization_done_by_master = 1;
if real_dataset_flag == 0
    mkdir(['../Initialization_Files'],['N_',num2str(N),'_K_',num2str(K)]);   % Create a specific folder for the current N and K
else
    mkdir(db_file_folder,'Preprocessed');   % Create a specific folder for the current N and K
end
%==========================================================================


%%
%======================IF THE DATASET IS SYNTHETIC=========================

%-------------------------Create the Generator Matrix----------------------
if real_dataset_flag == 0
    for index = 1:index_max

        %--------------Create the Pattern Generating Matrix----------------
        G = bipartite(N,K,deg_column_G,deg_row_G);
        %------------------------------------------------------------------

        %---------------Create the Training List Indices-------------------
        mu_list = 1+floor( (pattern_learn_number-1)*rand(1,pattern_learn_number));     
        %------------------------------------------------------------------
    
        %---------------------Save Training Sets---------------------------
        save(['../Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(index),'.mat'],'G','mu_list');
        %------------------------------------------------------------------
    end
    
    %--------------------Save Initialized Variables------------------------
    save(['../Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_parameters_v1_N_',num2str(N),'_K_',num2str(K),'.mat']);                    % Store all the variables in the appropriate file. 
    %----------------------------------------------------------------------

%--------------------------------------------------------------------------
%==========================================================================

%%
%==================IF THE DATASET IS REAL, E.G., IMAGES====================
else    
    %---------------Load the Real-Valued Training Dataset------------------
    load([db_file_folder,'/',db_file_name]);
    eval(['dataset_learn = ',db_name_in,';']);       
    %----------------------------------------------------------------------
    
    %-------------------Quantize the Dataset-----------------------------
    [C,N] = size(dataset_learn);            
    max_dataset = max(max(dataset_learn));
    dataset_learn = round(Q*dataset_learn/max_dataset);    
    %----------------------------------------------------------------------
        

    %--------------------Transform the Dataset to Bits---------------------
    q = (1+floor(log(Q)/log(2)));
    dataset_learn_binary = zeros(C,N*q);
    for i = 1:C
         temp = de2bi(dataset_learn(i,:),q);
         
         %............Set the Least Significant Values to Zero.............         
         temp(:,1) = 0;
         %temp(:,2) = 0;
         %.................................................................
         temp = reshape(temp,[1,N*q]);
         dataset_learn_binary(i,:) = temp;
    end
    %N = N*q;
    dataset_learn = dataset_learn_binary;
    %----------------------------------------------------------------------


    %----------------Store the Preprocessed Dataset------------------------
    save([db_file_folder,'/Preprocessed/','Preprocessed_',db_file_name],'dataset_learn');
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------
%==========================================================================



