%==========================================================================
%******************FUNCTION: image_database_extractor**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% source_database_address: The address of the folder containing the database of original images, stored as rows of a big matrix (each image is mapped from a 2D array to a one-dimensional vector)
% source_database_name: The name of the database which is stored in the specified address
% window_start_coord: A vector (x,y) with the coordinates of the top left corner of the window used to specifiy a portion of the image to be chopped
% window_end_coord: A vector (x,y) with the coordinates of the bottom right corner of the window used to specifiy a portion of the image to be chopped
% width: The width of the images in the database
% window_index: A vector (a,b) containing the index of the considered image  patch. The first element shows which plane (set of rows in the image) the patch belongs to and the second element represent the index of the patch within each plane.
% destination_depository: A string containing the address of the folder which will store the produced datasets. Example: '/home/hesam'
% options: A vector (c,p) of integers. If c = 1, chopping is done. If p = 1 or 2, projection will be done, based on two different methods to determine the dominant eigenvalues.
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% NONE
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function receives specifications of a database of image files in .mat format. The
% database is arranged such that every image is converted to a
% one-dimensional vector (instead of a 2D array) and then each vector is
% stored as one of the rows in the input database. Based on the parameters
% specified, the function take out the portion of the image specified by a
% window of given size and store it into a new file, which makes it easier
% for further processing.
% Finally, if needed, the function projects the chopped part to a smaller
% dimension by only keeping the dominant principal components.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function image_database_extractor(source_database_address,source_database_name,window_start_coord,window_end_coord,width,window_index,options,destination_depository)

%=============================INTIALIZATIONS===============================
l_h = window_index(1);
l_v = window_index(2);
fid = fopen(source_database_address, 'r');                                   
input_dataset = fscanf(fid,'%f',[width^2,inf]);
fclose(fid);
input_dataset = input_dataset';
[data_set_size,~] = size(input_dataset);
% load(source_database_address);
% eval(['[data_set_size,~] = size(',source_database_name,');']);
%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
%--------------------------------------------------------------------------

%==========================================================================

%=======================IF CHPPOING SHOULD BE DONE=========================
if (options(1) == 1)                                         % If the database should be chopped into smaller sized datasets
    
    dataset_chopped = [];

    for i = 1:data_set_size
    
%         eval(['IMG = vector2matrix(',source_database_name,'(i,:),width);']);            
        IMG =  vector2matrix(input_dataset(i,:),width);
        
        A = IMG(window_start_coord(1): window_end_coord(1),window_start_coord(2): window_end_coord(2));                
        dataset_chopped = [dataset_chopped;matrix2vector(A)];                
    end
    dataset_chopped = double(dataset_chopped);

    %----------------------------Save the Result---------------------------
    save([destination_depository,'/dataset_chopped_',num2str(l_h),'_',num2str(l_v),'.mat'],'dataset_chopped');
    %----------------------------------------------------------------------
    
end
%==========================================================================

%=====================IF PROJECTION SHOULD BE DONE=========================
if (options(2) > 0)                                        % If the limit on the latent varianbles are on INIDIVIDUAL elements...            
    
    %-----------------------Load the Chopped Dataset-----------------------
    load([destination_depository,'/dataset_chopped_',num2str(l_h),'_',num2str(l_v),'.mat']);
    %----------------------------Perform PCA-------------------------------
    
    %----------------------------Perform PCA-------------------------------
    [COEFF,SCORE,latent] = princomp(dataset_chopped);             
    %----------------------------------------------------------------------
    
    %--------Determine the Number of Dominant Princiapl Components---------
    for i = 1:length(latent)
        
        if (options(2) == 1)                            % If the limit on the latent varianbles are on the ABSOLUTE VALUE of eigenvalues...
            if (latent(i)<.005)        
                no_of_PC = i;
                break;            
            else
                no_of_PC = length(latent);
                display('Dataset is full rank!');
            end
        elseif (options(2) == 2)                        % If the limit on the latent varianbles are on the maximum ENERGY of the eigenvalues...
            if (sum(latent(1:i))/sum(latent)>.99)
                no_of_PC = i;
                break;
            end
        end
    end
    %----------------------------------------------------------------------
    
    %----------------------Perform the Projection--------------------------
    dataset_projected = SCORE(:,1:no_of_PC)*COEFF(:,1:no_of_PC)';
    features = SCORE;
    features_sparse= [];
    
    no_of_features = length(latent);
    for i = 1:data_set_size
        [s,ind] = sort(abs(SCORE(i,:)),'descend');
        for j = 1:length(s)
            if (sum(s(1:j).^2)/norm(s)^2>.95)
                no_of_features = j;
                break;
            end
        end
        
        no_of_features = max(10,no_of_features);
        feature_temp = zeros(1,length(latent));
        for i = 1:no_of_features
            feature_temp(ind(i)) = s(i);
        end
%         feature_temp = SCORE(i,:).*(abs(SCORE(i,:))/max(abs(SCORE(i,:)))>.25);        
        features_sparse = [features_sparse;feature_temp];
    end
    %     dataset_projected = dataset_projected + ones(data_set_size,1)*mean(dataset); 
    %     dataset_projected = round(dataset_projected);
    %----------------------------------------------------------------------
    
    %----------------------Singular Value Decomposition--------------------
    [U,S,V] = svd(dataset_chopped);
    dataset_svd = U*S;
    %----------------------------------------------------------------------
    
    
    %----------------------------Save the Result---------------------------
    save([destination_depository,'/dataset_projected_',num2str(l_h),'_',num2str(l_v),'.mat'],'dataset_projected');    
    save([destination_depository,'/dataset_features_',num2str(l_h),'_',num2str(l_v),'.mat'],'features');    
    save([destination_depository,'/dataset_features_sparse_',num2str(l_h),'_',num2str(l_v),'.mat'],'features_sparse');    
    save([destination_depository,'/dataset_svd_',num2str(l_h),'_',num2str(l_v),'.mat'],'dataset_svd');    
    %----------------------------------------------------------------------    
end
%========================================================================== 
        

% load([destination_depository,'/success_status.mat']);
% success_status(window_index(1),window_index(2)) = 1;
% save([destination_depository,'/success_status.mat'],'success_status');
    
 
