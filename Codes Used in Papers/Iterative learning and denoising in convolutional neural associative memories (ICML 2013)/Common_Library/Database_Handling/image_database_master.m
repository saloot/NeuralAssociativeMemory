%==========================================================================
%*********************CODE: image_database_master**************************
%==========================================================================


%---------------------------CODE DESCRIPTION-------------------------------
% This piece of code gets the specification of a repoistory of images and
% prepares an image database to be used later in neural applications.
% In particular, the code performs the following steps (if had not done
% before):
% 1) It transforms the images from color versions to gray-scale images.
% Furthermore, the images are rescaled if necessary.
% 2) Two dimesnional images are transformed to vectors and stored as rows 
% of a matrix. The matrix is then saved in a .mat file.
% 3) The large matrix is then broken into smaller pieces, corresponding to
% various patches of the images in the database. These smaller pathaces
% (and their projected versions) are stored separately. 
% 4) The code verifies if the submitted jobs to cluster have converged
% successfully.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================
%============================INITIALIZATION================================

%------------------------Repository Specifications-------------------------
% source_depository = ['/home/amir/Hesam/Research/Biocoding Project/Databases/STL-10'];       % The address of the database of original images, stored as rows of a big matrix (each image is mapped from a 2D array to a one-dimensional vector)
source_depository = '/scratch/amir/Databases/STL_10/Gabor/';
db_name = 'STL_10_Gray_DB';                                                    % Specifies the name of the MATLAB variable within the .mat file that holds the processed files up to this point
database_folder = [source_depository,'/Processed'];                                 % Specifies the folder containing the .mat file that holds the processed files up to this point
% database_address = [database_folder,'/',db_name,'.mat'];                                   % Specifies the name of .mat file that holds the processed files up to this point
database_address ='/scratch/amir/Databases/STL_10/Gabor/Filtered_image.txt';
file_name_base = ['/scratch/amir/Databases/STL_10/Gabor/Projected/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert),'/dataset_chopped_'];
%--------------------------------------------------------------------------

%-------------------------Database Specifications--------------------------
format = 'jpg';                                                                     % Specifies the format of the images which should be read
height = 96;                                                                        % The desired height of the rescaled images
width = 96;                                                                         % The desired wdith of the rescaled images
window_width = 24;                                                                   % The height of the rectangualr moving window which is used to chop the original image into smaller parts
window_height = 24;                                                                  % The width of the rectangualr moving window which is used to chop the original image into smaller parts
window_shift_horiz = 24;                                                             % The amount of the horizontal shift applied to the moving window in each iteration of the chopping process
window_shift_vert = 24;                                                              % The amount of the vertical shift applied to the moving window in each iteration of the chopping process
%--------------------------------------------------------------------------

%--------------------------Algorithm Parameters----------------------------
max_size = 100;                                                                     % It is the maximum number of files that should be processed by the each job submitted to the cluster 
option_img2mat = 1;                                                                 % If option_img2mat = 1, the function also rescales the images and map them to gray scale levels.
option_chop = [1 2];                                                                % A vector (c,p) of integers. If c = 1, chopping is done. If p = 1 or 2, projection will be done, based on two different methods to determine the dominant eigenvalues.
if (width-window_width>0)
    L_horiz = 1+ (width-window_width)/window_shift_horiz;                   % This is the maximum number of horizontal shifts
elseif (width-window_width == 0)
    L_horiz = 1;
else
    error('Invalid input width or window_width');
end

if (height-window_height>0)
    L_vert = 1+ (height-window_height)/window_shift_vert;                   % This is the maximum number of vertical shifts
elseif (height-window_height == 0)
    L_vert = 1;
else 
    error('Invalid input height or window_height');
end

destination_depository = [source_depository,'/Non_Overlap/Projected/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert)];                          % A string containing the address of the folder which will store the produced datasets. Example: '/home/hesam'
%--------------------------------------------------------------------------

%----------------------Create NEcessary Folders----------------------------
if (exist(database_folder,'dir') ~=7)
    mkdir(database_folder);
end
if (exist(destination_depository,'dir') ~=7)
    mkdir(destination_depository);
end
%--------------------------------------------------------------------------

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                     % Include the library of common functions in the search path
%--------------------------------------------------------------------------


% %--------------Check The Number of Files in The Reporistory----------------
% listing = dir([source_depository,'/*.',format]);                            % This is the list of images in the depository with the specified folder
% if (isempty(listing))
%     error('No files with these specifications in the repository');
% end
% %--------------------------------------------------------------------------
% 
% %==========================================================================
% 
% 
% %==============STORE THE IMAGES AS ROWS OF A MATLAB FILE===================
% 
% %---------------Check if Storing Has Already Been Done---------------------
% fid = fopen(database_address);
% if (fid > -1)
%     load(database_address);
%     eval(['[current_size,~] = size(',db_name,');']);
% else
%     current_size = 0;
% end
% %--------------------------------------------------------------------------
% 
% if (current_size ~= length(listing))
%     image2mat(source_depository,format,max_size,db_name,database_address,option_img2mat,height,width);
% end
% %==========================================================================


%================CHOP THE IMAGES INTO SMALLER PORTIONS=====================
success_index = image_db_verify(file_name_base,L_horiz,L_vert);
image_database_chopper(database_address,db_name,height,width,window_height,window_width,window_shift_horiz,window_shift_vert,destination_depository,option_chop,success_index)
%==========================================================================


%===============VERIFY THE SUCCESS OF THE CHOPPING PROCESS=================
success_index = image_db_verify(file_name_base,L_horiz,L_vert);
111
% pause
%-----------If Some Parts Have not Converged, Re-execute Them--------------
% if (sum(sum(success_index))<L_horiz*L_vert)
%     display('Some jobs have not been successful :( Resubmitting!');
%     image_database_chopper(database_address,db_name,height,width,window_height,window_width,window_shift_horiz,window_shift_vert,destination_depository,option_chop,success_index)
% end
%--------------------------------------------------------------------------

%==========================================================================