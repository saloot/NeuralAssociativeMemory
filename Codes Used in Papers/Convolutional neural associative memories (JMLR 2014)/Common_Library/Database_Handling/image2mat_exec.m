%==========================================================================
%**************************FUNCTION: image2mat****************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% source_depository: The address of the database of original images, stored as rows of a big matrix (each image is mapped from a 2D array to a one-dimensional vector)
% format: Specifies the format of the images which should be read
% offset: It is the relative place of the file which the function should start processing. It is measured from the firt file in the folder
% max_size: It is the maximum number of files that should be processed by the code
% destination_file: A string containing the address of the file which will store the produced datasets. Example: '/home/hesam/test_file.mat'
% destination_db_name: Specifies the name of the MATLAB variable within the .mat file that holds the processed files up to this point
% option: If option = 1, the function also rescales the images and map them to gray scale levels.
% height: The desired height of the rescaled images
% width: The desired widthof the rescaled images
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% NONE
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the address of an image depository and stores the set
% of images within the depository as rows of a .mat file.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function image2mat_exec(source_depository,format,offset,max_size,destinaion_file,destination_db_name,option,height,width)


%===========================INITIALIZATION=================================

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
%--------------------------------------------------------------------------

%--------------------------Image Specifications----------------------------
% max_level = 16;                                                       % This is the maximum number of gray-levels we would like to have in the produced patches
listing = dir([source_depository,'/*.',format]);                        % This is the ist of images in the depository with the specified folder

if (isempty(listing))
    fid = fopen('/home1/amir/cluster/Spatially_Coupled/img2mat_error.txt','w');
    fprintf(fid,'error');
    fclose(fid);
    error('Depository is empty!');
end

load(destinaion_file);
%--------------------------------------------------------------------------

%==========================================================================



%=======================READ FILES AND STORE THEM==========================
for i = offset+1:min(offset+max_size,length(listing))
    file_name = listing(i).name;
    file_name_full = [source_depository,'/',file_name];
    IMG = imread(file_name_full,format);                        
    
    %----------If Resizing/Mapping to Gray Scale is Needed as Well---------
    if (option == 1)
        [~,~,c] = size(IMG);            
        if (c>1)                
            gray_fig = rgb2gray(IMG);            
        else            
            gray_fig = IMG;
        end
        
        %--------------------------Resize the Image------------------------
        [m,n] = size(gray_fig);
        scale = min(height/m,width/n);

        IMG = imresize(gray_fig,scale);
        %------------------------------------------------------------------
        
        %--------------------Save the Final Image--------------------------
        mkdir(source_depository,'Processed');
        imwrite(IMG,[source_depository,'/Processed/',file_name,'.gif'],'gif','ScreenSize',[width,height],'BackgroundColor',255);
        %------------------------------------------------------------------
        
        %-------------------Resize the MATLAB Files Too--------------------
        TEMP_IMG = zeros(height,width);
        [m,n]=size(IMG);
        TEMP_IMG(1:m,1:n) = IMG;
        IMG = TEMP_IMG;
        %------------------------------------------------------------------
        
    end
    %----------------------------------------------------------------------
    
    eval([destination_db_name '= [',destination_db_name,';matrix2vector(IMG)];']);
            
end
save(destinaion_file,destination_db_name);        
fid = fopen('/home1/amir/cluster/Spatially_Coupled/img2mat_error.txt','w');
fprintf(fid,' ');
fclose(fid);
%==========================================================================

