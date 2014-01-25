%==========================================================================
%*******************FUNCTION: image_database_chopper***********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% source_database_address: The address of the database of original images, stored as rows of a big matrix (each image is mapped from a 2D array to a one-dimensional vector)
% source_database_name: The name of the database which is stored in the specified address
% height: the height of the images in the database
% width: The width of the images in the database
% widnow_height: The height of the rectangualr moving window which is used to chop the original image into smaller parts
% widnow_width: The width of the rectangualr moving window which is used to chop the original image into smaller parts
% window_shift_horiz: The amount of the horizontal shift applied to the moving window in each iteration of the chopping process
% window_shift_vert: The amount of the vertical shift applied to the moving window in each iteration of the chopping process
% destination_depository: A string containing the address of the folder which will store the produced datasets. Example: '/home/hesam'
% options: A vector (c,p) of integers. If c = 1, chopping is done. If p = 1 or 2, projection will be done, based on two different methods to determine the dominant eigenvalues.
% execute_index: A matrix where a 0 at position (i,j) will tell the function to run the chopping job for the patch (i,j) of the image.
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% NONE
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function receives the specification of a database of image files (in
% .mat format). The database is arranged such that every image is converted 
% to a one-dimensional vector (instead of a 2D array) and then each vector is
% stored as one of the rows in the input database. Based on the parameters
% specified, the function starts by moving a window horizontally and 
% vertically and storing that particular part of the image for all the
% images in the database in a separate file. 
% Finally, if needed, the function can also projects the chopped part to a 
% smaller dimension by only keeping the dominant principal components.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function image_database_chopper(source_database_address,source_database_name,height,width,window_height,window_width,window_shift_horiz,window_shift_vert,destination_depository,options,success_status)


%===========================INITIALIZATION=================================

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
%--------------------------------------------------------------------------

%--------------------------Image Specifications----------------------------
% max_level = 16;                                                       % This is the maximum number of gray-levels we would like to have in the produced patches

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

if (~exist('success_status','var'))
    success_status = zeros(L_horiz,L_vert);
end
%--------------------------------------------------------------------------
% 
% %------------------Create the Status Marix If Necessary--------------------
% fid = fopen([destination_depository,'/success_status.mat']);
% if (fid == -1)
%     success_status = zeros(L_horiz,L_vert);
%     save([destination_depository,'/success_status.mat'],'success_status');
% else
%     load([destination_depository,'/success_status.mat']);
% end
% %--------------------------------------------------------------------------

%--------------------Initialize the SSH Connection-------------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');
%--------------------------------------------------------------------------

%==========================================================================


%==================MOVE THE WINDOW AND CHOP THE PATCHES====================
for index_horiz = 1:L_horiz
    for index_vert = 1:L_vert
        if (success_status(index_horiz,index_vert) == 0)
            %-----------------Determine Window's Specifications----------------
            window_start_coord =[1+(index_vert-1)*window_shift_vert,1+(index_horiz-1)*window_shift_horiz,];        
            window_end_coord = [window_height+(index_vert-1)*window_shift_vert,window_width+(index_horiz-1)*window_shift_horiz];      
            window_index = [index_horiz,index_vert];
            %------------------------------------------------------------------           
            
%             %------------------Submit Each Jobs to the Cluster-----------------
%             command = ['cd /scratch/amir/Spatially_Coupled/Submit_Cluster;qsub -N "spatialy_s_database_chopp_',num2str(index_horiz),'_',num2str(index_vert),...
%                 '" -v source=''"''"''',source_database_address,'''"''"'',source_name=''"''"''',source_database_name,'''"''"'',start_coord="''[',num2str(window_start_coord),']''",end_coord="''[',num2str(window_end_coord),...
%                 ']''",width=',num2str(width),',win_index="''[',num2str(window_index),']''",options="''[',num2str(options),']''",destination=''"''"''',destination_depository,'''"''"'' data_image_chop.pbs'];
%             [channel, result]  =  sshfrommatlabissue(channel,command);     
%                                                                  
%             if (isequal(result, {''}))                             
%                 error('Unsubmitted chopping job! :(');                                                             
%             else                    
%                 display('Chopping job submitted successfully! :D');                    
% %                 pause(2);
%             end            
%             %------------------------------------------------------------------
            image_database_extractor(source_database_address,source_database_name,window_start_coord,window_end_coord,width,window_index,options,destination_depository)
        end        
    end
    index_horiz
end
%==========================================================================

%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
%==========================================================================