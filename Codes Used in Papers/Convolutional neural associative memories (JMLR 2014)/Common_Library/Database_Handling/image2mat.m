%==========================================================================
%**************************FUNCTION: image2mat****************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% source_depository: The address of the database of original images, stored as rows of a big matrix (each image is mapped from a 2D array to a one-dimensional vector)
% format: Specifies the format of the images which should be read
% max_size: It is the maximum number of files that should be processed by the each job submitted to the cluster 
% destinaion_file: Specifies the name of .mat file that holds the processed files up to this point
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


function image2mat(source_depository,format,max_size,destination_db_name,destinaion_file,option,height,width)


%===========================INITIALIZATION=================================

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
%--------------------------------------------------------------------------

%--------------------------Image Specifications----------------------------
% max_level = 16;                                                       % This is the maximum number of gray-levels we would like to have in the produced patches
listing = dir([source_depository,'/*.',format]);                        % This is the list of images in the depository with the specified folder
max_itr = floor(length(listing)/max_size);                              % This is the number of times the code should run to store all the images in the database
i_max = max_itr;                                                        % This is the number of times the code should run to store all the images in the database
%--------------------------------------------------------------------------

%---------------Create and Save The Empty Database File--------------------
eval([destination_db_name '= []']);
save(destinaion_file,destination_db_name);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------

%--------------------Initialize the SSH Connection-------------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');
%--------------------------------------------------------------------------

%==========================================================================


%=======================READ FILES AND STORE THEM==========================

for i = 0:i_max 
    offset = i*max_size;            
    command = ['cd /scratch/amir/Spatially_Coupled/Submit_Cluster;qsub -N "spatialy_image2mat_',num2str(i),...
        '" -v source=''"''"''',source_depository,'''"''"'',dest_file_name=''"''"''',destination_db_name,'''"''"'',format=''"''"''',format,'''"''"'',offset=',num2str(offset),',max_size=',num2str(max_size),...
        ',destination=''"''"''',destinaion_file,'''"''"'',option=',num2str(option),',height=',num2str(height),',width=',num2str(width),' image_to_mat.pbs'];
    
        
    %--------------------Submit Each Jobs to the Cluster-------------------
    [channel, result]  =  sshfrommatlabissue(channel,command);     
                                                               
    if (isequal(result, {''}))                             
        error('What the...! Unsubmitted image2mat job! :(');
    else                    
        display('Image2mat job submitted successfully! :D');
    end            
    %----------------------------------------------------------------------
    
    %----Verify Convergence of the Previous Job Before Submitting Next-----
    wait_flag = 1;
    while (wait_flag == 1)
        load(destinaion_file);
        eval(['[current_size,~] = size(',destination_db_name,');']);
        
        %-------------Determine the Desired Output File Size---------------
        if (i < i_max)
            desired_size = offset+max_size;
        else
            desired_size = length(listing);
        end
        %------------------------------------------------------------------
        
        if (current_size==desired_size)
            wait_flag = 0;
        else
            %-----------------Check for Possible Errors--------------------
            fid = fopen('/home1/amir/cluster/Spatially_Coupled/img2mat_error.txt','r');
                        
            if (fid>-1) 
                s = fscanf(fid,'%s');
                fclose(fid);
                if (strcmp(s,'error'))
                    error('Something has went wrong!');
                end
            end
            %--------------------------------------------------------------
            
            %-------------If No Errors, Wait Until Convergence-------------
            pause(200);
            %--------------------------------------------------------------
        end
    end
    %----------------------------------------------------------------------
end

%==========================================================================


        
