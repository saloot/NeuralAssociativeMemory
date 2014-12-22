%==========================================================================
%**************************FUNCTION: image2mat****************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% source_depository: The address of the database of original images, stored as rows of a big matrix (each image is mapped from a 2D array to a one-dimensional vector)
% max_size: It is the maximum number of files that should be processed by the each job submitted to the cluster 
% destinaion_file: Specifies the name of .mat file that holds the processed files up to this point
% destination_db_name: Specifies the name of the MATLAB variable within the .mat file that holds the processed files up to this point
% DS_rate: If DS_rate = 1, the function also rescales the images and map them to gray scale levels.
% height: The desired height of the rescaled images
% duration: The desired widthof the rescaled images
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


function wav2mat(source_depository,max_size,destination_db_name,destinaion_depository,DS_rate,duration,freq_threshold,N_out)


%===========================INITIALIZATION=================================

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
%--------------------------------------------------------------------------

%--------------------------Image Specifications----------------------------
% max_level = 16;                                                       % This is the maximum number of gray-levels we would like to have in the produced patches
listing = dir([source_depository,'/*.wav']);                        % This is the list of images in the depository with the specified folder
max_itr = floor(length(listing)/max_size);                              % This is the number of times the code should run to store all the images in the database
i_max = max_itr;                                                        % This is the number of times the code should run to store all the images in the database
if (i_max == 0)
    error('Folder is empty!');
end
mkdir(destinaion_depository,'DS');
mkdir(destinaion_depository,'FFT');
mkdir([destinaion_depository,'/DS'],'FFT');

%--------------------------------------------------------------------------

%---------------Create and Save The Empty Database File--------------------
eval([destination_db_name '= []']);

destinaion_file = [destinaion_depository,'/',destination_db_name,'.mat'];
save(destinaion_file,destination_db_name);
destinaion_file = [destinaion_depository,'/DS/',destination_db_name,'.mat'];
save(destinaion_file,destination_db_name);
destinaion_file = [destinaion_depository,'/FFT/',destination_db_name,'.mat'];
save(destinaion_file,destination_db_name);
destinaion_file = [destinaion_depository,'/DS/FFT/',destination_db_name,'.mat'];
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
    
%     command = ['cd /scratch/amir/Clustered_Neural/Submit_Cluster;qsub -N "wav2mat_cluster_',num2str(i),...
%         '" -v source=''"''"''',source_depository,'''"''"'',dest_file_name=''"''"''',destination_db_name,'''"''"'',offset=',num2str(offset),',max_size=',num2str(max_size),...
%         ',destination=''"''"''',destinaion_depository,'''"''"'',DS_rate=',num2str(DS_rate),',Nout=',num2str(N_out),',duration=',num2str(duration),',freq_threshold=',num2str(freq_threshold),' db_wave_to_mat.pbs'];
%     
%         
%     %--------------------Submit Each Jobs to the Cluster-------------------
%     [channel, result]  =  sshfrommatlabissue(channel,command);     
%                                                                
%     if (isequal(result, {''}))                             
%         error('What the...! Unsubmitted wav2mat job! :(');
%     else                    
%         display('Image2mat job submitted successfully! :D');
%     end            
%     %----------------------------------------------------------------------
%     
%     %----Verify Convergence of the Previous Job Before Submitting Next-----
%     wait_flag = 1;
%     while (wait_flag == 1)
%         load(destinaion_file);
%         eval(['[current_size,~] = size(',destination_db_name,');']);
%         
%         %-------------Determine the Desired Output File Size---------------
%         if (i < i_max)
%             desired_size = offset+max_size;
%         else
%             desired_size = length(listing);
%         end
%         %------------------------------------------------------------------
%         
%         if (current_size==desired_size)
%             wait_flag = 0;
%         else
% %             %-----------------Check for Possible Errors--------------------
% %             fid = fopen('/home1/amir/cluster/Clustered_Neural/wav2mat_error.txt','r');
% %                         
% %             if (fid>-1) 
% %                 s = fscanf(fid,'%s');
% %                 fclose(fid);
% %                 if (strcmp(s,'error'))
% %                     error('Something has went wrong!');
% %                 end
% %             end
% %             %--------------------------------------------------------------
%             
%             %-------------If No Errors, Wait Until Convergence-------------
%             pause(20);
%             %--------------------------------------------------------------
%         end
%     end
%     %----------------------------------------------------------------------
    wav2mat_exec(source_depository,offset,max_size,destinaion_depository,destination_db_name,DS_rate,N_out,freq_threshold,duration)
end

%==========================================================================


        
