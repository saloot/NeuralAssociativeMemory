%==========================================================================
%**************************FUNCTION: wav2mat****************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% source_depository: The address of the database of original images, stored as rows of a big matrix (each image is mapped from a 2D array to a one-dimensional vector)
% offset: It is the relative place of the file which the function should start processing. It is measured from the firt file in the folder
% max_size: It is the maximum number of files that should be processed by the code
% destination_file: A string containing the address of the file which will store the produced datasets. Example: '/home/hesam/test_file.mat'
% destination_db_name: Specifies the name of the MATLAB variable within the .mat file that holds the processed files up to this point
% DS_rate: If DS_rate = 1, the function also rescales the images and map them to gray scale levels.
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


function wav2mat_exec(source_depository,offset,max_size,destinaion_depository,destination_db_name,DS_rate,N_out,freq_threshold,duration)


%===========================INITIALIZATION=================================

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
%--------------------------------------------------------------------------

%--------------------------Image Specifications----------------------------
% max_level = 16;                                                       % This is the maximum number of gray-levels we would like to have in the produced patches
listing = dir([source_depository,'/*.wav']);                        % This is the ist of images in the depository with the specified folder

if (isempty(listing))
    fid = fopen('/home1/amir/cluster/Clustered_Neural/img2mat_error.txt','w');
    fprintf(fid,'error');
    fclose(fid);
    error('Depository is empty!');
end

if (~exist(destinaion_depository,'dir'))
    mkdir(destinaion_depository);
end

destinaion_file = [destinaion_depository,'/',destination_db_name,'.mat'];
destinaion_file_DS = [destinaion_depository,'/DS/',destination_db_name,'.mat'];
destinaion_file_FFT = [destinaion_depository,'/FFT/',destination_db_name,'.mat'];
destinaion_file_DS_FFT = [destinaion_depository,'/DS/FFT/',destination_db_name,'.mat'];
load(destinaion_file);
%--------------------------------------------------------------------------

temp_DS_db = [];
temp_fft_db = [];
temp_DS_fft_db = [];

%==========================================================================



%======================READ FILES AND CONVERT THEM=========================
for i = offset+1:min(offset+max_size,length(listing))
    file_name = listing(i).name;
    file_name_full = [source_depository,'/',file_name];
    WAV = wavread(file_name_full);
    WAV = [WAV',zeros(1,duration-length(WAV))];
    WAV = WAV(1:duration);
    
    DS_WAV = downsample(WAV,DS_rate);                                        
    
    temp_DS_db = [temp_DS_db;DS_WAV];                                   % Store the downsampled version
    temp_fft_db = [temp_fft_db; wav2freq(WAV,N_out,freq_threshold)];                               % Store the frequency response of the sounds
    temp_DS_fft_db = [temp_DS_fft_db; wav2freq(DS_WAV,N_out,freq_threshold)];                      % Store the frequency response of the downsampled sound
    
    eval([destination_db_name '= [',destination_db_name,';WAV];']);                
end
%==========================================================================

%===========================STORE THE RESULTS==============================
save(destinaion_file,destination_db_name);        
eval(['clear ' destination_db_name]);

load(destinaion_file_DS);
eval([destination_db_name '= [',destination_db_name,';temp_DS_db];']);                
save(destinaion_file_DS,destination_db_name);        
eval(['clear ' destination_db_name]);

load(destinaion_file_FFT);
eval([destination_db_name '= [',destination_db_name,';temp_fft_db];']);                
save(destinaion_file_FFT,destination_db_name);        
eval(['clear ' destination_db_name]);

load(destinaion_file_DS_FFT);
eval([destination_db_name '= [',destination_db_name,';temp_DS_fft_db];']);                
save(destinaion_file_DS_FFT,destination_db_name);        
%==========================================================================


% fid = fopen('/home1/amir/cluster/Clustered_Neural/wav2mat_error.txt','w');
% fprintf(fid,' ');
% fclose(fid);
%==========================================================================

