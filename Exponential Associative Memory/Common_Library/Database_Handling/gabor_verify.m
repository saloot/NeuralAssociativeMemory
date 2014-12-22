%==========================================================================
%*************************FUNCTION: gabor_verify***************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% file_name_base: The base of the file names, which are supposed to be indexed according to file_i_j. In this case, file_name_base = '/home1/.../file_'.
% L_horiz: The maximum number of horizontal shifts (the first index)
% L_VERT: The maximum number of vertical shifts (the second index)
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% success_index: a 0/1 matrix where a 1 in position (i,j) of the matrix indicaes that the file corresponding to the portion (i,j) of the image exists in the database.
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the base of image databases and checks if the
% corresponding .mat file corresponding to the various parts of the images
% exsit in the database or not. The function returns a 0/1 matrix where a 1
% in position (i,j) of the matrix indicaes that the file corresponding to
% the portion (i,j) of the image exists in the database.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function success_index = gabor_verify(file_name_base,L_horiz,L_vert)

%===========================INITIALIZATION=================================

%-------------------------Add Necessary Libraries--------------------------
addpath('/home1/amir/cluster/Common_Library');                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
%--------------------------------------------------------------------------

success_index = zeros(L_horiz,L_vert);

%==========================================================================


%==============================MAIN LOOP===================================
for index_horiz = 1:L_horiz
    for index_vert = 1:L_vert
%         for i = 1:32
%             file_name = [file_name_base,'/Patch_',num2str(index_horiz),'_',num2str(index_vert),'/Gabor_out_',num2str(i),'.txt'];        
% %         listing = dir([file_name_base]);  
% %         fid = fopen(file_name);
%             if (exist(file_name,'file'))
%                 success_index(index_horiz,index_vert) = 1;
%             else
%                 success_index(index_horiz,index_vert) = 0;
%                 break;
%             end
%         end
        file_name = [file_name_base,'/Patch_',num2str(index_horiz),'_',num2str(index_vert),'/Features_vector.txt'];
        fid = fopen(file_name,'r');
        if (fid>-1)
        dataset_learn = fscanf(fid, '%f',[32,inf]);
        dataset_learn = dataset_learn';        
        fclose(fid);
        
        [dataset_size,~] = size(dataset_learn);
        if (dataset_size == 500)
            success_index(index_horiz,index_vert) = 1;
        else
            success_index(index_horiz,index_vert) = 0;
        end
        end
    end
end
%==========================================================================