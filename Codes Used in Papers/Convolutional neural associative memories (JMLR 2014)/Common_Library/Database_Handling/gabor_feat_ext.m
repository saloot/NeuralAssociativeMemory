%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                 % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
%--------------------------------------------------------------------------

%--------------------Initialize the SSH Connection-------------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');
%--------------------------------------------------------------------------

L_horiz = 19;
L_vert = 19;
file_name_base = ['/scratch/amir/Databases/STL_10/L_horiz_',num2str(L_horiz),'_L_vert_',num2str(L_vert)];
success_index = gabor_verify(file_name_base,L_horiz,L_vert);

for l_h = 1 : L_horiz
    for l_v = 1 : L_vert
        if (success_index(l_h,l_v) == 0)
            l_h
            l_v
        command = ['cd /scratch/amir/Databases/Submit_Cluster;qsub -N "job_gabor_',num2str(l_h),'_',num2str(l_v),'" -v L_horiz=',...
                   num2str(L_horiz),',L_vert=',num2str(L_vert),',l_h=',num2str(l_h),',l_v=',num2str(l_v),' gabor_find.pbs'];        
            [channel, result]  =  sshfrommatlabissue(channel,command);             
            %----------Check the success of the submission---------
            if (isequal(result, {''}))     
                display('Unsubmitted Gabor job!');        
            else            
            %--------------------------------------------------                        
                display('Gabor job submitted successfully!');            
            end    
            %------------------------------------------------------
    
        end
    end
end

% file_name_base = '/scratch/amir/Databases/STL_10/L_horiz_22_L_vert_22/Class_1';
% success_index = gabor_verify(file_name_base,L_horiz,L_vert);

channel  =  sshfrommatlabclose(channel);


