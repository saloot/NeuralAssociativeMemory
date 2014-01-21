function qdel(start_id,end_id)

javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');

%--------------------Initialize the SSH Connection-------------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');
%--------------------------------------------------------------------------
for i = start_id:end_id
    command = ['qdel ',num2str(i)];            
    [channel, result]  =  sshfrommatlabissue(channel,command); 
%     if (isequal(result, {''}))                
%         display('Unsubmitted learning job!');        
%     else        
%         display('Learning job submitted successfully!');        
%     end
%     display(result)
end
