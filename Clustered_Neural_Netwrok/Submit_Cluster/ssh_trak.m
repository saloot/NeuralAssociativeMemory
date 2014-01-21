addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');

load('/scratch/amir/ssh_connections.mat');        
[L,~] = size(connection_list);

for l = 1:L
    channel = connection_list(l);
    channel  =  sshfrommatlabclose(channel);
end
connection_list = [];
save('/scratch/amir/ssh_connections.mat','connection_list')
        
