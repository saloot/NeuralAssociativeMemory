addpath(genpath('/home/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');

channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos');

num_train = 1000;
num_test = 1000;
training_method = 1;
test_method = 2;
no_dims = 50;


for err_bits = 1:4
    for pattern_neur_noise = 0:.1:.8
        command = ['cd /scratch/amir/Fault_Tolerant_Decoding/Submit_Cluster;qsub -N "P_C_',num2str(err_bits),'_',num2str(pattern_neur_noise),'" -v err_bits='...
            ,num2str(err_bits),',pattern_neur_noise=',num2str(pattern_neur_noise),', Calculate_P_c_Empirical.pbs'];

        [channel, result]  =  sshfrommatlabissue(channel,command);
        result
    end
end
channel  =  sshfrommatlabclose(channel);