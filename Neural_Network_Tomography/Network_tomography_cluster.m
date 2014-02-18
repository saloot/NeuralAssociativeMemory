
%===========================INITIALIZATION=================================
B_LLR_thr = -100000;                  % the threshold usedin updating belief LLRs
B_LLR_flag = 1;                 
theta = 0.05;                   % update threshold  
p = 0.3;                        % connection probability
network_size = 200;             % number of neurons
n_exc = 150;                    % number of excitatory neurons
n_inh = 30;                     % number of inhibotory neurons
tau = 5.02;                   % This is the rate according to which membrane potential drops

q = .75*(1-1/tau)*theta/p;      % stimulus firing probability    

no_averaging_itrs = 30;         % number of times we perform each simulation for the sake of averging

    
parameter_range = [700:500:5100];

addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
mode = 2;

%------------------------Initialize the SSH Connection---------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos'); 
%--------------------------------------------------------------------------

%==========================================================================



%======================LOOP OVER THE PARAMETERS============================
for itr = 1:length(parameter_range)
    no_samples = parameter_range(itr);       % number of sample recordings
    
    %--------------------Submit the Job to the Cluster---------------------    
    
%     command = ['cd /scratch/amir/Network_Tomography/Submit_Cluster;qsub -N "Net_tomo_',num2str(itr),...
%         '" -v nn=',num2str(network_size),',T=',num2str(no_samples),',tau=',num2str(tau),',theta=',num2str(theta),...                                        
%         ',p=',num2str(p),',q=',num2str(q),',E=',num2str(no_averaging_itrs),',B_LLR_thr=',num2str(B_LLR_thr),',B_LLR_flag=',num2str(B_LLR_flag),' net_tomography.pbs'];                                            
%     command = ['cd /scratch/amir/Network_Tomography/Submit_Cluster;qsub -N "Net_tomo_',num2str(itr),...
%         '" -v nn=',num2str(network_size),',T=',num2str(no_samples),',tau=',num2str(tau),',theta=',num2str(theta),...                                        
%         ',p=',num2str(p),',q=',num2str(q),',E=',num2str(no_averaging_itrs),',B_LLR_thr=',num2str(B_LLR_thr),',B_LLR_flag=',num2str(B_LLR_flag),' net_tomography_BP.pbs'];
    command = ['cd /scratch/amir/Network_Tomography/Submit_Cluster;qsub -N "Net_tomo_',num2str(itr),...
        '" -v nn_exc=',num2str(n_exc),',nn_inh=',num2str(n_inh),',T=',num2str(no_samples),',tau=',num2str(tau),',theta=',num2str(theta),...                                        
        ',p=',num2str(p),',q=',num2str(q),',E=',num2str(no_averaging_itrs),',B_LLR_thr=',num2str(B_LLR_thr),',B_LLR_flag=',num2str(B_LLR_flag),',mode=',num2str(mode),' net_tomography_BP_inhib.pbs'];
    [channel, result]  =  sshfrommatlabissue(channel,command);                                                 
    %----------------------------------------------------------------------
                                
                                        
    %------------------Check the success of the submission-----------------                                    
    if (isequal(result, {''}))                                         
        display('Unsubmitted job!');
    else        
        display('Job submitted successfully!');                                    
    end 
end
%==========================================================================
    
%=========================CLOSE THE SSH SESSION============================
channel  =  sshfrommatlabclose(channel);
%==========================================================================


% %==============================PLOT RESULTS================================
% figure
% plot(parameter_range,success_measure(1,:),'g-*')
% hold on
% plot(parameter_range,success_measure(2,:),'k-o')
% plot(parameter_range,acc_theory,'r')
% plot(parameter_range,err_theory,'k')
% title(['q=',num2str(q),' p=',num2str(p),' theta=',num2str(theta),' n=',num2str(n)])
% legend('No. correct edges','No. false edges')
% xlabel('T')
% %==========================================================================