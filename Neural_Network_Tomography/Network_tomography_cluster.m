
%===========================INITIALIZATION=================================
B_LLR_thr = 0.5;                  % the threshold usedin updating belief LLRs

theta = 0.02;                   % update threshold  
p = 0.15;                       % connection probability
network_size = 400;             % number of neurons
n_exc = 160;                    % number of excitatory neurons
n_inh = 40;                     % number of inhibotory neurons
n = n_exc + n_inh;
p_plus = p* n_exc/n;
p_minus = p * n_inh/n;
tau = 5.08;                   % This is the rate according to which membrane potential drops
Delta0 = 1;


q = .65*(1-1/tau)*theta/(p_plus-p_minus);      % stimulus firing probability    

no_averaging_itrs = 30;         % number of times we perform each simulation for the sake of averging

    
parameter_range = [500:500:10000];

addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
javaaddpath('/home1/amir/cluster/Common_Library/ganymed-ssh2-build250.jar');
mode = 1;                       % specifies how neural acitivty is generated
weight_rule = 1;                % specifies how the updates will work: additive (1) or multlipicative (2)
B_LLR_flag = 1;
%------------------------Initialize the SSH Connection---------------------
channel = sshfrommatlab('amir','lth.epfl.ch','arashmidos'); 
%--------------------------------------------------------------------------

%==========================================================================

%========================CALCULATING BEST THRESHOLDS=======================
t_fire = 1/(1-1/tau);                           % The average weighting time between spikes
c_1 = (1-(1/tau^(t_fire + 1)) )/(1-(1/tau));
c_2 = (1-(1/tau^(2*t_fire + 2)) )/(1-(1/tau^2));

% P_A_plus = q_function(B_LLR_thr,q*c_1,sqrt(q*(1-q)*c_2));            % probability that the net input from a particular neuron is larger than B_LLR_thr
% P_A_minus = q_function(B_LLR_thr,q*c_1,sqrt(q*(1-q)*c_2));          % probability that the net input from a particular neuron is larger than B_LLR_thr
P_A_plus = q;                                   % probability that the net input from a particular neuron is larger than B_LLR_thr
P_A_minus = q;                                  % probability that the net input from a particular neuron is larger than B_LLR_thr
c_1 = 1;
c_2 = 1;

mu = c_1*(n-1)*q*(p_plus - p_minus);                                  
sigma = sqrt(c_2*(n-1)*q*(p_plus*(1-q*p_plus)+p_minus*(1-q*p_minus)));

p_inc_p = P_A_plus*q_function(theta * n - B_LLR_thr,mu,sigma);             % Probability of increasing the weight on the condition that Gi = 1
p_inc_m = P_A_plus*q_function(theta * n + B_LLR_thr,mu,sigma);             % Probability of increasing the weight on the condition that Gi = -1
p_inc_0 = P_A_plus*q_function(theta * n ,mu,sigma);                         % Probability of increasing the weight on the condition that Gi = 0

p_dec_p = P_A_minus*(1-q_function(theta * n - B_LLR_thr,mu,sigma));       % Probability of decreasing the weight on the condition that Gi = 1
p_dec_m = P_A_minus*(1-q_function(theta * n + B_LLR_thr,mu,sigma));       % Probability of decreasing the weight on the condition that Gi = -1
p_dec_0 = P_A_minus*(1-q_function(theta * n ,mu,sigma));                    % Probability of decreasing the weight on the condition that Gi = 0

%==========================================================================


%======================LOOP OVER THE PARAMETERS============================
for itr = 1:length(parameter_range)
    no_samples = parameter_range(itr);       % number of sample recordings
    T = no_samples;
    Delta = Delta0/T;
    
    mu1 = Delta*T*(p_inc_p-p_dec_p);
    var1 = Delta*sqrt(T*(p_inc_p*(1-p_inc_p) + p_dec_p*(1-p_dec_p)));

    mu2 = Delta*T*(p_inc_m-p_dec_m);
    var2 = Delta*sqrt(T*(p_inc_m*(1-p_inc_m) + p_dec_m*(1-p_dec_m)));

    mu3 = Delta*T*(p_inc_0-p_dec_0);
    var3 = Delta*sqrt(T*(p_inc_0*(1-p_inc_0) + p_dec_0*(1-p_dec_0)));

    if (weight_rule == 1)
        epsilon_plus = (mu1+mu3)/2;
        epsilon_minus = abs(mu2+mu3)/2;
    else
        a = 0.5*(1+Delta)^(mu1/Delta);
        b = 0.5*(1+Delta)^(mu2/Delta);
        c = 0.5*(1+Delta)^(mu3/Delta);
        epsilon_plus = (a+c)/2;
        epsilon_minus = (b+c)/2;
    end
    %--------------------Submit the Job to the Cluster---------------------    
    
%     command = ['cd /scratch/amir/Network_Tomography/Submit_Cluster;qsub -N "Net_tomo_',num2str(itr),...
%         '" -v nn=',num2str(network_size),',T=',num2str(no_samples),',tau=',num2str(tau),',theta=',num2str(theta),...                                        
%         ',p=',num2str(p),',q=',num2str(q),',E=',num2str(no_averaging_itrs),',B_LLR_thr=',num2str(B_LLR_thr),',B_LLR_flag=',num2str(B_LLR_flag),' net_tomography.pbs'];                                            
%     command = ['cd /scratch/amir/Network_Tomography/Submit_Cluster;qsub -N "Net_tomo_',num2str(itr),...
%         '" -v nn=',num2str(network_size),',T=',num2str(no_samples),',tau=',num2str(tau),',theta=',num2str(theta),...                                        
%         ',p=',num2str(p),',q=',num2str(q),',E=',num2str(no_averaging_itrs),',B_LLR_thr=',num2str(B_LLR_thr),',B_LLR_flag=',num2str(B_LLR_flag),' net_tomography_BP.pbs'];
    command = ['cd /scratch/amir/Network_Tomography/Submit_Cluster;qsub -N "Net_tomo_',num2str(itr),...
        '" -v nn_exc=',num2str(n_exc),',nn_inh=',num2str(n_inh),',T=',num2str(no_samples),',tau=',num2str(tau),',theta=',num2str(theta),...                                        
        ',p=',num2str(p),',q=',num2str(q),',E=',num2str(no_averaging_itrs),',B_LLR_thr=',num2str(B_LLR_thr),',B_LLR_flag=',num2str(B_LLR_flag),...
        ',mode=',num2str(mode),',weight_rule=',num2str(weight_rule),',Delta=',num2str(Delta),',epsilon_plus=',num2str(epsilon_plus),',epsilon_minus=',num2str(epsilon_minus),' net_tomography_BP_inhib.pbs'];
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