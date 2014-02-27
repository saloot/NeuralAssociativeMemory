
%%
%==============================INITIALIZATION==============================
addpath(genpath('/home1/amir/cluster/Common_Library')); 
figure_legend = [];                         % This variable will be field with the legend for the graphs

h_error_PER = figure;
% h_error_BER = figure;
set(01,'defaulttextinterpreter','tex')
%---------------------------Fixed Parameters-------------------------------

no_averaging_itrs = 30;
n_exc = 160;
n_inh = 40;
n = n_exc + n_inh;
p = 0.15;
p_plus = p*n_exc/n;
p_minus = p*n_inh/n;
Delta0 = 1;

theta = 0.02;
B_LLR_thr = .5;
tau = 5.06;
q = .65*(1-1/tau)*theta/(p_plus-p_minus);

variable_parameter = 'no_samples';
BP_flag = 1;
inhib_flag = 1;
weight_rule = 1;
B_LLR_flag = 0;
%--------------------------------------------------------------------------

%==========================================================================

%%
%============================PROCESS THE RESULTS===========================
                   
%-----------------------------Read Recall Results--------------------------
if (inhib_flag)
    if (BP_flag)
        if (B_LLR_flag)
            fid = fopen(['Simulation_Results/Belief_LLR_leaky_inhib_n_',num2str(n_exc),'_',num2str(n_inh),'_no_averaging_itrs_',num2str(no_averaging_itrs),'_w_rule_',num2str(weight_rule),'.txt'], 'r');
        else
            fid = fopen(['Simulation_Results/Belief_LLR_leaky_inhib_Delta_1_n_',num2str(n_exc),'_',num2str(n_inh),'_no_averaging_itrs_',num2str(no_averaging_itrs),'_w_rule_',num2str(weight_rule),'.txt'], 'r');
        end
    end
else
    if (BP_flag)
        if (B_LLR_flag)
            fid = fopen(['Simulation_Results/Belief_LLR_leaky_n_',num2str(n),'_no_averaging_itrs_',num2str(no_averaging_itrs),'.txt'], 'r');
        else
            fid = fopen(['Simulation_Results/Belief_LLR_leaky_Delta_1_n_',num2str(n),'_no_averaging_itrs_',num2str(no_averaging_itrs),'.txt'], 'r');
        end
    else
        if (B_LLR_flag)
            fid = fopen(['Simulation_Results/Belief_LLR_n_',num2str(n),'_no_averaging_itrs_',num2str(no_averaging_itrs),'.txt'], 'r');
        else
            fid = fopen(['Simulation_Results/Belief_LLR_Delta_1_n_',num2str(n),'_no_averaging_itrs_',num2str(no_averaging_itrs),'.txt'], 'r');
        end
    end
end
%  TdthetafpfqftaufBBLR_thrfaccfacc_0f
if (fid > -1)            
    results = fscanf(fid, '%s %d %s %f %s %f %s %f %s %f %s %f %s %f %s %f %s %f',[50,inf]);
    fclose(fid);       
else    
    error('Undefined input file!')   
end
%--------------------------------------------------------------------------        

%-----------------------Process the Results--------------------------------    
unprocessed_no_samples = results(2,:);    
unprocessed_theta = results(8,:);    
unprocessed_p = results(10,:);
unprocessed_q = results(12,:);
unprocessed_tau = results(16,:);
unprocessed_BLLR_thr = results(25,:);
unprocessed_acc_plus = results(34,:);
unprocessed_acc_minus = results(44,:);
unprocessed_err = results(50,:);
    



processed_no_samples = [];    
processed_theta = [];
processed_p = [];    
processed_q = [];
processed_tau = [];
processed_acc_plus = [];
processed_acc_minus = [];
processed_err = [];
processed_count = [];
    
for i = 1:length(unprocessed_no_samples)
    if ( (unprocessed_theta(i) == theta) && (unprocessed_p(i) == p) && (norm(unprocessed_q(i)-q)<0.0001) && ( norm(unprocessed_tau(i) - tau)<0.0001) && (unprocessed_BLLR_thr(i) == B_LLR_thr) )
        processed_flag = 0;        
    else
        continue;
    end
        
    for j = 1:length(processed_no_samples)                          
        if (unprocessed_no_samples(i) == processed_no_samples(j))                                
            processed_flag = 1;                                
            break;                        
        end
    end
        
    if (processed_flag == 0)
        processed_no_samples = [processed_no_samples,unprocessed_no_samples(i)];            
        processed_acc_plus = [processed_acc_plus,unprocessed_acc_plus(i)];            
        processed_acc_minus = [processed_acc_minus,unprocessed_acc_minus(i)];            
        processed_err = [processed_err,unprocessed_err(i)];            
        processed_count = [processed_count,1];                
    else        
        processed_acc_plus(j) = processed_acc_plus(j) + unprocessed_acc_plus(i);
        processed_acc_minus(j) = processed_acc_minus(j) + unprocessed_acc_minus(i);
        processed_err(j) = processed_err(j) + unprocessed_err(i);        
        processed_count(j) = processed_count(j) + 1;        
    end   
end

processed_acc_plus = processed_acc_plus./processed_count;
processed_acc_minus = processed_acc_minus./processed_count;
processed_err = processed_err./processed_count;
%--------------------------------------------------------------------------

%----------------------Theoretical Estimations-----------------------------
P_C_p_tot = zeros(1,length(processed_no_samples));
P_C_m_tot = zeros(1,length(processed_no_samples));
P_C_0_tot = zeros(1,length(processed_no_samples));
P_C_tot = zeros(1,length(processed_no_samples));
P_E_tot = zeros(1,length(processed_no_samples));

for itr = 1:1:length(processed_no_samples)    
    TT = processed_no_samples(itr);
    Delta = Delta0/TT;
    [P_C,P_C_plus,P_C_minus,P_C_0,P_E] =  BP_rule_leaky_inhibitory_theory_v2(n,tau,theta,TT,p_plus,p_minus,q,B_LLR_thr,B_LLR_thr,Delta);
    P_C_p_tot(itr) = P_C_plus;
    P_C_m_tot(itr) = P_C_minus;
    P_C_0_tot(itr) = P_C_0;
    P_C_tot(itr) = P_C;
    P_E_tot(itr) = P_E;
end
%--------------------------------------------------------------------------

%-----------------------Constructu the Legend------------------------------
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Sim-Acc+'];    
figure_legend = [figure_legend;temp_leg];    
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Sim-Acc-'];    
figure_legend = [figure_legend;temp_leg];    
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Sim-Acc0'];    
figure_legend = [figure_legend;temp_leg];    
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Thr-Acc+'];    
figure_legend = [figure_legend;temp_leg];    
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Thr-Acc-'];    
figure_legend = [figure_legend;temp_leg];    
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Thr-Acc0'];    
figure_legend = [figure_legend;temp_leg];    
%--------------------------------------------------------------------------
    
%------------------------Display The Results on Graph----------------------
figure(h_error_PER);
% plot(sort(processed_error_bits),sort(PER_tehroy_ave),'b','LineWidth',2);    
[val,ind] = sort(processed_no_samples);

plot(val,processed_acc_plus(ind),'b','LineWidth',2);    
hold on
plot(val,processed_acc_minus(ind),'b--','LineWidth',2);    
hold on
plot(val,processed_err(ind),'r','LineWidth',2);    
hold on
plot(val,P_C_p_tot(ind),'b*','LineWidth',2);    
hold on
plot(val,P_C_m_tot(ind),'b*-','LineWidth',2);    
hold on
plot(val,P_C_0_tot(ind),'r*','LineWidth',2);    
hold on

if (B_LLR_flag)
    title('Belief Propagation')
else
    title('\Delta = 1')    
end
ylabel({'Accuracy'},'interpreter', 'latex','fontsize',24)
xlabel('T','interpreter', 'latex','fontsize',24)
%----------------------------------------------------------------------

%-----------------------Constructu the Legend------------------------------
figure
plot(val,P_C_tot(ind),'b','LineWidth',2);    
hold on
P_C_sim = processed_acc_plus * p_plus + processed_acc_minus * p_minus + (1-p_plus - p_minus) * processed_err;
plot(val,P_C_sim(ind),'b--','LineWidth',2);    
legend('P_C - Thr','P_C - Sim')
%--------------------------------------------------------------------------

%==========================================================================

%======================ADD THE LEGEND TO THE FIGURES=======================
figure(h_error_PER)
legend({figure_legend},'interpreter', 'latex','fontsize',20)
%==========================================================================

