
%%
%==============================INITIALIZATION==============================
addpath(genpath('/home1/amir/cluster/Common_Library')); 
figure_legend = [];                         % This variable will be field with the legend for the graphs

h_error_PER = figure;
% h_error_BER = figure;
set(01,'defaulttextinterpreter','tex')
%---------------------------Fixed Parameters-------------------------------
n = 200;
no_averaging_itrs = 30;

p = 0.3;

theta = 0.05;
B_LLR_thr = 0;
tau = 500000;
q = .65*(1-1/tau)*theta/p;

variable_parameter = 'no_samples';
BP_flag = 1;
%--------------------------------------------------------------------------

%==========================================================================

%%
%============================PROCESS THE RESULTS===========================
                   
%-----------------------------Read Recall Results--------------------------
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
%  TdthetafpfqftaufBBLR_thrfaccferrf
if (fid > -1)            
    results = fscanf(fid, '%s %d %s %f %s %f %s %f %s %f %s %f %s %f %s %f',[33,inf]);                
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
unprocessed_acc = results(29,:);
unprocessed_err = results(33,:);
    



processed_no_samples = [];    
processed_theta = [];
processed_p = [];    
processed_q = [];
processed_tau = [];
processed_acc = [];
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
        processed_acc = [processed_acc,unprocessed_acc(i)];            
        processed_err = [processed_err,unprocessed_err(i)];            
        processed_count = [processed_count,1];                
    else        
        processed_acc(j) = processed_acc(j) + unprocessed_acc(i);
        processed_err(j) = processed_err(j) + unprocessed_err(i);        
        processed_count(j) = processed_count(j) + 1;        
    end   
end

processed_acc = processed_acc./processed_count;
processed_err = processed_err./processed_count;
%----------------------------------------------------------------------                       

%-----------------------Constructu the Legend------------------------------
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Acc'];    
figure_legend = [figure_legend;temp_leg];    
temp_leg = ['$n=',num2str(n),', E=',num2str(no_averaging_itrs),', p=',num2str(p),', q=',num2str(q),', \theta=',num2str(theta),', \tau=',num2str(tau),'$-Err'];    
figure_legend = [figure_legend;temp_leg];    
%--------------------------------------------------------------------------
    
%----------------------Display The Results on Graph--------------------   
figure(h_error_PER);
% plot(sort(processed_error_bits),sort(PER_tehroy_ave),'b','LineWidth',2);    
[val,ind] = sort(processed_no_samples);

plot(val,processed_acc(ind),'b','LineWidth',2);    
hold on
plot(val,processed_err(ind),'r','LineWidth',2);    
hold on

if (B_LLR_flag)
    title('Belief Propagation')
else
    title('\Delta = 1')    
end
ylabel({'Accuracy'},'interpreter', 'latex','fontsize',24)
xlabel('T','interpreter', 'latex','fontsize',24)
%----------------------------------------------------------------------

%==========================================================================

%======================ADD THE LEGEND TO THE FIGURES=======================
figure(h_error_PER)
legend({figure_legend},'interpreter', 'latex','fontsize',20)
%==========================================================================

