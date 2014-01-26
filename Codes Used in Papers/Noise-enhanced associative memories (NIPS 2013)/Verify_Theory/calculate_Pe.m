%==========================================================================
%******************FUNCTION: calculate_Pe_theoretical**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% error_bits_in: The number of external errors to the network
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% const_neur_noise: The maximum amount of noise a constraint neuron will "suffer" from
% gamma: The update threshold in the original bit-flipping recall algorithm 
% index: The index of the simulation setup among various random scenarios

%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% f: The upper bound on BER of the inter-cluster error correcing algorithm
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates a theoretical upper bound on the probability of
% error in a clustered neural networks built out of faulty neurons. The
% bound is explained in our ISIT 2013 paper.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

pattern_neur_noise_step = 0.01;
pattern_neur_noise_range = [0:0.001:.01,.02:pattern_neur_noise_step:1];
L1 = length(pattern_neur_noise_range);
n = 40;
m = 20;
cd '..'
[~,deg_column,lambda,rho] = read_cluster_degree(n,m,50,0.95,0.75,0.05);
cd './Verify_Theory'

Pai2_1 = 0;
% 
gamma_step = 0.01;
gamma_range = [0:gamma_step:1];
L2 = length(gamma_range);

P2_2 = zeros(L2,L1);
P2_1 = zeros(L2,L1);
for i = 1:L1
    pattern_neur_noise = pattern_neur_noise_range(i);
    for j = 1:L2
        gamma = gamma_range(j);
        temp1 = 0;
        temp2 = 0;
        for k = 1:length(lambda)
            if (lambda(k) == 0)
                continue;
            else
                d1 = deg_column(k);
            end
            
            temp2_1 = 0;            
            p_neighbor = d1/m;
            
            temp1_1 = P2_2_calculation(gamma,pattern_neur_noise,d1,Pai2_1);
            
            for l = 1:length(lambda)
                if (lambda(l) == 0)
                    continue;
                else
                    deg = deg_column(l);
                end
               
                temp2_1 = temp2_1 + lambda(l) * P2_1_calculation(gamma,pattern_neur_noise,deg,d1,Pai2_1,p_neighbor);
            end
            temp1 = temp1 + lambda(k) * temp1_1;
            temp2 = temp2 + lambda(k) * temp2_1;
        end
        P2_2(j,i) = temp1;
        P2_1(j,i) = temp2;
    end
end

% 
% plot(gamma_range,P2_2(:,100));
% hold on
% plot(gamma_range,P2_2(:,90),'r');
% plot(gamma_range,P2_2(:,60),'g');
% plot(gamma_range,P2_2(:,40),'k');
% plot(gamma_range,P2_2(:,20),'c');
% plot(gamma_range,P2_2(:,10),'y');
% 
% xlabel('\gamma')
% ylabel('Pe^{(2)}(\pattern_neur_noise)');



Pe_tot = 1- (1-P2_1).^(n-1) .* (1-P2_2);
Pc_tot = 1-Pe_tot;
% Pe_tot = (1-1/n)*P2_1 + (1/n) * P2_2;


figure


plot(gamma_range,Pc_tot(:,105),'r');
hold on
plot(gamma_range,Pc_tot(:,80),'g');
plot(gamma_range,Pc_tot(:,50),'b');
plot(gamma_range,Pc_tot(:,30),'k');
plot(gamma_range,Pc_tot(:,20),'c');
plot(gamma_range,Pc_tot(:,6),'y');
xlabel('\gamma')
ylabel('P_e(\pattern_neur_noise)');
legend_str = ['\upsilon =  ',num2str(pattern_neur_noise_range(105))];
legend_str = [legend_str;'\upsilon =   ',num2str(pattern_neur_noise_range(80))];    
legend_str = [legend_str;'\upsilon =   ',num2str(pattern_neur_noise_range(50))];    
legend_str = [legend_str;'\upsilon =   ',num2str(pattern_neur_noise_range(30))];    
legend_str = [legend_str;'\upsilon =   ',num2str(pattern_neur_noise_range(20))];    
legend_str = [legend_str;'\upsilon = ',num2str(pattern_neur_noise_range(6))];    
legend(legend_str)

gamma_min = zeros(1,L1);
Pe_min = zeros(1,L1);

for i = 1:L1
    Pe = Pe_tot(:,i);
    [mm,ind] = min(Pe);
    gamma_min(i) = gamma_range(ind);
    Pe_min(i)=mm;
end

figure
plot(pattern_neur_noise_range,gamma_min)

figure
plot(pattern_neur_noise_range,Pe_min)
% 
% save(['Pro/(2*eb_e_pai_',num2str(Pai2_1),'.mat']);

% P1_bound = zeros(1,L1);
% P1_exact = zeros(1,L1);
% for i = 1:L1
%     
%     pattern_neur_noise = pattern_neur_noise_range(i);
%     gamma = gamma_min(i);
%     if (gamma+pattern_neur_noise<1)
%         continue;
%     end
%     temp = 2*(1-d_ave/m)*(((d_ave/m)/(1-d_ave/m))^(gamma-pattern_neur_noise));
%     P1_bound(i) = (1-gamma+pattern_neur_noise)*((min(1,temp))^d1)/(2*pattern_neur_noise);   
% %     temp =(1-d_ave/m)*(((d_ave/m)/(1-d_ave/m))^(gamma-pattern_neur_noise))*(exp(1)/max(.5,gamma-pattern_neur_noise))^max(.5,gamma-pattern_neur_noise);
% %     P1_bound(i) = d1*(1-gamma+pattern_neur_noise)^2*((min(1,temp))^d1)/(2*pattern_neur_noise);   
%     P1_exact(i) = P2_1_calculation(gamma,pattern_neur_noise,d1,Pai2_1,d_ave,m);
% end
% 
% 
% figure
% plot(pattern_neur_noise_range,P1_bound);
% hold on
% plot(pattern_neur_noise_range,P1_exact,'r');