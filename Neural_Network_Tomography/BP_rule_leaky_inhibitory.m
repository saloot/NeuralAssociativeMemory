function BP_rule_leaky_inhibitory(n_exc,n_inh,T,tau,theta,p,q,no_averaging_itrs,B_LLR_thr,B_LLR_flag)

%===========================INITIALIZATION=================================
n = n_exc + n_inh;
acc_ave_plus = 0;                    % average number of correctly infered excitatory edges
acc_ave_minus = 0;                    % average number of correctly infered inhibitory edges
err_ave = 0;                    % average number of mistakenly infered edges + missed edges
        
avg_no_edge = round(p*n);       % estimated number of edges in the graph
no_edges_fire = ceil(theta*n);  % estimated number of edges required for a neuron to fire

success_measure = zeros(1,3);
acc_theory = 0;
err_theory = 0;
% fit_coeff_0 = [-0.0000    0.0002   -0.0041    0.0362   -0.1104];
% fit_coeff_1 = [0.0001   -0.0073    0.2498   -3.8715   24.0011];

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%==========================================================================


%=======================AVERAGE OVER AN ENSEMBLE===========================
for av_itr = 1:no_averaging_itrs        
    
    t_last = 0;                     % The last time the output neuron fired
        
    %----------------------Create an Erdos Reny Random Graph---------------    
    G = erdos_reny_inhibitory(1,n_exc,n_inh,p);        
    %----------------------------------------------------------------------
                
    %------------------------Record Sample States--------------------------
    S = zeros(T,n);                 % recorded stimulus states         
    R = zeros(T,1);                 % recorder output state            
                
    for i = 1:T    
        s = max(rand(1,n) < q,0);   % stimulus applied to n neurons            
        S(i,:) = s';                                
        denom = tau.^(0:1:i-t_last-1);                        
        h = sum(G*(S(i:-1:t_last+1,:)')./denom);
            
        
        R(i) = 1*(h/n > theta);        
        if (R(i))
            t_last = i;
        end
    end    
    %----------------------------------------------------------------------
                
    %--------------------Infer the Connectivity Matrix---------------------        
    W = zeros(1,n);
    pp = p *(n_exc-n_inh)/n;
    p_plus = p * n_exc/n;
    p_minus = p * n_inh/n;
    p_var = p_plus * (1-p_plus) + p_minus * (1-p_minus) + 2 *p_plus * p_minus;
    for i = 1:n       
        t_last = 0;
        for jj = 1:T                            
            temp = 0;            
            temp_L = 0;
            temp_var = 0;
            temp_var_L = 0;
            for j = jj:-1:t_last+1                    
                temp = temp + S(j,i)/tau^(1*(jj-j));
                L = sum(S(j,:));                
                temp_L = temp_L + L/tau^(1*(jj-j));
                temp_var = temp_var + S(j,i)/tau^(2*(jj-j));
                temp_var_L = temp_var_L + L/tau^(2*(jj-j));
            end
                
            mu_plus = (1-pp)*temp + pp * temp_L;
            mu_zero = (-pp)*temp + pp * temp_L;
            mu_minus = (-pp-1)*temp + pp * temp_L;
            
            sigma_plus = sqrt(p_var*(temp_var_L - temp_var));

            % sigma_zero = sigma_plus;
            sigma_zero = sqrt(p_var*(temp_var_L));

            % sigma_minus = sigma_plus;
            sigma_minus = sqrt(p_var*(temp_var_L + temp_var));
            
            Q_plus = q_function(no_edges_fire,mu_plus,sigma_plus);
            Q_zero = q_function(no_edges_fire,mu_zero,sigma_zero);
            Q_minus = q_function(no_edges_fire,mu_minus,sigma_minus);
            
            Q_plus = R(jj)*Q_plus+(1-R(jj))*(1-Q_plus);
            Q_minus = R(jj)*Q_minus+(1-R(jj))*(1-Q_minus);
            Q_zero = R(jj)*Q_zero+(1-R(jj))*(1-Q_zero);
            
            
            Q = Q_plus + Q_zero + Q_minus;
            B_LLR1 = log(Q_plus/((Q-Q_plus)/2));
            B_LLR2 = log(Q_zero/((Q-Q_zero)/2));
            B_LLR3 = log(Q_minus/((Q-Q_minus)/2));
            B_LLR = [B_LLR1,B_LLR2,B_LLR3];
                
            if ( (Q_plus > Q_zero) && (Q_plus > Q_minus) ) 
                val = B_LLR1;
            elseif ( ( Q_minus > Q_plus) && ( Q_minus > Q_zero) ) 
                val = -abs(B_LLR3);
            else
                val = -abs(B_LLR2)*sign(W(i))/2;
            end
                
            if (val>B_LLR_thr)                            
                if (B_LLR_flag)                                
                    W(i) = W(i) + val * temp;                            
                else                        
                    W(i) = W(i) + 1 * temp;
                end                    
            end                
%             else      
%                 Q_plus = 1-Q_plus;
%                 Q_minus = 1-Q_minus;
%                 Q_zero = 1-Q_zero;
%                 Q = Q_plus + Q_zero + Q_minus;
%                 
%                 B_LLR1 = log(Q_plus/((Q-Q_plus)/2));
%                 B_LLR2 = log(Q_zero/((Q-Q_zero)/2));
%                 B_LLR3 = log(Q_minus/((Q-Q_minus)/2));
%                 B_LLR = [B_LLR1,B_LLR2,B_LLR3];
%                 
%                 if ( (Q_plus > Q_zero) && (Q_plus > Q_minus) ) 
%                     val = B_LLR1;
%                 elseif ( ( Q_minus > Q_plus) && ( Q_minus > Q_zero) ) 
%                     val = -abs(B_LLR3);
%                 else
%                     val = 0;
%                 end
%                     
%                 
%                 if (val>B_LLR_thr)                            
%                     if (B_LLR_flag)                                
%                         W(i) = W(i) + val * temp;                            
%                     else                        
%                         W(i) = W(i) + 1 * temp;
%                     end                    
%                 end                
%             end
            if (R(jj))                    
                t_last = jj;                            
            end
        end        
    end    
    W2 = W;        
    %----------------------------------------------------------------------

    %--------------------Determine the Weight Threshold--------------------        
    [W_sorted,ind] = sort(W);
     % w_thr = determine_weight_thr_leaky(n,p,q,theta);
     % w_thr = q*(1/avg_T_fire+1/avg_T_fire2)/2;
     % W = (W>w_thr);        
     
     W = zeros(1,n);     
     W(ind(end-sum(G>0)+1:end)) = abs(sign(W_sorted(end-sum(G>0)+1:end)));        
     W(ind(1:sum(G<0))) = -abs(sign(W_sorted(1:sum(G<0))));        
     %---------------------------------------------------------------------
                
     %----------------------Calculate the Accuracy-------------------------
     acc_plus = sum((G>0).*(W>0))/sum(G>0);
     acc_minus = sum((G<0).*(W<0))/sum(G<0);
     % err = abs(sum(sign(abs(W)+abs(G))-((G==W).*(abs(G)>0))));
     err = sum(W~=G);
     acc_ave_plus = acc_ave_plus + acc_plus;
     acc_ave_minus = acc_ave_minus + acc_minus;
     
     err_ave = err_ave + err/(n-sum(abs(G)));        
     %---------------------------------------------------------------------
    
end
%==========================================================================


%===========================SAVE THE RESULTS===============================
success_measure(1) = acc_ave_plus/av_itr;
success_measure(2) = acc_ave_minus/av_itr;
success_measure(3) = err_ave/av_itr;   

% [p1_T,p2_T] = simple_hebb_rule_theory_leaky(n,p,q,T,theta,tau,avg_T_fire);
% simple_hebb_rule_theory(n,p,q,T,theta,0);    
% acc_theory(itr) =p1_T;    
% err_theory(itr) = p2_T;

if (B_LLR_flag)
    fid = fopen(['Simulation_Results/Belief_LLR_leaky_inhib_n_',num2str(n_exc),'_',num2str(n_inh),'_no_averaging_itrs_',num2str(no_averaging_itrs),'.txt'], 'a+');        
else
    fid = fopen(['Simulation_Results/Belief_LLR_leaky_inhib_Delta_1_n_',num2str(n_exc),'_',num2str(n_inh),'_no_averaging_itrs_',num2str(no_averaging_itrs),'.txt'], 'a+');        
end
fprintf(fid, 'T \t %d \t theta \t %f \t p \t %f \t q \t %f \t tau \t %f \t BBLR_thr \t %f \t acc_plus \t %f \t acc_minus \t %f \t err \t %f \t',T,theta,p,q,tau,B_LLR_thr,success_measure(1),success_measure(2),success_measure(3));
fprintf(fid,'\n');
fclose(fid);
    
%==========================================================================

