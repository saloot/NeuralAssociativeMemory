function simple_Hebb_rule_leaky(n,T,tau,theta,p,q,no_averaging_itrs

%===========================INITIALIZATION=================================
acc_ave = 0;                    % average number of correctly infered edges
err_ave = 0;                    % average number of mistakenly infered edges + missed edges
        
avg_no_edge = round(p*n);       % estimated number of edges in the graph
no_edges_fire = ceil(theta*n);  % estimated number of edges required for a neuron to fire

success_measure = zeros(1,2);
acc_theory = 0;
err_theory = 0;

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%==========================================================================


%=======================AVERAGE OVER AN ENSEMBLE===========================
for av_itr = 1:no_averaging_itrs        
    
    t_last = 0;                     % The last time the output neuron fired
        
    %----------------------Create an Erdos Reny Random Graph---------------    
    G = erdos_reny(1,n,p);        
    %----------------------------------------------------------------------
                
    %------------------------Record Sample States--------------------------
    S = zeros(T,n);                 % recorded stimulus states     
    S2 = zeros(T,n);                 % recorded stimulus states     
    R = zeros(T,1);                 % recorder output state    
    R2 = zeros(T,1);                 % recorder output state    
    t_fire = zeros(T,1);            % recorder fie time for the output neuron
                
    for i = 1:T    
        s = max(rand(1,n) < q,0);   % stimulus applied to n neurons            
        S(i,:) = s';                                
        denom = tau.^(0:1:i-t_last-1);                        
        h = sum(G*(S(i:-1:t_last+1,:)')./denom);
            
        R(i) = 1*(h/n > theta);
        R2(i) = R(i);                                                                                
    end    
    %----------------------------------------------------------------------
        
        %------------------Infer the Connectivity Matrix-------------------
        W = zeros(1,n);
        R2 = R;
        ss = sum(S');
        ss = ss-no_edges_fire;
        ss = ss.*(ss>0);
        ss = ss + .2;
        
        R = R - (R==0);
%         R2 = R./ss';
                
%         S = S - (S==0)./(sum(S')'*ones(1,n))/p;        
%         for i=1:n            
%             W(i) = sum(S(:,i)'*R)/T;
% %             temp = sum(S(:,i)'*t_fire)/T;
%         end
    for i = 1:n
        t_last = 0;
        for jj = 1:T
            if (S(jj,i))
                temp = 0;
                for j = jj:-1:t_last+1
%                     temp = temp+ ((-1)^(R2(jj)+S(j,i)))/tau2^(1*(jj-j));
                    temp = temp + S(j,i)/tau^(1*(jj-j));
                end
                L = sum(S(jj,:));
                if (L > no_edges_fire)                    
%                     B_LLR = calculate_belief_LLR(no_dges_fire,L,p);
                    
                    
                    if (R2(jj))
                        B_LLR = polyval(fit_coeff_1,L);
%                         if (B_LLR>1.5)
                            W(i) = W(i) + temp;
%                         end
                    else
                        B_LLR = polyval(fit_coeff_0,L);
%                         if (B_LLR>.15)
                            W(i) = W(i) - temp;
%                         end
                    end
                    
                end
                if (R2(jj))
                    t_last = jj;            
                end
            end
        end
        111;
    end
        W2 = W;
        %------------------------------------------------------------------


        %------------------Determine the Weight Threshold------------------
        [W_sorted,ind] = sort(W);
%         w_thr = determine_weight_thr_leaky(n,p,q,theta);
%         w_thr = q*(1/avg_T_fire+1/avg_T_fire2)/2;
%         W = (W>w_thr);
        W = zeros(1,n);
        W(ind(end-sum(G):end)) = abs(sign(W_sorted(end-sum(G):end)));
        %------------------------------------------------------------------

        
        %--------------------Calculate the Accuracy-------------------------
        acc = sum(W.*G);
        err = abs(sum(sign(W+G)-W.*G));
        acc_ave = acc_ave + acc/sum(G);
        err_ave = err_ave + err/(n-sum(G));
        %------------------------------------------------------------------

    end
    
    %-----------Calculate Theoretical Bounds and Store Results-------------
    success_measure(1,itr) = acc_ave/av_itr;
    success_measure(2,itr) = err_ave/av_itr;    
%     [p1_T,p2_T] = simple_hebb_rule_theory_leaky(n,p,q,T,theta,tau,avg_T_fire);
    simple_hebb_rule_theory(n,p,q,T,theta,0);
    acc_theory(itr) =p1_T;
    err_theory(itr) = p2_T;
    %----------------------------------------------------------------------    
end
%==========================================================================

%==============================PLOT RESULTS================================
figure
plot(parameter_range,success_measure(1,:),'g-*')
hold on
plot(parameter_range,success_measure(2,:),'k-o')
plot(parameter_range,acc_theory,'r')
plot(parameter_range,err_theory,'k')
title(['q=',num2str(q),' p=',num2str(p),' theta=',num2str(theta),' n=',num2str(n)])
legend('No. correct edges','No. false edges')
xlabel('T')
%==========================================================================