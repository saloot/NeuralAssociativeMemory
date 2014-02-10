function P_error_av = error_theory(N,K,gamma_BFO,processed_error_bits_BFO,alpha0,beta0,theta0,index_max,method)

P_error_av = zeros(1,length(processed_error_bits_BFO));
index_itr = 0;
for index_in = 1:index_max
    [deg_row,deg_column,lambda,rho] = deg_dribution_ind(N,K,alpha0,beta0,theta0,index_in);    
    if ((method==4)||(method == 7))
        eval(['[Pe,P_tot_itr] = Error_prob_ireg',num2str(method),'(N,N-K,gamma_BFO,deg_column,lambda,rho,processed_error_bits_BFO,10);']);
    else
        eval(['Pe = Error_prob_ireg',num2str(method),'(N,N-K,gamma_BFO,deg_column,lambda,rho,processed_error_bits_BFO);']);
    end
    
    if (sum(Pe)<length(Pe)-1)
        P_error_av = P_error_av + Pe;
        index_itr = index_itr+1;
    end
end
P_error_av = P_error_av/index_itr;