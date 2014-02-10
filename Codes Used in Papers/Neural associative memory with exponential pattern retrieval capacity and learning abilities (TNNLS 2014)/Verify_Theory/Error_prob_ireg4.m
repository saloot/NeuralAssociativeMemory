%==========================================================================
%*********FUNCTION: Error_prob_ireg2(n,m,gamma,d,lambda,error_bits)*********
%==========================================================================

%--------------------------------INPUTS------------------------------------
% n: The number of variable (left) nodes.
% m: The number of check (right) nodes.
% gamma: The update threshold for the bit-flipping algorithm.
% d: The degree distribution of variable nodes.
% lambda: The distribution of degrees for variable nodes.
% error_bits: number of initial errors
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe: The average probability of error in the first iteration
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function calculates the probability of making error when we
% have a given number of input errors. The correction algorithm is based on
% the bit-flipping method introduced in the ITW paper. The formulas based 
% on which this code works can be found in the progress report of 17-30
% March 2012.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function Pe = Error_prob_ireg4(n,m,gamma,d,lambda,rho,error_bits,try_max)

%-----------------------------Initialization-------------------------------
Pe = [];
Pb = [];
S = 10;
rho = rho/sum(rho);
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
alpha = 1-deg_ave/m;            
%--------------------------------------------------------------------------

for k = 1:length(error_bits)
    x = error_bits(k)/n;                  % Number of input errors    
    e = zeros(1,S+1); 
    e(2) = x;
    e(1) = 1-sum(e(2:S+1));
    end_flag = 0;    
    try_itr = 0;
    x_itr = [];
    while (end_flag == 0)        
        try_itr = try_itr + 1;
        q_c_w = 0;
        q_e_c = 0;
        rho_pol_val = rho_poly(1-x,rho);
        pii = 1-(1-x)*rho_pol_val;    
        tt = max(0,round(x*n)-1);
        p_cor = 1*rho_pol_val + .5*(1-rho_pol_val);
        if (p_cor>1)
            111;
        end
        %----Find Average Error Probability Over the Degree Distribution---
       
        for j = 1:length(d)        
            p_temp = 0;                
            dp = d(j);                                           
            for i = ceil(gamma*dp):dp
                p_temp = p_temp + nchoosek(dp,i)*((pii)^i)*((1-pii)^(dp-i));
            end        
            q_c_w = q_c_w + lambda(j)*p_temp;            
            if (round(x*n)>1)
                p_temp = 0;
                for i = ceil(0.5*dp):dp
                    p_temp = p_temp + nchoosek(dp,i)*((p_cor)^i)*((1-p_cor)^(dp-i));
%                     p_temp = p_temp + nchoosek(dp,i) * ((alpha^tt)^i) * ((1-alpha)^((tt*(dp-i)))); 
                end
                q_e_c = q_e_c + lambda(j)*p_temp;            
            else
                q_e_c = 1;
            end
        end    
        
        %------------------Density Evolution for---------------------------
        P2 = 1-q_e_c;
        P1 = q_c_w;
        
        temp = e;
        for j = 3:S                
            e(j) = temp(j+1)*(1-P2) + temp(j-1)*(P2);                        
        end        
        e(2) = temp(3)*(1-P2) + temp(1)*P1;            
        e(S+1) = temp(S)*(P2)+temp(S+1)*(P2);            
        e(1) = temp(1)*(1-P1)+temp(2)*(1-P2);            
        e = min(e,ones(1,length(e)));            
        e = max(e,zeros(1,length(e)));
        %------------------------------------------------------------------
        
        %--------------------Update Error Probability----------------------
%         x = x + (1-x)*q_c_w -x*q_e_c;
%         x = x + e(1)*q_c_w -e(2)*q_e_c;
        x = 1-e(1);
        x_itr = [x_itr,x];
        %------------------------------------------------------------------        
        
        %--------------------Check for Convergence-------------------------
        if ((try_itr > try_max) || (x == 1) || (x<=1e-5))
            end_flag = 1;
        end
        %------------------------------------------------------------------        
    end
    Pb = [Pb,x];
    % Pex = 1-(x*P3_tot/n)-(1-x/n)*P2_tot;    
    Pe = [Pe,1-(1-x)^n];
    %----------------------------------------------------------------------
end
111;
% plot(Pe,'r*')
