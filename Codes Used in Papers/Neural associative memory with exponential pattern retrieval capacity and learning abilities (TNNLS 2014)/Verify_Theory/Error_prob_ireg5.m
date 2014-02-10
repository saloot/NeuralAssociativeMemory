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


function [Pe,Pb] = Error_prob_ireg5(n,m,gamma,d,lambda,error_bits,try_max)

%-----------------------------Initialization-------------------------------
Pe = [];
Pb = [];
S = 10;
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
%--------------------------------------------------------------------------

for k = 1:length(error_bits)
            
    %-----Find Average Error Probability Over the Degree Distribution------
    error_bits(k);
    x = error_bits(k)/n;
    e = zeros(1,S+1); 
    e(2) = x;
    e(1) = 1-sum(e(2:S+1));
    end_flag = 0;
    alpha = 1-deg_ave/m;            
    try_itr = 0;
    xx = [];
    while (end_flag == 0)
        
        try_itr = try_itr + 1;
        P_tot = 0;
                    
        for t = 0:n
%             if (try_itr==1)
%                 t = error_bits(k);
%             end
%             t = round(x*n);
            tt = max(t-1,0);
            P1 = 0;
            P2 = 0;
            for j = 1:length(d)                        
                dp = d(j);                                                       
                P1_x = 0;
                P2_x = 0;
                for i = ceil(gamma*dp):dp
                    P1_x = P1_x + nchoosek(dp,i) * ((1-alpha^t)^i) * (alpha^(t*(dp-i)));                            
                end        
            
                for i = ceil(dp*.5):dp
                    P2_x = P2_x + nchoosek(dp,i) * ((1-alpha^tt)^i) * (alpha^((tt*(dp-i))));                            
                end
                P1 = P1 + lambda(j)*P1_x;
                P2 = P2 + lambda(j)*P2_x;
            end
         
            %-----------------Calculate Error Evolution--------------------
            temp = e;
            for j = 3:S
                e(j) = temp(j+1)*(1-P2) + temp(j-1)*(P2);            
            end
            e(2) = temp(3)*(1-P2) + temp(1)*P1;
            e(S+1) = temp(S)*(P2)+temp(S+1)*(P2);
            e(1) = temp(1)*(1-P1)+temp(2)*(1-P2);
            e = min(e,ones(1,length(e)));
            e = max(e,zeros(1,length(e)));
%             if (try_itr==1)
%                 P_tot = 1-e(1);
                
%                 break;
%             else
            P_tot = P_tot + nchoosek(n,t)*(x^t)*((1-x)^(n-t))*(x+(e(1))*P1-(e(2))*(1-P2));                                    
%             P_tot = P_tot + nchoosek(n,t)*(x^t)*((1-x)^(n-t))*(1-e(1));            
%             P_tot = x+(e(1))*P1-e(2)*(1-P2);
%             end
        end
        x = max(P_tot,0);
        x = min(x,1);
%         x
%         e
%         P1
%         P2
%         x = P_tot
%         x = 1-e(1);
        xx = [xx,x];
        
        if ((try_itr > try_max) || (x == 1) || (x<=1e-5))
            end_flag = 1;
        end
    end
    % Pex = 1-(t*P3_tot/n)-(1-t/n)*P2_tot;    
    Pe = [Pe,1-(1-x)^n];
    Pb = [Pb,x];
    %----------------------------------------------------------------------
end
% plot(Pe,'r*')
