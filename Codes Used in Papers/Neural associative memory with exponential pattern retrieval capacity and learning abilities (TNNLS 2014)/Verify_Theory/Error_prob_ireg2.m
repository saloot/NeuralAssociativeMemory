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


function [Pe,q_e_c,q_c_c] = Error_prob_ireg2(n,m,gamma,d,lambda,rho,error_bits)

%-----------------------------Initialization-------------------------------
Pe = [];
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
alpha = 1-deg_ave/m;     

%--------------------------------------------------------------------------

for k = 1:length(error_bits)
    t = error_bits(k);                  % Number of input errors
    q_e_c = 0;
    q_c_c = 0;
    P2_tot = 0;                         
    P3_tot = 0;        
    P_tot = 1;
    %-----Find Average Error Probability Over the Degree Distribution------
    rho_pol_val = rho_poly(1-t/n,rho);
    pii = 1-rho_pol_val;    
    pii2 = pii;%1-(1-t/n)*rho_pol_val;    
    if (t > 0)
    for j = 1:length(d)        
                
        if (lambda(j) ==0)
            continue;
        else
            P2 = 0;        
            P3 = 0;        
            dp = d(j);                                           
        end
        if (dp > 0)
        for i = ceil(gamma*dp):dp
%             P2 = P2 + nchoosek(dp,i) * ((1-alpha^t)^i) * (alpha^(t*(dp-i)));                            
             P2 = P2 + nchoosek(dp,i) * (pii2^i) * ((1-pii2)^(dp-i));                            
        end        
            P2 = 1-P2;
        else
            P2 = 1;
        end
        
        p_error = pii;
%         p_error = (1-alpha^(t-1));
        if (t>1)
            for i = 0:floor(0.5*dp)
%                 P3 = P3 + nchoosek(dp,i) * ((1-alpha^(t-1))^i) * (alpha^((t-1)*(dp-i)));                             
                  P3 = P3 + nchoosek(dp,i) * (p_error^i) * ((1-p_error)^(dp-i));                             
            end            

%             P3  = P3 + .5*(1-P3);
        else
            P3 = 1;
        end
        if (dp == 0)
            P3 = 0;
        end
        P_tot = P_tot * (((P3*t/n)+(P2*(1-t/n)))^(n*lambda(j)));            
        q_c_c = q_c_c + lambda(j)*P2;                        
        q_e_c = q_e_c + lambda(j)*P3;
    end    
    else
        P_tot = 1;
    end
 
    % Pex = 1-(t*P3_tot/n)-(1-t/n)*P2_tot;    
%     x = 1-P_tot;
%     Pe = [Pe,1-(1-x)^n];
    Pe = [Pe,1-P_tot];
    %----------------------------------------------------------------------
end
% plot(Pe,'r*')
