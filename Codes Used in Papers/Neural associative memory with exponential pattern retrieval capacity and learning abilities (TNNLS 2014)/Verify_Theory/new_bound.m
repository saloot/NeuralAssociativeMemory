%==========================================================================
%***********FUNCTION: new_bound(n,m,gamma,d,lambda,error_bits)*************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% n: The number of variable (left) nodes.
% m: The number of check (right) nodes.
% gamma: The update threshold for the bit-flipping algorithm.
% d: The degree distribution of variable nodes.
% lambda: The distribution of degrees for variable nodes.
% error_bits: number of initial errors
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe: The average pattern probability of error in the first iteration
% Pb: The average bit probability of error in the first iteration
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function calculates a bound on the probability of making 
% error when we have a given number of input errors. 
% The correction algorithm is based on the bit-flipping method introduced 
% in the ITW paper. The formulas based on which this code works are the 
% approximation of those in the progress report of 17-30 March 2012.
% However, as one expects, these formulas are too loose and usually results
% in trivial bounds.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================

function [Pe,Pb] = new_bound(n,m,gamma,d,lambda,error_bits)
%-----------------------------Initialization-------------------------------
Pe = [];
Pb = [];
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
alpha = 1-deg_ave/m;       
%--------------------------------------------------------------------------

for k = 1:length(error_bits)
    t = error_bits(k);                  % Number of input errors
    
    
        
    
    P1_tot = 0;                         
    P2_tot = 0;        
    
    %-----Find Average Error Probability Over the Degree Distribution------
    if (t == 0)
        Pe = [Pe,0];
        Pb = [Pb,0];
    else
        
    for j = 1:length(d) 
        dp = d(j);
        if (t>= 1+log(2)/log(m/(m-deg_ave)))
            P2 = ((1-alpha^(t-1))^dp) *2^(dp-1);
            P1 = (1-alpha^t)^(dp) *((exp(1)/gamma)^(gamma*dp))*(1-gamma)*dp;
        elseif (t == log(2)/log(m/(m-deg_ave)))
            P2 = (alpha^((t-1)*(dp/2)))*(1-alpha^(t-1))^(dp/2) *2^(dp-1);
            P1 = (1-alpha^t)^(dp) *((exp(1)/gamma)^(gamma*dp))*(1-gamma)*dp;
        else 
            P2 = (alpha^((t-1)*(dp/2)))*(1-alpha^(t-1))^(dp/2) *2^(dp-1);
            P1 = (alpha^(t*dp*(1-gamma)) ) *((1-alpha^t)^(dp*gamma)) *((exp(1)/gamma)^(gamma*dp))*(1-gamma)*dp;
            if(P1>1)
                111;
            end
        end
        P1_tot = P1_tot + lambda(j)*P1;                       
        P2_tot = P2_tot + lambda(j)*P2;                       
                                        
    end 
    P_tot = (P2*t/n)+(P1*(1-t/n));            
 
    Pex = (t*P2_tot/n)+(1-t/n)*P1_tot;    
    Pe = [Pe,1-(1-Pex)^n];
    Pb = [Pb,Pex];
    end
    
    %----------------------------------------------------------------------
end
plot(Pe,'r*')
% hold on
% plot(Pe_v2,'g*')
% Pe = Pe_v1;