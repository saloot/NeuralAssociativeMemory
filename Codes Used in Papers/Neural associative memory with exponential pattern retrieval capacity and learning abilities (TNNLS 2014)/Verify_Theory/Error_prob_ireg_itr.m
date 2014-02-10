%==========================================================================
%*********FUNCTION: Error_prob_ireg(n,m,gamma,d,lambda,error_bits)*********
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
% The following function is part of the code to assess the conditions under
% which the number of erroneous nodes decrease in each iteration of the 
% bit-flipping method introduced in the ITW paper. The idea is to find two
% functions (hz and gz) which are functions of the number of initial input
% errors (z). Then if gz < hz for all values of z, the norm-1 of error
% decreases in each iteration. The formulas based 
% on which this code works can be found in the progress reports of 17-30 
% March and 2-15 April 2012.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================

function [fz,gz,hz] = Error_prob_ireg_itr(n,m,gamma,d,lambda)

%----------------------------Initialization--------------------------------
deg_ave = sum(d.*lambda);
fz = [];
gz = [];
hz = [];
%--------------------------------------------------------------------------

%-----------------------Calculate hz and gz for Each z---------------------
for z=0:n    
    if (z ==0)
        fz = [fz,0];
        gz = [gz,0];
        hz = [hz,0];
    else
       
        P2_tot = 0;
        P1_tot = 0;
        for j = 1:length(d)
            P2 = 0;
            P1 = 0;
            dp = d(j);
            alpha = 1-deg_ave/m;
    
            for i = ceil(gamma*dp):dp
                P1 = P1 + nchoosek(dp,i) * ((1-alpha^z)^i) * (alpha^(z*(dp-i)));                                
            end
            if (z>1)
                for i = ceil(0.5*dp):dp
                    P2 = P2 + nchoosek(dp,i) * ((1-alpha^(z-1))^i) * (alpha^((z-1)*(dp-i)));                
                end               
            else
                P2 = 0;
            end
            P1_tot = P1_tot + lambda(j)*P1;
            P2_tot = P2_tot + lambda(j)*P2;
        end
        gz = [gz,P2_tot];        
        fz = [fz,P1_tot]; 
        hz = [hz,(n-z)*P1_tot+z*(1*P2_tot-1)];
    end
    
end
%--------------------------------------------------------------------------


%----------------------------Plot Results----------------------------------
plot(hz,'b')
hold on
% plot(gz,'r')
% legend('hz','gz')
title('We are looking for values of  for which gz<hz');
%--------------------------------------------------------------------------