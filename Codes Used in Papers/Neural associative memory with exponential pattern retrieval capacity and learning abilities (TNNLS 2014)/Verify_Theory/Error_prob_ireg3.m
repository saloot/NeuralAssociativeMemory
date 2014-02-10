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


function Pe = Error_prob_ireg3(n,m,gamma,d,lambda,rho,error_bits)

%-----------------------------Initialization-------------------------------
Pe = [];
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
%--------------------------------------------------------------------------

for k = 1:length(error_bits)
    t = error_bits(k)/n;                  % Number of input errors
    q_c_w = lambda_poly(1-rho_poly2(1-t,rho),lambda);
    q_e_c = 0;
    %-----Find Average Error Probability Over the Degree Distribution------
    for j = 1:length(d)        
        p_temp = 0;                
        dp = d(j);                                           
        for i = ceil(gamma*dp):dp
            p_temp = p_temp + ((1+rho_poly(1-t,rho))^i) *((1-rho_poly(1-t,rho))^(dp-i));
        end        
        q_e_c = q_e_c + (0.5)^dp *lambda(j)*p_temp;            
    end    
    Pb = t + (1-t)*q_c_w -t*q_e_c;
    % Pex = 1-(t*P3_tot/n)-(1-t/n)*P2_tot;    
    Pe = [Pe,1-(1-Pb)^n];
    %----------------------------------------------------------------------
end
% plot(Pe,'r*')
