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
% The following function calculates the probability of making error when we
% have a given number of input errors. The correction algorithm is based on
% the bit-flipping method introduced in the ITW paper. The formulas based 
% on which this code works can be found in the progress report of 17-30
% March 2012.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function [Pe,Pe_bit] = Error_prob_ireg(n,m,gamma,d,lambda,error_bits)


%-----------------------------Initialization-------------------------------
Pe = [];
Pe_bit = [];
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
d_max = max(d);
if (norm(sum(lambda)-1)>1e-5)
    error('Invalid degree distribution');
end
%--------------------------------------------------------------------------

for k = 1:length(error_bits)
    t = error_bits(k);                  % Number of input errors
    
    
    
    %-----Find Average Error Probability Over the Degree Distribution------
    if (t>0)
        temp3 = 0;
        for jj = 1:length(d)
            if (lambda(jj) == 0)
                continue;
            end
            alpha = 1-d(jj)/m;            
            P2_tot = 0;                         
            P3_tot = 0;        
            for j = 1:length(d)        
                if (lambda(j) == 0)
                    continue;
                end
                P2 = 0;        
                P3 = 0;        
                dp = d(j);                       
                        
                for i = 0:ceil(gamma*dp)-1
                    P2 = P2 + nchoosek(dp,i) * ((1-alpha^t)^i) * (alpha^(t*(dp-i)));                            
                    if ((i <= floor(0.5*dp))&&(t>1))              
                        P3 = P3 + nchoosek(dp,i) * ((1-alpha^(t-1))^i) * (alpha^((t-1)*(dp-i)));                              
                    end            
                end
    
                P2_tot = P2_tot + lambda(j)*P2;                
                P3_tot = P3_tot + lambda(j)*P3;             
            end    
            if (t == 1)
                P3_tot = 1;
            end
            if (gamma ==1)
                111;
            end
            temp3 = temp3 + lambda(jj)*((t*P3_tot/n)+(1-t/n)*P2_tot);
        end
        Pex = 1-temp3;    
        Pe_bit = [Pe_bit,Pex];
        Pe = [Pe,1-(1-Pex)^n];
    else
        Pe = [Pe,0];
        Pe_bit = [Pe_bit,0];
    end
    %----------------------------------------------------------------------
end
% plot(Pe,'r*')