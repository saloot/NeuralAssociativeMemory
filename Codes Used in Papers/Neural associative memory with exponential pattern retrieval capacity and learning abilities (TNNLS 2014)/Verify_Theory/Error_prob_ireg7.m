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


function [Pe,P_tot_itr] = Error_prob_ireg7(n,m,gamma,d,lambda,rho,error_bits,try_max)

%-----------------------------Initialization-------------------------------
Pe = [];
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
S = 10;
t_itr = [];
P_tot_itr = [];
alpha = 1-deg_ave/m;            
%--------------------------------------------------------------------------

for k = 1:length(error_bits)
    P_tot_itr = [];
    x = error_bits(k)/n;                  % Number of input errors
    if (x == 0)
        P_total = 1;
        Pe = [Pe,1-P_total];
        continue;
    end
    e = zeros(1,S+1); 
    e(2) = x;
    e(1) = 1-sum(e(2:S+1));
    end_flag = 0;    
    try_itr = 0;
    x_itr = [];
    P_tot = 1;
    q_e_c = 1-x;
    q_c_w = x;
    q_c_c = 1;
    q_c_e = 1;
    multiplier = 1;
    P_total = 0;    
    
    while (end_flag == 0)        
        try_itr = try_itr + 1;
        
        
        t = ceil(x*n);                
        [P_tot,q_e_c,q_c_c] = Error_prob_ireg2(n,m,gamma,d,lambda,rho,t);
        P_tot = 1-P_tot;
                        
        %------------------Density Evolution for---------------------------
%         q_c_c = 1-q_c_e;                
        temp = e;
        for j = 3:S                
            e(j) = temp(j+1)*(q_e_c) + temp(j-1)*(1-q_e_c);                        
        end        
        e(2) = temp(3)*(q_e_c) + temp(1)*(1-q_c_c);            
        e(S+1) = temp(S)*(1-q_e_c)+temp(S+1)*(1-q_e_c);            
        e(1) = temp(1)*(q_c_c)+temp(2)*(q_e_c);            
        e = min(e,ones(1,length(e)));            
        e = max(e,zeros(1,length(e)));
        %------------------------------------------------------------------
        
        
        P_total = P_total + P_tot *multiplier;
        multiplier = multiplier*(1-P_tot);
        x = 1-e(1);
%         x = x - e(2)*q_e_c+(1-x)*(1-q_e_c);
        t_itr = [t_itr,x];
        
        P_tot_itr = [P_tot_itr,1-P_tot];            % Probability of SUCCESS in each round
        
        %--------------------Check for Convergence-------------------------
        if ((try_itr > try_max) || (x == 1) || (x<=1e-5))
            end_flag = 1;
        end
        %------------------------------------------------------------------        
        
    end
    
    % Pex = 1-(t*q_e_c/n)-(1-t/n)*q_c_c;    
%     Pe = [Pe,1-P_total];
    
    Pe = [Pe,1-(1-x)^n];
    %----------------------------------------------------------------------
end
% plot(Pe,'r*')
111
