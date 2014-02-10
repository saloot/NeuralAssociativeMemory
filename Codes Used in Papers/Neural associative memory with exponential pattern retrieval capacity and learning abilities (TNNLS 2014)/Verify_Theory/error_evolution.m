%==========================================================================
%**********************FUNCTION: error_evolution***************************
%==========================================================================

%--------------------------------INPUTS------------------------------------

% lambda: the distribution of column degrees. lambda(i) represents the fraction of pattern nodes with degree i.
% deg_col: the degree of pattern nodes, from 1 to N-K
% rho: the distribution of row degrees. rho(i) represents the fraction of constraint nodes with degree i.
% deg_row: the degree of constraint nodes, from 1 to N
% x: The initial probability of a neuron being NOT noisy.
% S: The maximum firing rate (state) of a neuron
% itr_max: The maximum number of iterations spent on tracking the probability of error
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% success_flag: Indicates if the probability of error has converged to 0
% zz: The vector containing the probability of error in each iteration
% e: The probability of a neuron being in the state of having a noise with amplitude (absolute value) varying from 0 to S.
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This piece of code tracks the evolution of the probability of error in a 
% single-cluster neural network with non-binary noise. The code tracks the 
% probability of each pattern neuron having a noise value from 1 to a 
% maximum S. If in the end of simulation all these probabilities converged 
% to zero, a success is declared. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function [success_flag,zz,e] = error_evolution(lambda,deg_column,rho,deg_row,x,S,itr_max)

%=============================INITITIALIZATION=============================
z = 1-x;                                    % This is the initial probability of error
success_flag = 0;                           % Indicates if the probability of error has converged to 0
zz=[z];                                     % The vector containing the probability of error in each iteration

e = zeros(1,S+1);                           % The probability of a neuron being in the state of having a noise with amplitude (absolute value) varying from 0 to S 

%--------------------Determine the Distribution of Noise-------------------
e(2) = 1*z;
e(1) = 1-sum(e(2:S+1));
% e(2) = 7*z;
% e(3) = .2*z;
%--------------------------------------------------------------------------

%==========================================================================
for itr=1:itr_max
   
    %---------Calculate the Probability of Error in This Iteration---------
    P2_tot = zeros(1,length(x));                
    
    alpha = rho_poly(x,rho);    
    
    for j = 1:length(deg_column)    
    
        P2 = zeros(1,length(x));                                    
        dp = deg_column(j);                            
        for i = ceil(0.5*dp):dp                    
            P2 = P2 + nchoosek(dp,i) * ((1+alpha).^i) .* ((1-alpha).^(dp-i)) ;                      
        end        
        P2_tot = P2_tot + lambda(j)*P2*(.5)^dp;        
    end    
    g = P2_tot;            
    f_rho = rho_poly(x,rho);    
    f1 = lambda_poly(1-x.*f_rho,lambda);
    %----------------------------------------------------------------------
    
    %--------------Calculate the Evolution of Noise Values-----------------    
    temp = e;        
    for j = 3:S
        e(j) = temp(j+1)*g + temp(j-1)*(1-g);            
    end
    e(2) = temp(3)*g + temp(1)*f1;
    e(S+1) = temp(S)*(1-g)+temp(S+1)*(1-g);
    e(1) = temp(1)*(1-f1)+temp(2)*g;
    e = min(e,ones(1,length(e)));
    e = max(e,zeros(1,length(e)));
    %     e = e/sum(e);
    z = 1-e(1);
    %----------------------------------------------------------------------
    
    %------------------Update the Probability of Error---------------------
    % z = e(1)*f1 - e(2)*g+z;
    z = max(z,0);
    z = min(z,1);
    
    x = 1-z;        
    zz= [zz,z];
    %----------------------------------------------------------------------
    
    %------------------------Check for Convergence-------------------------
    if (z <= 1e-6)        
        success_flag = 1;           
        break;
    end
    %----------------------------------------------------------------------

end