%==========================================================================
%**************************FUNCTION: rho_poly******************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% x: The point we would like to evaluate the function
% rho: The degree distribution coefficients
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% f: The value of the function
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates the value of the pecified degree distibution
% polynomial.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function f = rho_poly2(x,rho)
f = zeros(1,length(x));
for i = 1:length(rho)
    f = f+ rho(i)*x.^(i);                                 % Note the power
end