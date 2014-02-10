%==========================================================================
%************************FUNCTION: spatial_potential***********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% lambda: The pattern (variable) nodes degree distribution
% rho: The constraint/cluster nodes degree distribution
% e: The number of errors each constraint/cluster can correct (0 or 1)
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% th_const: The threshold for the constrained coupled system, i.e. the one in which the boundaries are fixed
% th_unconst: The threshold for the unconstrained coupled system, i.e. the one in which the boundaries are not fixed
%--------------------------------------------------------------------------

%---------------------------FUNCTION DESCRIPTION---------------------------
% This piece of code computes the error correction threshold of a spatially-coupled
% neural associative memory which is capable of correcting e errors (e = 0
% or 1). It is based on the formulas specified in the paper "A simple proof
% of threshold saturation for coupled scalar recursions" by Yedla et al. 
%
% The code first computes the potential function of the system for various
% initial error probabilities and then finding the maximum such probability
% for which he potential function is always positive (th_const). The function
% also calculates the same threshold for the unconstrained coupled system,
% i.e. the one in which the boundaries are not fixed (th_unconst). This correspond to
% the maximum initial probability for which the potential has positive
% derivative everywhere so that in each step of recursion, the error
% probability is reduced. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

% function [th_const,th_unconst] = spatial_potential(lambda,rho,e)
e = 1;
version = 2;

%==============================INTIALIZATION===============================
addpath(genpath('/home1/amir/cluster/Common_Library'));                          % Include the library of common functions in the search path
x = [0:.01:1];                                                          % The domain of the potential function
flag = zeros(1,2);                                                      % The flag is used to decide when thresholds are found
p_e_range = [0.00001:.001:1];                                           % The domain of initial channel bit error probability
%==========================================================================


%=================================MAIN LOOP================================
for itr = 2:length(p_e_range)
    p_e = p_e_range(itr);                                               % Set the initial error probability
    
    %-------------------Calculate the Potential Function-------------------
    if (e == 0)
        if (version == 2)
            g = 1-rho_poly(1-x,rho);
            GG = x+rho_poly_int(1-x,rho)-rho_poly_int(1,rho);
        else
            g = 1;
            GG = x;
        end
    elseif (e == 1)
        if (version == 2)
            g = 1-rho_poly(1-x,rho)-x.*rho_poly_der(1-x,rho);
            GG = x+rho_poly_int(1-x,rho)-rho_poly_int(1,rho)+x.*rho_poly(1-x,rho) +rho_poly_int(1-x,rho)./(1-x)-rho_poly_int(1,rho);
        else
            g = 1-rho_poly(1-x,rho);
            GG = x+rho_poly_int(1-x,rho)-rho_poly_int(1,rho);
        end
    else
        error ('This mode is not supported yet!');
    end

    U = x.*g-GG-p_e*lambda_poly_int(g,lambda);
    EE = x-p_e*lambda_poly(g,lambda);
    %----------------------------------------------------------------------

    %----------Check If The Constrained Threshold Has Been Found-----------
    if (min(U)<0)
        if (flag(1,1) ==0)
            th_const = p_e_range(itr-1);
            flag(1,1) = 1;
        end
    end
    %----------------------------------------------------------------------
    
    %---------Check If The Unconstrained Threshold Has Been Found----------
    if (min(EE)<0)
        if (flag(1,2) ==0)
            th_unconst = p_e_range(itr-1);
            flag(1,2) = 1;
        end
    end
    %----------------------------------------------------------------------
    
    %-------------Exit the Loop When Both Thresholds Are found-------------
    if (norm(flag) > 1)
        break;
    end
    %----------------------------------------------------------------------
end
%==========================================================================
111
% plot(x,U)