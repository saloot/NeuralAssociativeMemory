%==========================================================================
%*********************FUNCTION: P2_2_calculation***************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% gamma: The update threshold in the original bit-flipping recall algorithm 
% deg: The degree of the corrupted pattern neuron
% Pai2_1: the probability that a constrain neuron make a mistake a not return a +1/0/-1 message when it should
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe: The upper bound on the probability that a corrupted pattern neuron mistakenly does not update its state
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates a theoretical upper bound on the probability
% that a corrupted pattern neuron does not update its state by mistake in a
% given cluster. The bound is explained in our ISIT 2013 paper.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function Pe = P2_2_calculation(gamma,pattern_neur_noise,deg,Pai2_1)

if (pattern_neur_noise > 0 )
    P2_2 = 0;
    for j = max(0,ceil(1-pattern_neur_noise-gamma)*deg):min(floor(1+pattern_neur_noise-gamma)*deg,deg)
        P2_2 = P2_2 + nchoosek(deg,j) * ((Pai2_1)^j) * ((1-Pai2_1)^(deg-j)) * (pattern_neur_noise+gamma-( (deg-j)/deg) )/(2*pattern_neur_noise);
    end
    for j = ceil(1+pattern_neur_noise-gamma)*deg:deg
        P2_2 = P2_2 + nchoosek(deg,j) * ((Pai2_1)^j) * ((1-Pai2_1)^(deg-j));
    end
elseif (pattern_neur_noise == 0)
    P2_2 = 0;
    for j = ceil(1-gamma)*deg:deg
        P2_2 = P2_2 + nchoosek(deg,j) * ((Pai2_1)^j) * ((1-Pai2_1)^(deg-j));
    end
else
    error('Invalid pattern_neur_noise! "pattern_neur_noise" should be positive.');
end

Pe = P2_2;