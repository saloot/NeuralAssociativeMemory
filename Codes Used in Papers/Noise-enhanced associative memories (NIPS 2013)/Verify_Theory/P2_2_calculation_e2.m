%==========================================================================
%*********************FUNCTION: P2_2_calculation***************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% gamma: The update threshold in the original bit-flipping recall algorithm 
% deg: The degree of the corrupted pattern neuron
% Pai2_1: The probability that a constrain neuron make a mistake a not return a +1/0/-1 message when it should
% P_common_neighbor: The probability that the two noisy neurons share a common neighbor
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe: The upper bound on the probability that a corrupted pattern neuron mistakenly does not update its state
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates a theoretical upper bound on the probability
% that a corrupted pattern neuron does not update its state by mistake in a
% given cluster WHEN there are TWO external erros in the cluster AND the
% CONSTRAINT neurons are NOISELESS.

% NOTE: The reason we have chosen constraint neurons to be noiseless is
% that we would like to obtain tighter theoretical bounds for BER.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function Pe = P2_2_calculation_e2(gamma,pattern_neur_noise,deg,Pai2_1,P_common_neighbor)

if (pattern_neur_noise > 0 )
    Pe = 0;
    for d_bar = 0:floor(1*deg)            
        temp = 0;
        for e = 0:d_bar
            x = deg+2*e-2*d_bar;
            x = max(0,gamma+pattern_neur_noise-x/deg)/(2*pattern_neur_noise);
            temp = temp + nchoosek(d_bar,e) * min (x,1);
        end
        P_neighbor = ((P_common_neighbor)^d_bar)* ((1-P_common_neighbor)^(deg-d_bar));
        Pe = Pe + nchoosek(deg,d_bar) * P_neighbor *(.5^d_bar) * temp;
    end
    
    
elseif (pattern_neur_noise == 0)
    Pe = 0;
    for d_bar = 0:deg
        temp = 0;
        for e = 0:d_bar            
            x = deg+2*e-2*d_bar;
            if (x<gamma*deg)
                x = 1;
            else
                x = 0;
            end            
            temp = temp + nchoosek(d_bar,e) * x;
        end
        P_neighbor = ((P_common_neighbor)^d_bar)* ((1-P_common_neighbor)^(deg-d_bar));
        Pe = Pe + nchoosek(deg,d_bar) * P_neighbor *(.5^d_bar) * temp;
    end
else
    error('Invalid pattern_neur_noise!pattern_neur_noise should be positive.');
end
