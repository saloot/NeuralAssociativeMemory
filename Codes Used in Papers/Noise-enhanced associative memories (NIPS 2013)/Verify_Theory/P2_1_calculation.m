%==========================================================================
%**********************FUNCTION: P2_1_calculation**************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% gamma: The update threshold in the original bit-flipping recall algorithm 
% deg: The degree of the non-corrupted pattern neuron
% Pai2_1: the probability that a constrain neuron make a mistake a not return a +1/0/-1 message when it should
% P_neighbor: the probability that the considered pattern neuron has a common neighbor with the corrupted pattern neuron
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe: The upper bound on the probability that a non-corrupted pattern neuron mistakenly updates its state
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates a theoretical upper bound on the probability
% that a non-corrupted pattern neuron updates its state by mistake in a
% given cluster. The bound is explained in our ISIT 2013 paper.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function Pe = P2_1_calculation(gamma,pattern_neur_noise,deg,deg_noisy,Pai2_1,P_neighbor)
                                


Pe = 0;

for b = 0:min(deg,deg_noisy)                                               % b is the number of common neighbors ,d1)
    temp_ave = 0;
    for b1 = 0:b                                            % b1 is the number of non-zero messages received from the common neighbors
        temp = 0;
        for e = 0:b1                                        % e is the number of messages that share the same sign with the edges they are coming from
            temp = temp + nchoosek(b1,e) *p_1_2_j(2*e-b1,pattern_neur_noise,gamma,deg);
        end
        temp = temp/(2^b1);
    
        temp_ave = temp_ave + nchoosek(b,b1)*((1-Pai2_1)^b1)*(Pai2_1^(b-b1))*temp;
    end
    
    Pe = Pe + nchoosek(deg,b)*((P_neighbor)^b)*((1-P_neighbor)^(deg-b))*temp_ave;
end
