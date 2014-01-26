%==========================================================================
%*********************FUNCTION: P2_2_calculation***************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% gamma: The update threshold in the original bit-flipping recall algorithm 
% deg: The degree of the corrupted pattern neuron
% Pai2_1: the probability that a constrain neuron make a mistake a not return a +1/0/-1 message when it should
% P_common_neighbor: The probability that the two noisy neurons share a common neighbor
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe1_no_update: The upper bound on the probability that the first corrupted neurons is not updated at all
% Pe1_wrong_dir: The upper bound on the probability that the first corrupted neurons is updated in the wrong direction
% Pe2_no_update: The upper bound on the probability that the second corrupted neurons is not updated at all
% Pe2_wrong_dir: The upper bound on the probability that the second corrupted neurons is updated in the wrong direction
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates a theoretical upper bound on the probability
% that a corrupted pattern neuron does not update its state by mistake in a
% given cluster WHEN there are TWO external erros in the cluster AND the
% CONSTRAINT neurons are NOISELESS. The reason we have chosen constraint 
% neurons to be noiseless is that we would like to obtain tighter 
% theoretical bounds for BER.

% NOTE: This version differs from the previous one in that it calculates
% all expectations NUMERICALLY so it is more time consuming but, hopefuly,
% more accurate. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function [Pe1_no_update,Pe1_wrong_dir,Pe2_no_update,Pe2_wrong_dir] = P2_2_calculation_e2_v2(gamma,pattern_neur_noise,deg1,deg2,Pai2_1,d_bar)

if (pattern_neur_noise > 0 )          
    temp1_wrong_dir = 0;    
    temp2_wrong_dir = 0;    
    
    temp1_no_update = 0;    
    temp2_no_update = 0;  
    for e = 0:d_bar            
        x1 = deg1+2*e-2*d_bar;                    
        x2 = deg2-2*e;                    
        
        x1_wrong_dir = max(0,pattern_neur_noise-(gamma+x1/deg1))/(2*pattern_neur_noise);                    
        x2_wrong_dir = max(0,pattern_neur_noise-(gamma+x2/deg2))/(2*pattern_neur_noise);            
        
        temp1 = min(gamma-x1/deg1,pattern_neur_noise);    
        temp2 = min(gamma+x1/deg1,pattern_neur_noise);            
        x1_no_update = max(0,temp1+temp2)/(2*pattern_neur_noise);                    
        
        temp1 = min(gamma-x2/deg2,pattern_neur_noise);    
        temp2 = min(gamma+x2/deg2,pattern_neur_noise);    
        x2_no_update = max(0,temp1+temp2)/(2*pattern_neur_noise);                    
        
        temp1_no_update = temp1_no_update + nchoosek(d_bar,e) * min (x1_no_update,1);        
        temp2_no_update = temp2_no_update + nchoosek(d_bar,e) * min (x2_no_update,1);        
        
        temp1_wrong_dir = temp1_wrong_dir + nchoosek(d_bar,e) * min (x1_wrong_dir,1);        
        temp2_wrong_dir = temp2_wrong_dir + nchoosek(d_bar,e) * min (x2_wrong_dir,1);        
    end
    
    Pe1_wrong_dir = temp1_wrong_dir*(.5^d_bar);
    Pe2_wrong_dir = temp2_wrong_dir*(.5^d_bar);
                    
    Pe1_no_update = temp1_no_update*(.5^d_bar);
    Pe2_no_update = temp2_no_update*(.5^d_bar);
    
elseif (pattern_neur_noise == 0)    
    temp1 = 0;    
    temp2 = 0;    
    for e = 0:d_bar                        
        x1 = deg1+2*e-2*d_bar;            
        x2 = deg2-2*e;                    
        if (x1<gamma*deg1)                
            x1 = 1;            
        else            
            x1 = 0;
        end        
        if (x2<gamma*deg2)                
            x2 = 1;            
        else            
            x2 = 0;
        end        
        x1_wrong_dir = max(0,pattern_neur_noise-(gamma+x1/deg1))/(2*pattern_neur_noise);                    
        x2_wrong_dir = max(0,pattern_neur_noise-(gamma+x2/deg2))/(2*pattern_neur_noise);            
        
        temp1 = min(gamma-x1,pattern_neur_noise);    
        temp2 = min(gamma+x1,pattern_neur_noise);            
        x1_no_update = max(0,temp1+temp2)/(2*pattern_neur_noise);                    
        
        temp1 = min(gamma-x2,pattern_neur_noise);    
        temp2 = min(gamma+x2,pattern_neur_noise);    
        x2_no_update = max(0,temp1+temp2)/(2*pattern_neur_noise);                    
        
        temp1_no_update = temp1_no_update + nchoosek(d_bar,e) * min (x1_no_update,1);        
        temp2_no_update = temp2_no_update + nchoosek(d_bar,e) * min (x2_no_update,1);        
        
        temp1_wrong_dir = temp1_wrong_dir + nchoosek(d_bar,e) * min (x1_wrong_dir,1);        
        temp2_wrong_dir = temp2_wrong_dir + nchoosek(d_bar,e) * min (x2_wrong_dir,1); 
    end    
    Pe1_wrong_dir = temp1_wrong_dir*(.5^d_bar);
    Pe2_wrong_dir = temp2_wrong_dir*(.5^d_bar);
                    
    Pe1_no_update = temp1_no_update*(.5^d_bar);
    Pe2_no_update = temp2_no_update*(.5^d_bar);
    
else
    error('Invalid pattern_neur_noise!pattern_neur_noise should be positive.');
end

