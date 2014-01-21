%==========================================================================
%*******************FUNCTION: recall_step_clustered************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% W_total: The connectivity matrix of the neural graph
% x: The corrupted pattern
% N: The number of pattern nodes in the graph
% gamma_BFO: The threshold that determines the margin for majority voting in the original bit-flipping algoroithm 
% gamma_BFS: The threshold that determines the margin for majority voting in the simiplified bit-flipping algoroithm 
% max_noise_amp: The maximum amplitude of the generated noise. It should be an integer.
% err_bits: The number of initial erroneous nodes 
% y_min: The minimum value that a pattern node can have
% y_max: The minimum value that a pattern node can have
% recall_algorithm_option: Selects between bit-flipping or winner-take-all algorithm in the recall phase
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% x_out: The output pattern of the recall process
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function gets the specifying parameters of a bipartite
% neural network as well as a corrupted patterns and performs a recall step
% of our proposed algorithm in the NIPS 2012 paper. The step is defined as 
% the process of iteratively correcting the error until either there is no
% errors left or a maximum number of iterations has been executed.

% There are three options for the recall algorithm: bit flipping (original 
% and simplified) as well as the winner-take-all approach. The selection 
% between the two is determined by the recall_algorithm_option flag.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function [x_out,updated_flag] = recall_step_clustered_v2(W_total,x,varphi,psi,y_min,y_max,Q)

%============================INITIALIZATION================================
exit_flag = 0;        
itr = 0;                
[N_const,N] = size(W_total);
updated_flag= 1;
W_total_deg_dist_0 = sum(abs(W_total));
W_total_deg_dist = sum(abs(sign(W_total)));
%==========================================================================
        


while (exit_flag == 0)
    updated_flag = 1;
    itr = itr + 1;           
    
    %------------Calculate the Feedback and Decision Parameters------------                
    cost = (W_total*x') + .00 * (rand(N_const,1)-.5);                
    c = (abs(cost) > psi).* (sign(cost));                
    feedback = (W_total)'*c;                            
%     feedback = sign(W_total)'*c;                            
    normalized_feedback = .01*(rand(N,1)-.5) + feedback./(W_total_deg_dist_0');                
    normalized_feedback = normalized_feedback.*(W_total_deg_dist'>2);
%     normalized_feedback = feedback./(W_total_deg_dist');                
    %----------------------------------------------------------------------
                    
    %----------------------Update the Dataset------------------------------
    [val,ind] = max(abs(normalized_feedback));                
    if ( (val > varphi) ) ;% && (sum(abs(normalized_feedback)>varphi) == 1) )        
            x(ind) =  min(max(x(ind)-sign(normalized_feedback(ind)),y_min),y_max);                            
        updated_flag = 1;                
    else
        updated_flag = 0;
        break;  
        
    end

%     updated_flag = sign(temp_var).*(abs(temp_var) > varphi)/Q; 
%     x = x - updated_flag;
%     x =  min(max(x-sign(normalized_feedback').*(abs(normalized_feedback')>varphi),y_min),y_max);
    
    %----------------------------------------------------------------------
                                       
    %-----------------------Check for Convergence--------------------------
    if ( itr > 10)
        break;
    end    
    %----------------------------------------------------------------------                                                                       
               
    
end

%--------------------------Generate Output---------------------------------
x_out = x;
%--------------------------------------------------------------------------