%==========================================================================
%*******************FUNCTION: recall_step_clustered************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% W: The connectivity matrix of the neural graph
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


function [x_out] = recall_step_clustered_real_value(W,x,varphi,y_min,y_max,max_y_threshold)

%============================INITIALIZATION================================
exit_flag = 0;        
itr = 0;                
[N_const,N] = size(W);
% W = W./(ones(N_const,1)*sqrt(sum(W.*W)));
%==========================================================================
        


while (exit_flag == 0)
    itr = itr + 1;           
                                                                                                   
    %----------------------Update Constraint Nodes-------------------------
    c_temp = W*x' ;                                      % Find which constraint nodes are violated.                                           
    if (itr > 1)
        c_old = c;
    end
    c = c_temp.*(abs(c_temp)>max_y_threshold);
    %----------------------------------------------------------------------
    
    %-----------------------Check for Convergence--------------------------
    if ( (norm(c) < 0.0001)||(itr > 4000))
        exit_flag = 1;                                                  
        break;
    end
    
%     if (itr > 1)
%         if (norm(c)-norm(c_old)> .0001)
%             break;
%         end
%     end
    %----------------------------------------------------------------------                                                                       

    %---------------------Update Pattern Nodes-----------------------------
    x_temp = c'*W./sqrt(sum(W.*W));                                      % x_temp is the raw feedback received by each pattern node.                                                    
    
    [val,ind] = max(abs(x_temp));    
    x = x - x_temp(ind);                
%     x = x - x_temp.*(abs(x_temp)<varphi);
%     x = min(max(x,y_min),y_max);
    %----------------------------------------------------------------------
                        
end
x_out = x;

