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


% function [x_out] = recall_step_clustered_real_value_v2(W,x_orig,varphi,y_min,y_max,max_y_threshold)

%============================INITIALIZATION================================
exit_flag = 0;        
itr = 0;                
[N_const,N] = size(W);
x = x_orig;
% W = W./(ones(N_const,1)*sqrt(sum(W.*W)));
alph0 = 0.95/norm(x)^2;
alph_itr= 0;
cost = [];
%==========================================================================
        


while (exit_flag == 0)
    itr = itr + 1;           
                                                                                                   
    %----------------------Update Constraint Nodes-------------------------
    y = W*x';                                      % Find which constraint nodes are violated.                                           
    
    if (mod(itr,10) == 1)
        alph_itr = alph_itr + 1;
        alph = 8*alph0/(8+alph_itr);
    end
%     c = c_temp.*(abs(c_temp)>max_y_threshold);
    %----------------------------------------------------------------------
    
    %-----------------------Check for Convergence--------------------------
    if ( (norm(y) < 0.0001)||(itr > 40000))
        exit_flag = 1;                                                  
        break;
    end
    %----------------------------------------------------------------------                                                                       

    %---------------------Update Pattern Nodes-----------------------------
    x = x - alph * (tanh(100*y'*W));
    cost = [cost,norm(x-projected_pattern)];
    norm(x-projected_pattern)
    %----------------------------------------------------------------------
                        
end
x_out = x;

