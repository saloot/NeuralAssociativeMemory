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


function [x_out] = recall_step_clustered(W,x,gamma_BFO,gamma_BFS,y_min,y_max,recall_algorithm_option,max_y_threshold)

%============================INITIALIZATION================================
exit_flag = 0;        
itr = 0;                
x_WTA = x;
x_BFO = x;
x_BFS = x;
[N_const,N] = size(W);
%==========================================================================
        


while (exit_flag == 0)
    itr = itr + 1;           
                                                                                                   
    %----------------------Update Constraint Nodes-------------------------
    c_temp_WTA = W*x_WTA' ;                                      % Find which constraint nodes are violated.                                       
    c_temp_BFO = W*x_BFO' ;                                      % Find which constraint nodes are violated.                                       
    c_temp_BFS = W*x_BFS' ;                                      % Find which constraint nodes are violated.                                       
    
    for iii = 1:N_const                
        if (abs(c_temp_WTA(iii))<max_y_threshold)                        % Take care of numerical issues                        
            c_temp_WTA(iii) = 0;                                
        end        
        
        if (abs(c_temp_BFO(iii))<max_y_threshold)                        % Take care of numerical issues                        
            c_temp_BFO(iii) = 0;                                
        end        
        
        if (abs(c_temp_BFS(iii))<max_y_threshold)                        % Take care of numerical issues                        
            c_temp_BFS(iii) = 0;                                
        end        
    end
    
    c_WTA = sign(c_temp_WTA);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    c_BFO = sign(c_temp_BFO);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    c_BFS = sign(c_temp_BFS);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    
    %----------------------------------------------------------------------
    
    %-----------------------Check for Convergence--------------------------
    if ( ((norm(c_WTA) < 0.0001)&&(norm(c_BFO) < 0.0001)&&(norm(c_BFS) < 0.0001) )||(itr > 200))                                                            
        exit_flag = 1;                                                  
        break;
    end    
    %----------------------------------------------------------------------                                                                       

    %---------------------Update Pattern Nodes-----------------------------
    x_temp_WTA = c_WTA'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    x_temp_BFO = c_BFO'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    x_temp_BFS = c_BFS'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    
    
    %-----------------------Winner-Take-All--------------------------------
    temp_var = abs(c_WTA')*abs(sign(W))./sum(abs(sign(W)));                % Normalize the feedback by the norm-1 of the outgoing edges    
    [~,inde] = max(temp_var);                           % Find the node(s) with the maximum input sum.                                                                                                                    
    x_WTA(inde)=x_WTA(inde)-sign(x_temp_WTA(inde));             % Update the code nodes.                                        
    %----------------------------------------------------------------------
                                       
    %-------------------Bit Flipping Original------------------------------                
    x_temp2 = abs(c_BFO')*abs(sign(W))./sum(abs(sign(W)));               % x_temp2 is the raw feedback received by each pattern node.                                                       
        
    tempp = zeros(1,N);                                        
    for jjj = 1:N                                    
        if ( x_temp2(jjj)>= gamma_BFO)          % Update all the nodes that receive a lot of feedback from their neighbors.                                                                              
            tempp(jjj) = sign(x_temp_BFO(jjj));                                                        
        else                        
            tempp(jjj) = 0;                            
        end        
    end    
    
    x_BFO = x_BFO - tempp;                
    %----------------------------------------------------------------------
                       
                
    %-------------------Bit Flipping Simplified----------------------------
    x_temp2 = abs(c_BFS'*W)./sum(abs(W));               % x_temp2 is the raw feedback received by each pattern node.                                                       
        
    tempp = zeros(1,N);                                        
    for jjj = 1:N                                    
        if ( x_temp2(jjj)>= gamma_BFS)          % Update all the nodes that receive a lot of feedback from their neighbors.                                                                              
            tempp(jjj) = sign(x_temp_BFS(jjj));                                                        
        else                        
            tempp(jjj) = 0;                            
        end        
    end    
    
    x_BFS = x_BFS - tempp;           
    %----------------------------------------------------------------------
        
        
    x_WTA = max(x_WTA,y_min*ones(1,N));         % Wacth for saturations            
    x_WTA = min(x_WTA,y_max*ones(1,N));         % Wacth for saturations                    
    x_BFO = max(x_BFO,y_min*ones(1,N));         % Wacth for saturations            
    x_BFO = min(x_BFO,y_max*ones(1,N));         % Wacth for saturations                    
    x_BFS = max(x_BFS,y_min*ones(1,N));         % Wacth for saturations            
    x_BFS = min(x_BFS,y_max*ones(1,N));         % Wacth for saturations                    
    
    %----------------------------------------------------------------------
                    
    
end

%--------------------------Generate Output---------------------------------
if (recall_algorithm_option == 0)
    x_out = x_WTA;
elseif (recall_algorithm_option == 1)
    x_out = x_BFO; 
elseif (recall_algorithm_option == 2)
    x_out = x_BFS;
else
    error('Unknown recall algorithm');
end
%--------------------------------------------------------------------------