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


function [x_out] = recall_step_clustered_scheduled(W,x,N,gamma_BFO,y_min,y_max)

%============================INITIALIZATION================================
exit_flag = 0;        
itr = 1;                
[N_const,~] = size(W);
x_out = x;
%==========================================================================
        

%----------------------Update Constraint Nodes-------------------------    
c_temp = W*x' ;                                      % Find which constraint nodes are violated.                                       
            
for iii = 1:N_const                        
    if (abs(c_temp(iii))<.075)                        % Take care of numerical issues                            
        c_temp(iii) = 0;                                        
    end   
end

c = sign(c_temp);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    
if ( (norm(c) < 0.0001) || (itr > N))                                                                
    return
end    
    
while (exit_flag == 0)

    itr = itr + 1;           
                                                                                                       

    %---------------------Update Pattern Nodes-----------------------------    
    x_temp_BFO = c'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                        
    changed_flag = 0;
    x_temp2 = (c')*((W))./(abs((W)));               % x_temp2 is the raw feedback received by each pattern node.                                                       
                
    if ( x_temp2(itr)>= gamma_BFO)         % Update all the nodes that receive a lot of feedback from their neighbors.                                                                              
       x(itr) = x(itr)-sign(x_temp_BFO(itr));                                                                
       changed_flag = 1;
    end            
    %----------------------------------------------------------------------
                       
                                      
                    
    
    
    %----------------------Update Constraint Nodes-------------------------    
    c_temp = W*x' ;                                      % Find which constraint nodes are violated.                                       
        
    for iii = 1:N_const                        
        if (abs(c_temp(iii))<.075)                        % Take care of numerical issues                        
            c_temp(iii) = 0;                                
        end                        
    end
        
    c = sign(c_temp);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.                       
    %----------------------------------------------------------------------
    
    %-----------------------Check for Convergence--------------------------
    if ( (norm(c) < 0.0001) || (itr > N))                                                            
        exit_flag = 1;                                                  
        break;
    else
        if (changed_flag)
            x(itr) = x(itr)+sign(x_temp_BFO(itr)); 
        end
    end        
    %----------------------------------------------------------------------                                                                       
    
end

%--------------------------Generate Output---------------------------------
% x = max(x,y_min*ones(1,N+1));         % Wacth for saturations               
% x = min(x,y_max*ones(1,N+1));         % Wacth for saturations                    
x_out = x;
%--------------------------------------------------------------------------