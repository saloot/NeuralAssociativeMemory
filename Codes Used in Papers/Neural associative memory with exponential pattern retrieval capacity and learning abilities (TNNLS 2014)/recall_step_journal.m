%==========================================================================
%***********************FUNCTION: recall_step******************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% W: The connectivity matrix of the neural graph
% x: The corrupted pattern
% b: The constraints needed to be satisfied
% N: The number of pattern nodes in the graph
% N_const: The number of constraint nodes in the graph 
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
% max_noise_amp: The maximum amplitude of the generated noise. It should be an integer.
% err_bits: The number of initial erroneous nodes 
% y_min: The minimum value that a pattern node can have
% y_max: The minimum value that a pattern node can have
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% x_WTA: The output of the Winner-Take-All algorithm
% x_BFO: The output of the original bit-flipping algorithm
% x_BFS: The output of the simplified bit-flipping algorithm
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function gets the specifying parameters of a bipartite
% neural network as well as a corrupted patterns and performs a recall step
% with three different algorithms the winner-take-all, the original bit 
% flipping (introduced in our ITW 2011) paper and the simplified 
% bit-flipping (introduced in our ISIT 2012 paper). All three outputs are
% generated at the output for the comparison. 

%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function [x_WTA,x_BFO,x_BFS,results_first_itr] = recall_step_journal(W,x,b,N,N_const,gamma_BFO,gamma_BFS,max_noise_amp,err_bits,y_min,y_max)

%%
%============================INITIALIZATION================================
exit_flag = 0;        
itr = 0;                
x_WTA = x;
x_BFO = x;
x_BFS = x;
b = zeros(N_const,1);
% results_first_itr = zeros(1,3);
%==========================================================================
        

%%
while (exit_flag == 0)
    itr = itr + 1;           
                                                                                                   
    %----------------------Update Constraint Nodes-------------------------
    c_temp_WTA = W*x_WTA'-b;                                      % Find which constraint nodes are violated.                                       
    c_temp_BFO = W*x_BFO'-b;                                      % Find which constraint nodes are violated.                                       
    c_temp_BFS = W*x_BFS'-b;                                      % Find which constraint nodes are violated.                                       
    
    for iii = 1:N_const                
        if (abs(c_temp_WTA(iii))<.001)                        % Take care of numerical issues                        
            c_temp_WTA(iii) = 0;                                
        end        
        
        if (abs(c_temp_BFO(iii))<.001)                        % Take care of numerical issues                        
            c_temp_BFO(iii) = 0;                                
        end        
        
        if (abs(c_temp_BFS(iii))<.001)                        % Take care of numerical issues                        
            c_temp_BFS(iii) = 0;                                
        end        
    end
    
    c_WTA = sign(c_temp_WTA);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    c_BFO = sign(c_temp_BFO);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    c_BFS = sign(c_temp_BFS);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.    
    %----------------------------------------------------------------------
      
    %-----------------------Check for Convergence--------------------------
    if ( ((norm(c_WTA) < 0.001)&&(norm(c_BFO) < 0.001)&&(norm(c_BFS) < 0.001) )||(itr > 20*max_noise_amp*err_bits))                                                            
        exit_flag = 1;                                                  
        break;
    end    
    %----------------------------------------------------------------------

    %-------------------------Update Pattern Nodes-------------------------
    x_temp_WTA = c_WTA'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    x_temp_BFO = c_BFO'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    x_temp_BFS = c_BFS'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    
    
    %---------------------------Winner-Take-All----------------------------
    temp_var = abs(c_WTA'*sign(W))./sum(abs(sign(W)));                % Normalize the feedback by the norm-1 of the outgoing edges    
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
                       
                
    %---------------------Bit Flipping Simplified--------------------------                
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
    
    if (itr == 1)
        results_first_itr = [x_WTA;x_BFO;x_BFS];
    end
    
    %----------------------------------------------------------------------
                    
                     
end
