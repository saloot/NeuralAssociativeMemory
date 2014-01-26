%==========================================================================
%*********************FUNCTION: faulty_recall_step*************************
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
% pattern_neur_noise: The maximum amount of noise a pattern neuron will "suffer" from
% const_neur_noise: The maximum amount of noise a constraint neuron will "suffer" from
% const_update_threshold: The update threshold for the constraint neurons
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% x_out: The output pattern of the recall process
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function gets the specifying parameters of a bipartite
% neural network as well as a corrupted patterns and performs a recall step
% of our proposed algorithm in the ISIT 2013 paper. The step is defined as 
% the process of iteratively correcting the error until either there is no
% errors left or a maximum number of iterations has been executed.

% The recall process utilizes neurons that are contaminated with internal
% noise, i.e. a noise parameter affects their decision making during the
% iterative decoding process. The noise is a uniformly distributed random 
% variable between [-a,a], where "a" is the noise parameter and is 
% specified by "pattern_neur_noise" and "const_neur_noise" for pattern and 
% constraint neurons, respectively.

% There are three options for the recall algorithm: bit flipping (original 
% and simplified) as well as the winner-take-all approach. The selection 
% between the two is determined by the recall_algorithm_option flag.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function [x_out,itr_total] = faulty_recall_step(W,x,N,gamma_BFO,gamma_BFS,y_min,y_max,recall_algorithm_option,pattern_neur_noise,const_neur_noise,const_update_threshold)

%============================INITIALIZATION================================
exit_flag = 0;        
itr = 0;                
x_WTA = x;
x_BFO = x;
x_BFS = x;
[N_const,~] = size(W);
%==========================================================================
        


while (exit_flag == 0)
    itr = itr + 1;           
                                                                                                   
    %----------------------Update Constraint Nodes-------------------------
    v = -const_neur_noise + 2*const_neur_noise*rand(N_const,1);     % The noise that affects constraint neurons
    
    c_temp_WTA = W*x_WTA' + v;                                      % Find which constraint nodes are violated.                                       
    c_temp_BFO = W*x_BFO' + v;                                      % Find which constraint nodes are violated.                                       
    c_temp_BFS = W*x_BFS' + v;                                      % Find which constraint nodes are violated.                                       
    
    for iii = 1:N_const                
        if (abs(c_temp_WTA(iii))<const_update_threshold)                        % Take care of numerical issues                        
            c_temp_WTA(iii) = 0;                                
        end        
        
        if (abs(c_temp_BFO(iii))<const_update_threshold)                        % Take care of numerical issues                        
            c_temp_BFO(iii) = 0;                                
        end        
        
        if (abs(c_temp_BFS(iii))<const_update_threshold)                        % Take care of numerical issues                        
            c_temp_BFS(iii) = 0;                                
        end        
    end
    
    c_WTA = sign(c_temp_WTA);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    c_BFO = sign(c_temp_BFO);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    c_BFS = sign(c_temp_BFS);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.
    
    %----------------------------------------------------------------------
    
    %-----------------------Check for Convergence--------------------------
    if ( ((norm(c_WTA) < 0.0001)&&(norm(c_BFO) < 0.0001)&&(norm(c_BFS) < 0.0001) )||(itr > 80))                                                            
        if (itr<=40)
            itr_total = itr;
        else
            itr_total = 10000;
        end
        
        exit_flag = 1;                                                  
        break;
    end    
    %----------------------------------------------------------------------                                                                       

    %---------------------Update Pattern Nodes-----------------------------
    x_temp_WTA = c_WTA'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    x_temp_BFO = c_BFO'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    x_temp_BFS = c_BFS'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                    
    %----------------------------------------------------------------------                                                                       
    
    %-------------------Bit Flipping Original------------------------------                
    u = -pattern_neur_noise + 2*pattern_neur_noise*rand(1,N);   % The noise that affects pattern neurons    
    
    x_temp2 = u + c_BFO'*sign(W)./sum(abs(sign(W)));            % x_temp2 is the raw feedback received by each pattern node.                                                       
        
    tempp = zeros(1,N);                                        
    for jjj = 1:N                                    
        if ( abs(x_temp2(jjj)) >= gamma_BFO)                    % Update all the nodes that receive a lot of feedback from their neighbors.                                                                              
            tempp(jjj) = sign(x_temp2(jjj));                                                        
        else                        
            tempp(jjj) = 0;                            
        end        
    end    
    
    x_BFO = x_BFO - tempp;                
    %----------------------------------------------------------------------
                       
    %-----------------------Winner-Take-All--------------------------------
    temp_var = u + abs(c_WTA')*abs(sign(W))./sum(abs(sign(W)));                % Normalize the feedback by the norm-1 of the outgoing edges    
    [~,inde] = max(temp_var);                           % Find the node(s) with the maximum input sum.                                                                                                                    
    x_WTA(inde)=x_WTA(inde)-sign(x_temp_WTA(inde));             % Update the code nodes.                                        
    %----------------------------------------------------------------------
    
    
    %-------------------Bit Flipping Simplified----------------------------
    x_temp2 = u + (c_BFS'*W)./sum(abs(W));               % x_temp2 is the raw feedback received by each pattern node.                                                       
        
    tempp = zeros(1,N);                                        
    for jjj = 1:N                                    
        if ( abs(x_temp2(jjj))>= gamma_BFS)          % Update all the nodes that receive a lot of feedback from their neighbors.                                                                              
            tempp(jjj) = sign(x_temp_BFS(jjj));                                                        
        else                        
            tempp(jjj) = 0;                            
        end        
    end    
    
    x_BFS = x_BFS - tempp;           
    %----------------------------------------------------------------------
        
        
    %--------------Saturate the Values from Above and Below----------------
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