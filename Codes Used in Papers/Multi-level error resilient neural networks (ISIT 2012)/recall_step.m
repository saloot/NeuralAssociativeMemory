%==========================================================================
%***********************FUNCTION: recall_step******************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% W: The connectivity matrix of the neural graph
% x: The corrupted pattern
% b: The constraints needed to be satisfied
% x_s: The special solution of the network. We must have W*x_s = b
% learn_thr: A threshold that determines the given pattern is a new pattern or it has been learned before and its deviation of constraints is just due to noise.
% N: The number of pattern nodes in the graph
% N_const: The number of constraint nodes in the graph 
% alph: The learning step for reducing the projection of patterns on the weight matrix during the simultaneous learn and recall operation.
% betta: It is obsolete and not used anymore.
% lambda: The Lagrange multiplier for the sparsity penalty function
% gamma_bit: The threshold that determines the margin for majority voting in the bit-flipping algoroithm 
% sigm: The constant which determines how sharp the sparsity penalty function is
% algorithm_option: Selects between bit-flipping or winner-take-all algorithm in the recall phase
% convergence_flag: It is a flag that determines if the algorithm stops only after a maximum number of iterations or stops as soons as convergence is reached
% max_noise_amp: The maximum amplitude of the generated noise. It should be an integer.
% err_bits: The number of initial erroneous nodes 
% y_min: The minimum value that a pattern node can have
% y_max: The minimum value that a pattern node can have
% nois: The added noise to tge given pattern
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% x_final: The output pattern of the recall process
% learn_flag: Will be set to one if a learning round has been executed instead of the normal recall procedure as the given pattern had not been learned before.
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function gets the specifying parameters of a bipartite
% neural network as well as a corrupted patterns and performs a recall step
% of our proposed algorithm in the ITW paper. The step is defined as the
% process of iteratively correcting the error until either there is no
% errors left or a maximum number of iterations has been executed.

% There are two options for the recall algorithm: bit flipping and
% winner-take-all approach. The selection between the two is determined by
% the algorithm_option flag.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function [x_final,learn_flag] = recall_step(W,x,b,x_s,learn_thr,N,N_const,alph,betta,lambda,gamma_bit,sigm,algorithm_option,convergence_flag,max_noise_amp,err_bits,y_min,y_max,nois)

%============================INITIALIZATION================================
exit_flag = 0;        
itr = 0;                
%==========================================================================
        
while (exit_flag == 0)
    itr = itr + 1;           
                                                                                                   
    %----------------------Update Constraint Nodes-------------------------
    c_temp = W*x'-b;                                      % Find which constraint nodes are violated.                                       
    for iii = 1:N_const                
        if (abs(c_temp(iii))<.001)                        % Take care of numerical issues                        
            c_temp(iii) = 0;                                
        end        
    end    
    c = sign(c_temp);                                     % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.    
    %----------------------------------------------------------------------

            
            
    learn_flag = 0;                                       % Reset the learning flag
            
            
    %----------------------Short Learning Phase----------------------------            
    if ((norm(W*x'-b)>learn_thr)&(itr == 1))                   % If the given pattern is very far away from the subspace, modify the weights in its favor
    
        %------------------Some Necessary Initializations------------------
        learn_flag = 1;                    
        % mu_itr = 0;
        % xx = x-x_s;
        %------------------------------------------------------------------

        %----------------Gradually Adjust the Weight Matrix----------------
        % while ((norm (W*xx') > learn_thr)&(mu_itr<100))
        %                     mu_itr = mu_itr + 1;
        %                     
        %                     for j = 1:N_const                                   % For each constraint node, do the following steps to update the weights  attached to this constraint node.
        %                         y(j) = W(j,:)*xx';                              % Calculate the projection of x on the vector W(j,:).
        %                         
        %                         %---Update the Weights in the Correct Direction----
        %                         if (norm(W(j,:)) > 0)
        %                             W(j,:) = W(j,:) - alph * y(j) * ( xx - (y(j)*W(j,:)/(norm(W(j,:)))^2) );  % Update the corresponding weight vector.                               
        %                         else
        %                             W(j,:) = W(j,:) - alph * y(j) * (xx);
        %                         end
        %                         %--------------------------------------------------
        %                         
        %                         %----------Update Towards Sparsity As Well---------
        %                         W(j,:) = W(j,:) - lambda* 2*tanh(sigm*W(j,:)).*(1-(tanh(sigm*W(j,:))).^2);
        %                         %--------------------------------------------------                                                
        % 
        %                         xx = xx - betta*y(j)*W(j,:);      % Update variable nodes    
        %                     end
        %                 
        %                 end
        %             %--------------------------------------------------------------
    
        
    else                        % Else, de-noise the noisy already-learnt pattern                
        
        %---------------------Update Pattern Nodes---------------------                            
        x_temp = c'*W;                                      % x_temp is the raw feedback received by each pattern node.                                                
        temp_var = abs(x_temp)./sum(abs(W));                % Normalize the feedback by the norm-1 of the outgoing edges
        [~,inde] = max(temp_var);                           % Find the node(s) with the maximum input sum.
                                                                        
        if (algorithm_option == 0)                    
            
            %-----------------------Winner-Take-All------------------------
            x(inde)=x(inde)-sign(x_temp(inde));             % Update the code nodes.                            
            %--------------------------------------------------------------                
            
        else                        
            %-------------------Bit Flipping-------------------------------
            x_temp2 = abs(c'*W)./sum(abs(W));               % x_temp2 is the raw feedback received by each pattern node.                                   
            tempp = zeros(1,N);                            
            for jjj = 1:N                            
                if ( abs(x_temp2(jjj))> gamma_bit)          % Update all the nodes that receive a lot of feedback from their neighbors.                                                       
                    tempp(jjj) = sign(x_temp(jjj));                    
                else
                        tempp(jjj) = 0;                    
                end                
            end            
            x = x - tempp;            
            %--------------------------------------------------------------    
                       
        end
        
        x = max(x,y_min*ones(1,N));         % Wacth for saturations        
        x = min(x,y_max*ones(1,N));         % Wacth for saturations                    
        %------------------------------------------------------------------
    end
                
    %-----------------------Check for Convergence--------------------------
    if (convergence_flag == 0)              % If 0, do max_itr iterations and quite                                                    
        if (itr >=max_itr)                                                       
            exit_flag = 1;                                                            
        end        
    else                                    % Else, remain in the loop until convergence or exceeding 1000 iterations.    
        if ((learn_flag == 1)||(norm(c) < 0.001)||(itr > 20*max_noise_amp*err_bits))                                                
            exit_flag = 1;                                               
        end        
    end    
    %------------------------------------------------------------------                       
end

x_final = x;