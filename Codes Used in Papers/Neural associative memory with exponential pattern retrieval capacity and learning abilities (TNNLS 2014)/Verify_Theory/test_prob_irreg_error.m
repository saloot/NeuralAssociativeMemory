%==========================================================================
%******FUNCTION: test_prob_irreg_error(N,N_const,d,lambda,error_bits)******
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of variable (left) nodes.
% N_const: The number of check (right) nodes.
% gamma: The update threshold for the bit-flipping algorithm.
% d: The degree distribution of variable nodes.
% lambda: The distribution of degrees for variable nodes.
% error_bits: number of initial errors
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% error_rate_theory: The average probability of error in the first iteration derived theoretically
% error_rate_simul: The average probability of error overall the whole iterations obtained from simulations
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function verifies the correctness of the bounds derived
% for the error correction performance of a bipartite irregular neural
% graph when used in the bit-flipping error correction method suggested in
% the ITW paper. It generates random ensembles to simulate the error
% correction procedure and calculate the probability of recall error in
% simulations. Then, the code compares this value with an upper bound
% derived in the progress report of 17-30 March 2012.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================

function [error_rate_theory,error_rate_simul] = test_prob_irreg_error(N,N_const,d,lambda,error_bits)

%----------------Check the Vailidity of the Input Parameters---------------
if (norm(sum(lambda) - 1)>0.001)
    error('Probability distribution must add up to 1');
end
if (length(lambda) ~= length(d))
    error('Invalid input!');
end
%--------------------------------------------------------------------------

%-------------------------------Initialization-----------------------------
deg_ave = sum(d.*lambda);                       % The average degree of variable nodes
algorithm_option = 1;                           % Choose the recall algorithm: 1 indicates bit-flipping and 0 shows winner-take-all
gamma_bit = 0.8;                                % The update threshold in the bit-flipping algorithm
max_simulated_instances = 1000;                 % The number of simulated noisy patterns in each ensemble
max_no_ensembles = 10;                          % The number of random ensembles generated by the algorithm
max_noise_amp = 1;                              % The maximum noise value
error_count = zeros(1,length(error_bits));      % Number of recall errors for each noise iteration
%--------------------------------------------------------------------------


for noise_itr = 1:length(error_bits)
    err_bits = error_bits(noise_itr);

    for test_itr = 1:max_no_ensembles
        
        %----------Generate a Weighted Irregular bipartite Graph-----------
        [H,~] = bipartite_right_irregular(N,N_const,d,lambda);
        H = rand(N_const,N).*H;                     % Adjust random weights
        %------------------------------------------------------------------
        
        for net_simul_itr = 1:max_simulated_instances   % Simulate the error correction procedure for the given ensemble.
                 
            %-----------------------Generate Noise-------------------------    
            nois = zeros(1,N);        
            pp = 1+floor((N-1)*rand(1,err_bits));                
            for h = 1:err_bits        
                nois(pp(h)) =(1+floor((max_noise_amp)*rand))*((-1)^randint);            
            end                
            %--------------------------------------------------------------
    
                  
            %--------Initialize the Network with a Noisy Pattern-----------                  
            x = zeros(1,N);                             % Initialize the network with the zero-pattern
            syndrome = H*nois';                         % Syndrome is the result of linear summation of noise input at the constraint nodes.         
            %--------------------------------------------------------------
        
            %------------------Iterate Until Convergence-------------------
            exit_flag = 0;        
            itr = 0;        
            while (exit_flag == 0)                                                                                
                itr = itr + 1;
                                
                %------------------Update Constraint Nodes-----------------
                c_temp = syndrome-H*x';     % Find which constraint nodes are violated.
                
                    
                for iii = 1:N_const            
                    if (abs(c_temp(iii))<.001)      % Take care of numerical issues                
                        c_temp(iii) = 0;                    
                    end                
                end                        
                c = sign(c_temp);       % c is the vector of constraint nodes and non-zero positions indicate a constraint violation.            
                %----------------------------------------------------------

            
                %-------------------Update Pattern Nodes-------------------
            
                x_temp = c'*H;    %x_temp is the raw feedback received by each pattern node.            
                temp_var = zeros(1,N);
                                
                for jjj = 1:length(temp_var)            
                    temp_var(jjj) = abs(x_temp(jjj))/sum(H(:,jjj));     % Divide the number of feedbacks received by the out-going degree                                            
                end            
                [~,inde] = max(temp_var);           % Find the node(s) with the maximum input sum.                                                                
            
                if (algorithm_option == 0)           
                    %------------------Winner-Take-All---------------------            
                    x(inde)=x(inde)+sign(x_temp(inde)); % Update the code nodes.                
                    %------------------------------------------------------            
                else                
                    %---------------------Bit Flipping---------------------
                    x_temp2 = abs(c)'*H;                
                    tempp = zeros(1,N);
                
                    for jjj = 1:N                
                        if ( abs(x_temp2(jjj))> gamma_bit*sum(H(:,jjj)))           % Update all the nodes that receive a lot of feedback from their neighbors.                                        
                            tempp(jjj) = sign(x_temp(jjj));                                                    
                        end                    
                    end                
                    x = x + tempp;               
                    %------------------------------------------------------
                end
                                                            
                %--------------------Check for Convergence-----------------                                  
                if ((norm(c) < 0.0001)||(itr > 2*max_noise_amp*err_bits))            
                    exit_flag = 1;                        
                end                                    
                %----------------------------------------------------------
            end
            %---------Calculate Error Rate and Display Progress------------
                
            if (norm(c)>.001)
                error_count(noise_itr) = error_count(noise_itr) + 1;
            end        
            %--------------------------------------------------------------
        end
%         error_count(noise_itr)
    end
end

%-------------------------Wrap up the Results------------------------------
error_rate_theory = Error_prob_ireg(N,N_const,gamma_bit,d,lambda,error_bits);
error_rate_simul = error_count/test_itr/net_simul_itr;
%--------------------------------------------------------------------------
% norm(error_rate_theory-error_rate_simul)


%------------------------Plot the Results----------------------------------
plot(error_bits,log(error_rate_theory))
hold on 
plot(error_bits,log(error_rate_simul),'r')
legend('Theoretical upper bound','Simulation results');
%--------------------------------------------------------------------------