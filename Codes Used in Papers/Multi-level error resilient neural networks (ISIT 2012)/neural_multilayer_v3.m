%%=========================================================================
%********************************READ ME***********************************
%==========================================================================

%--------------------------------Summary-----------------------------------
% This file contains the code for simulating the error correcting procedure
% in a two layer bipartite non-binary neural associative memory. This network
% enforces two sets of constraints: one set on sub-patterns of a long
% pattern and the other set is enforced globally on the collection of
% sub-patterns.
% As input, the code gets the connectivity matrix of the graph (both layers),
% the number of initial erroneous bits, the maxium amplitude of noise, 
% the maximum and minimum values for the firing rate of pattern neurons. 
% The code then evaluate the performance of the given network by performing 
% simulations by generating a number of randomly generated patterns (that 
% satisfy the set of constraints specified by the connectivity matrix) and
% then initializing the network with these patterns plus some noise. 

% The network first corrects any errors in the sub-patterns and then, make
% an attemp to correct the remaining errors at a higher level (similar to
% natural languages). 

%-------------Simultaneous Learn and Recall-------------
%-------------------------------------------------------
% In this version we only focus on recalling the set of patterns that the
% network has been trained with before. Becuase now that we have L smaller
% graphs for the sub-patterns and a big global network, it is very unlikely
% to find a set of patterns that all the networks have been trained with,
% specially if the number of patterns used in the training set is small.
% And since we aim to investigate the gain we get in error correcting
% performance, we tend to stick to the patterns in the training set.

% However, in general it is possible to choose input patterns arbitrily and
% in reponse to the given patterns, if they are one of the patterns that
% the network had memorized in the learning phase, i.e. the weight matrix
% times this patterns has a small value, then the algorithm tries to
% eliminate noise and find the pattern which is orthogonal to the weight
% matrix. 

% However, if it is a "new" pattern, then the network do a short learning
% phase and adapts the weight matrix so that it becomes orthogonal to the
% given pattern.
%-------------------------------------------------------

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%----------------------HOW THE CODE WORKS IN DETAILS-----------------------
%--------------------------------------------------------------------------
% The code first loads some of the already initialized variables from a file.
% The code also reads the weight matrix for the corresponding neural graph
% from the specified files. 
% Then, in order to evaluate the performance of this network, a random 
% intger-valued noise is produced in which some of its elements (depending
% on the number of erroneous nodes we want to simulate) is an integer in 
% the interval [-max_noise_amp,max_noise_amp] excluding zero and all other
% elements are equal to zero.

% In reponse to the given patterns, if they are one of the patterns that
% the network had memorized in the learning phase, i.e. the weight matrix
% times this patterns has a small value, then the algorithm tries to
% eliminate noise and find the pattern which is orthogonal to the weight
% matrix. This is done by  iterating and updating the state of neurons
% according to the rules mentioned in the paper: "Exponential Pattern 
% Retrieval Capacity with Non-Binary Associative Memory" by R. K. Kumar,
% A. H. Salavati and A. Shokrollahi. 

% However, if it is a "new" pattern, then the network do a short learning
% phase and adapts the weight matrix so that it becomes orthogonal to the
% given pattern. The learning rules are based on the paper 
% "Neural Nets for Dual Subspace Pattern Recognition Method". The short
% learning phase is continued until either a maximum number of iterations
% has been done or the matrix is fully adjusted to be orthogonal to the new
% pattern.

% If the convergence flag is set to 1, then the code iterates until
% convergence is reached or max_itr iterations is done. If the flag is set
% to zero, then the code does max_itr iterations, no matter if the
% convergence is reached or not. After finishing the iterations, the code
% checks if the final state of neurons is the all-zero vector. If not, an
% error is declared and pattern_error_count is increased. The overall pattern_error_count
% divided by the total number of iterations then represent the performance
% of the error correcting algorithm

% The results will be stored in a file on the specified location in the end
% of the code. The progress is also displayed on the screen.




%--------------------------------------------------------------------------
%-------------------------------FURTHER NOTES------------------------------
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function neural_multilayer_v3(error_bits_in,max_noise_amp,max_simulated_instances,N,K,L)
%=============================INITIALIZATION===============================

%-----------------Load Previously Initialized Variables--------------------
load(['neural_learn_parameters_v3_N_',num2str(N),'_K_',num2str(K),'.mat']);
%--------------------------------------------------------------------------


%-------------------------Simulation Parameters----------------------------    
convergence_flag = 1;                   % 0 means you repeat for max_itr times. 1 means you wait untill convergence or do max_itr iteration
max_itr = 10;                           % This is the maximum number of iteration for the simulation in convergence_flag = 0
algorithm_option = 1;                   % If it is equal to zero, we use winner-take-all algorithm. If 1, bit flipping is used.     
%--------------------------------------------------------------------------
  
%-----------------------------Neural Parameters----------------------------
y_max = 200;                            % This is the maximum value of code nodes.
y_min = -200;                           % This is the minimum value of code nodes.   
gamma_bit = 0.8;                        % gamma is the update threshold for the bit-flipping algorithm.                
learn_thr_0 = 10;                       % learn_thr_0 is the initial threshold used in the simultaneous learn and recall phase. 
%--------------------------------------------------------------------------


%------------------------Other Initializations-----------------------------
subpattern_error_count = zeros(1,length(error_bits_in));     % This tracks the number of sub-pattern errors (to compute the block error rate).
pattern_error_count = zeros(1,length(error_bits_in));        % This tracks the number of pattern errors (to compute the block error rate).
pattern_error_count_primary = zeros(1,length(error_bits_in));% This tracks the number of pattern errors (to compute the block error rate).
sub_bit_error_count = zeros(1,length(error_bits_in));        % This tracks the number of bit errors (to compute the bit error rate).
bit_error_count = zeros(1,length(error_bits_in));            % This tracks the number of bit errors (to compute the bit error rate).
bit_error_count_primary = zeros(1,length(error_bits_in));    % This tracks the number of pattern errors (to compute the block error rate).

a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 

N_tot = N*L;                            % This is the total length of the patterns
N_const_tot = N_tot-K;                  % This is the total number of constraints
learn_no = 0;                           % This is te number of times the learning phase has been executed during the recallprocess. 
learn_no_tot = 0;
success_count = zeros(1,length(error_bits_in));
tic                                     % Reset the timer for display purposes.
%--------------------------------------------------------------------------

%==========================================================================


%%
%=======================READ THE WEIGHT MATRIX=============================
W_sub = zeros(N_const*L,N);

%-------------Load the SubPattern Weight Matrix from the File--------------
for l = 1:L
    fid = fopen(['neural_learn_weight_matrix_multilayer_v3_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(l),'.txt'], 'r');
    W = fscanf(fid, '%f',[N,N_const]);
    W = W';
    fclose(fid);
    
    %----------Check Whether the Initial Learning Phase Is Done------------
    [m,~] = size(W);
    if (m < N_const)
        error('Learning phase is not complete yet!');
    else    
        W_sub( (l-1)*N_const+1:l*N_const,:) = W;
        display(' ');
        display('---------------------------------------------------------');
        display('Learning phase finished successfully');
        display('---------------------------------------------------------');
        display(' '); 
    end
    %----------------------------------------------------------------------

end
    %--------------------------------------------------------------------------

%------------Load the Total Pattern Weight Matrix from the File------------
fid = fopen(['neural_learn_total_weight_matrix_multilayer_v3_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'.txt'], 'r');
W_tot = fscanf(fid, '%f',[N_tot,N_const_tot]);
W_tot = W_tot';
fclose(fid);

    %----------Check Whether the Initial Learning Phase Is Done------------
    [m,~] = size(W_tot);
    if (m < N_const_tot)
        error('Total learning phase is not complete yet!');
    else    
        display(' ');
        display('---------------------------------------------------------');
        display('Total learning phase finished successfully');
        display('---------------------------------------------------------');
        display(' '); 
    end
    %----------------------------------------------------------------------


%--------------------------------------------------------------------------



%----------------------Calculate the Special Answers-----------------------
x_s = 1+floor(2*rand(1,N));
b_sub = W_sub*x_s';        

x_s_tot = 1+floor(2*rand(1,N_tot));
b_tot = W_tot*x_s_tot';        
%--------------------------------------------------------------------------

%==========================================================================




%%
%================================MAIN LOOP=================================
for noise_itr = 1:length(error_bits_in)    
    
    err_bits = error_bits_in(noise_itr);   % Determine the number of erroneous bits.
           
    learn_thr = learn_thr_0 *err_bits;
    
    
    for net_simul_itr = 1:max_simulated_instances   % Simulate the error correction procedure for the given ensemble.                                
        
        %--------------------------Initialization--------------------------
        x_tot = [];                         % Initialize the total estimated pattern vector
        pattern_tot = [];                   % Initialize the total pattern vector
        %------------------------------------------------------------------
        
        %---------------------------Generate Noise-------------------------
        nois_tot = zeros(1,N_tot);          % This is the noise added to the whole pattern of length N*L        
        pp = 1+floor((N_tot-1)*rand(1,err_bits));                
        for h = 1:err_bits        
            nois_tot(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randint);            
        end        
        %------------------------------------------------------------------                                
                
        %--------------------Generate the Pattern Index--------------------
        p = 1+floor((pattern_learn_number-1)*rand);                 % Pick a pattern index at random
        temp = dec2bin(mu_list(p),N-N_const);                       % Transform the index into binary
        message = zeros(1,K);                                       % Generate the message from the index                    
        for j = 1:N-N_const
            message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;        % Generate the message from the index. 48 is the ASCII code for 0.                
        end
        %------------------------------------------------------------------                                
        
        %---------------------Correct Subpatterns First--------------------
        error_flag = 0;                                             % This flag is triggered when we have an uncorrected error at the output of the first level. 
        learn_flag = zeros(1,L);                                       % This flag determines if we executed a learning phase during the recall procedure.
        for l = 1:L
            
            nois = nois_tot((l-1)*N+1:l*N);                         % Find the sub-noise for this subpattern
            G_sub = G((l-1)*K+1:l*K,:);                             % Retrieve the corresponding generator matrix for this particular sub-network
            W = W_sub( (l-1)*N_const+1:l*N_const,:);                % Retrieve the corresponding weight matrix for this particular sub-network
            b = b_sub( (l-1)*N_const+1:l*N_const,:);                % Retrieve the corresponding constraints for this particular sub-network
                        
            pattern = x_s+message*G_sub;                            % Generate the pattern from the message and the special answer
            
            %--------------------------------------------------------------
            
            %--------Initialize the Network with a Noisy SubPattern-------- 
            pattern_tot = [pattern_tot,pattern-x_s];                % Store the subpatterns to make the bigger total pattern
            x = pattern+nois;                                       % Initialize the network with a nosiy version of the subpattern        
            %------------------------------------------------------------------
            
        
        
            %---------------------Iterate Until Convergence--------------------                
            [x_fin,learn_flag(l)] = recall_step(W,x,b,x_s,learn_thr,N,N_const,alph,betta,lambda,gamma_bit,sigm,algorithm_option,convergence_flag,max_noise_amp,err_bits,y_min,y_max,nois);       
            %--------------------------------------------------------------
        
            x_tot = [x_tot,x_fin-x_s];                              % Update the total estimation vector
            
                
            %-------------Calculate Subpatterns Error Rate-----------------
            if ((norm(x_fin-pattern)>.01)&&(learn_flag(l) == 0))       % If there is an error due to the recall process (and not an unlearned pattern) then update error rates.                   
                subpattern_error_count(noise_itr) = subpattern_error_count(noise_itr) + 1;                                                    
                sub_bit_error_count(noise_itr) = sub_bit_error_count(noise_itr) + sum(abs(sign(x_fin-pattern)) );      
                error_flag = 1;                
            end        
            %--------------------------------------------------------------
            
        end
        
        %------------Calculate Error Rate After the First Phase------------
        if ( (error_flag == 1) & (norm(learn_flag) == 0) )
            pattern_error_count_primary(noise_itr) = pattern_error_count_primary(noise_itr) + 1;
            bit_error_count_primary(noise_itr) = bit_error_count_primary(noise_itr) + sum(abs(sign(x_tot-pattern_tot)));
        end   
        if (norm(learn_flag) > 0)                                      % If a learning phase has been executed...         
            learn_no = learn_no + 1;
        end
            
        %------------------------------------------------------------------
               
       
        
        %--------------------Now Correct the Whole Pattern-----------------       
        second_layer_flag = 0;
        if ((norm(learn_flag) == 0) && (error_flag == 1))              % If there is an error at the end of the first level and it is not due to a pattern not in the training set then...
            second_layer_flag = 1;
            x_tot = x_tot + x_s_tot;
            [x_tot_fin,mu_flag_tot] = recall_step(W_tot,x_tot,b_tot,x_s_tot,learn_thr,N_tot,N_const_tot,alph,betta,lambda,gamma_bit,sigm,algorithm_option,convergence_flag,max_noise_amp,err_bits,y_min,y_max,nois_tot);                                                       
        end
        %------------------------------------------------------------------        
        
        %--------------------Calculate Error Rate--------------------------
        if (second_layer_flag == 1) 
            if ((norm(x_tot_fin-x_s_tot-pattern_tot)>.001)&& (mu_flag_tot == 0))
                pattern_error_count(noise_itr) = pattern_error_count(noise_itr) + 1;                                                    
                bit_error_count(noise_itr) = bit_error_count(noise_itr) + sum(abs(sign(x_tot_fin-x_s_tot-pattern_tot)) );              
            elseif (mu_flag_tot == 0)         
                success_count(noise_itr) = success_count(noise_itr) + 1;      
            end
            if (norm(mu_flag_tot) > 0);
                learn_no_tot = learn_no_tot + 1;
            end
        end
           
        %------------------------------------------------------------------
        
        
        %----------------Display Results Every 5 Seconds-------------------              
        if (toc >= 5)                                                             
            display(' ');
            display([' Iteration:',num2str(net_simul_itr)]);
            display([' Block error count after the first phase so far: ',num2str(pattern_error_count_primary(noise_itr))]);
            display([' Block error count after the second phase so far: ',num2str(pattern_error_count(noise_itr))]);            
            display(' ');   
            display([' Bit error rate after the first phase so far: ',num2str(bit_error_count_primary(noise_itr)/(N_tot*net_simul_itr))]);
            display([' Bit error rate after the second phase so far: ',num2str(bit_error_count(noise_itr)/(N_tot*net_simul_itr))]);
            display(' ');   
            display([' Sub-Block error count so far: ',num2str(subpattern_error_count(noise_itr))]);
            display([' Sub-Bit error rate so far: ',num2str(sub_bit_error_count(noise_itr)/(N*net_simul_itr))]);            
            display('-------------------------------');               
            tic                
        end          
        %------------------------------------------------------------------
                            
    end
        
    
    
    display(' ');
    display('---------------------------------------------------------');
    display(['Noise itr = ',num2str(noise_itr),' out of ',num2str(length(error_bits_in))]);
    display(['Block error rate = ',num2str(subpattern_error_count(noise_itr)/(net_simul_itr))]);    
    display(['Bit error rate = ',num2str(sub_bit_error_count(noise_itr)/(N*net_simul_itr))]);    
    display('---------------------------------------------------------');
    display(' ');
    
    
    %------------------Transform Error Count to Error Rate---------------------
    subpattern_error_rate = subpattern_error_count(noise_itr)/(net_simul_itr);
    pattern_error_rate_primary = pattern_error_count_primary(noise_itr)/(net_simul_itr);
    bit_error_rate_primary = bit_error_count_primary(noise_itr)/(N_tot*net_simul_itr);
    pattern_error_rate = pattern_error_count(noise_itr)/net_simul_itr;
    bit_error_rate = bit_error_count(noise_itr)/(N_tot*net_simul_itr);
    %--------------------------------------------------------------------------
    
    %----------------Store the Bit and Pattern Error Rates-----------------
    fid = fopen(['neural_multi_level_results_v3_N_',num2str(N),'_K_',num2str(K),'.txt'], 'a+');  
    fprintf(fid, 'err_bits \t %d \t pattern_error_rate_primary \t %f \t bit_error_rate_primary \t %f \t subpattern_error_rate \t %f \t pattern_error_rate_second_layer \t %f \t bit_error_rate_second_layer \t %f \t success_count \t %d \n',err_bits, pattern_error_rate_primary,bit_error_rate_primary,subpattern_error_rate,pattern_error_rate,bit_error_rate,success_count(noise_itr));
    fclose(fid);
    %----------------------------------------------------------------------
end
                 


%==========================================================================





%%
%=================================SAVE RESULTS=============================

%-----Store the Total Number of Times the Learn Phase Has Been Executed----
fid = fopen(['learn_no_N_',num2str(N),'_K_',num2str(K),'.txt'], 'a+');
fprintf(fid,'err_bits \t %d learn_number \t %d \t learn_number_tot %d \n',error_bits_in(i),learn_no,learn_no_tot);
fclose(fid);
%--------------------------------------------------------------------------

%==========================================================================
