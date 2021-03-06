%==========================================================================
%********************************READ ME***********************************
%==========================================================================

%--------------------------------Summary-----------------------------------
% This piece of code implements the idea mentioned in the paper 
% "Neural Nets for Dual Subspace Pattern Recognition Method" in order to
% learn the matrix W whose columns define the null basis of a set of given
% patterns. There is one modification though, as compared to the original
% idea in the paper, and that is the trick we use to make W sparser. 

% Furthermore, since this code is desined to implement the learning phase
% of multi-level neural network, it gets as input an index corresponding to
% its position in the first level of the whole network architecture. This
% only affects the set of patterns generated for the training set as
% everything else in the learning algorithm is the same.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%----------------------HOW THE CODE WORKS IN DETAILS-----------------------
%--------------------------------------------------------------------------
% This piece of code implements the idea mentioned in the paper 
% "Neural Nets for Dual Subspace Pattern Recognition Method" in order to
% learn the matrix W whose columns define the null basis of a set of given
% patterns. There is one modification though, as compared to the original
% idea in the paper, and that is the trick we use to make W sparser. 

% The code starts by reading the initialized variables from a file. This
% includes the generator matrix as well as the optimization variables.
% Then, it goes over some of the patterns (generated by multiplying G with a
% binary vector of proper length). For each pattern, the weights in W,
% initialized randomly, are updated towards the direction that in the end
% they converge to the dual basis of the subspace defined by G (read the 
% section on parallelization below as well. 

% During the process, we also try to find a sparse solution by assigning a 
% cost to the non-sparse answers. This process is repeated several times in 
% order to gradually converge to the correct answer.

% The code stops once a maximum number of iterations is passed or
% convergence is reached. 

% Once finished, the learnt vectors are added to the end of the appropriate
% file so that it constitute a new row in the constraint matrix. The
% corresponding file will be used later as the connectivity matrix resulted
% by the learning phase. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%----------------------------PARALLELIZATION-------------------------------
%--------------------------------------------------------------------------
% In itself, the idea mentioned above is computationally heavy. So we make 
% the procedure parallel in the sense that in this version, the code only
% learns the connectivity of one constraint node. Thus, if we run several
% instances of the code over several machines, we manage to learn the
% constraint matrix in parallel. 

% However, it is not as easy as it seems for three important issues:
% 1) The case of all-zero solution
% 2) Simultaneous file access problems
% 3) The case of non-zero betta, i.e. orthogonal and non-redudant constraints

% To overcome the first issue, we have seen in the serial version that once
% you start from a sparse bipartite graph without any zero columns or rows,
% the serial version of the algorithm yields a weight matrix without any
% all-zero rows or columns as well. So in the parallel version, we do the
% same. 

% However, one major issue is the way in which we pass this initialized 
% weight matrix to the learning functions. We use a file to store the
% initial version of the weight matrix. Then, each instance of the running
% code will pick one of the rows and use this vector as the initial
% estimation of the constraint node's connectivity. This estimation is
% gradually improved over time according to the method explained above. In
% order to make sure that different instances of the code use different
% rows and not the same row always, we store the index of the processed rows
% in a separate file. Whenever an instance of the code loads the initial
% weight matrix, it goesover this file and look at the last row index
% stored in the file. The code increases this number, add the new index to
% the end of the row index file and use the corresponding row as its
% initial estimation of the connectivity.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%-------------------------------FURTHER NOTES------------------------------
%--------------------------------------------------------------------------
% The second point above is still an issue because some times, two
% instances access the same file at the same time so what they write on the
% output file would be mis-formatted. Furthermore, if they access the row
% index simultaneously, they both process the same initial estimation and
% the final weight vector will most probably be the same (redundant
% constraints).

% The third issue above corresponds to the case that we want to remove the
% projection of the data vector over the basis vectors in order to update
% the new one. For the parallel version, so far we have assumed that this is
% not the case and different initial weight vectors will result in
% different basis vectors. Simulations also show that with good
% probability, most of the basis vectors will be different.
% However, it would be nice to test the algorithm with non-zero values of
% betta as well, i.e. to consider the case of projection removal. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function learn_weights_dual_space_multilayer_v3(N_in,K_in,cluster_index)

%=============================INITIALIZATION===============================

%----------------Load the Saved Initialized Parameters---------------------
load(['neural_learn_parameters_v3_N_',num2str(N_in),'_K_',num2str(K_in),'.mat']);
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
element_thr = 0.001;                    % Weights less than this threshold are set to zero.
MSE = zeros(1,learn_itr_max);           % This stores the Mean Square Error in each learn iteration.
a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------


%------------------Randomly Initialize the Weight Matrix-------------------
W = random_vector(N,2*round(log(N)));
%--------------------------------------------------------------------------

%==========================================================================


%=============================LEARNING PHASE===============================

lambda = 0.005;
%--------------------------Main Learning Loop------------------------------
for learn_itr = 1:learn_itr_max         % Repeat the steps below several times in order to enhance the quality of the learning phase    
    exit_flag = 1;
    
    for k = 1:pattern_learn_number           % For each pattern, repeat the following steps
        
        %----------------------Generate the Patterns-----------------------
        mu = mu_list(k);                                        % Pick a pattern from the training set and learn it.
       
        temp = dec2bin(mu,N-N_const);                           % Transform the pattern index to a binary sequence        
        message = zeros(1,K);                                   
        for j = 1:N-N_const
            message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;    % Generate the message from the index. 48 is the ASCII code for 0.                             
        end
        G_cluster = G((cluster_index-1)*K+1:(cluster_index)*K,:);
        pattern = message*G_cluster;                            % Generate the pattern from the message and the generator matrix.
        %------------------------------------------------------------------
        
        x = pattern;                    % Initialize the variable nodes with the pattern
                             
        y = W*x';                       % Calculate the projection of x on the vector W(j,:).
       
        %-----------Update the Weights in the Correct Direction------------
        if (norm(W) > 0)                
            W = W -alph * y * ( x - (y*W/(norm(W))^2) );  % Update the corresponding weight vector.            
        end        
        %------------------------------------------------------------------
                        
        %----------------Update Towards Sparsity As Well-------------------
%         W = W - lambda* 2*W.*exp(W/10);
        W = W - lambda* 2*W.*(1-(tanh(sigm*W)).^2);% - lambda * 4*W(j,:).*tanh(sigm*W(j,:).^2).*(1-(tanh(sigm*W(j,:).^2)).^2);
%         lambda = max(lambda + 0.001*(sum(abs(W))-deg_row_W),0);
%         lambda = min(lambda,0.01);
        %------------------------------------------------------------------                                                                    
                    
        %--------------------Check for Numerical Errors--------------------
        if (sum(isnan(W))>0)                  
            error('W is not defined!');                            
        end        
        %------------------------------------------------------------------                                    

        %-----------------------Calculate the MSE--------------------------
        MSE(learn_itr) = MSE(learn_itr)+norm(W*pattern');
        %------------------------------------------------------------------
        
        %--------------------In Case of Not Convergence--------------------        
        if (norm (W*pattern') > epsilon)
            exit_flag = 0;
        end
        %------------------------------------------------------------------
        
    end
     
    %-------------Remove Small Elements From the Weight Vector-------------
    for j = 1:N
        if (abs(W(j))<element_thr)            
            W(j) = 0;        
        end
    end
    %----------------------------------------------------------------------
                 
    %--------------------------Display Progress----------------------------
    if (mod(learn_itr,1) == 0)        
        norm(W)
        MSE(learn_itr)              % The projection of the patterns over the weight vector.
        sum(abs(sign(W)))           % Number of non-zero elements in the weight vector        
    end
    %----------------------------------------------------------------------
        
    %-----------------------Check for Convergence--------------------------
    if (exit_flag == 1)
        MSE = MSE(1:learn_itr);
        break;
    end
    %----------------------------------------------------------------------

    
end
%--------------------------------------------------------------------------

%==========================================================================



%%
%==============================SAVE RESULTS================================
% plot(MSE)
fid = fopen(['neural_learn_weight_matrix_multilayer_v3_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(cluster_index),'.txt'], 'a');
fprintf(fid, '%f \t',W');
fprintf(fid, '\n');
fclose(fid);

fid = fopen(['neural_learn_multilayer_v3_MSE_N_',num2str(N),'_K_',num2str(K),'.txt'], 'a');
fprintf(fid, '%f \t',MSE);
fprintf(fid, '\n');
fclose(fid);

fid = fopen(['neural_multilayer_v3_learn_itr_N_',num2str(N),'_K_',num2str(K),'.txt'], 'a');
fprintf(fid, '%d \t',learn_itr);
fprintf(fid, '\n');
fclose(fid);
%==========================================================================
