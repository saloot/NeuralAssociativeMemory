%==========================================================================
%**********************FUNCTION: deg_dribution*****************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% alpha0: the learning step in the learning algorithm
% beta0: the sparsity penalty coefficient in the learning algorithm
% theta0: sparsity threshold in the learning algorithm
% index_in: the index of the weight matrix to be analyzed
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% deg_row: the degree of constraint nodes, from 1 to N
% deg_col: the degree of pattern nodes, from 1 to N-K
% lambda: the distribution of column degrees. lambda(i) represents the fraction of pattern nodes with degree i.
% rho: the distribution of row degrees. rho(i) represents the fraction of constraint nodes with degree i.
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a neural network and then reads
% the neural networks learned by the any of the proposed learning
% algorithms from a file. It then calculates the row and column degree
% distributions. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function [deg_row,deg_column,lambda,rho] = deg_dribution_ind(N,K,alpha0,beta0,theta0,index_in)

%------------------------------INITIALIZATION------------------------------
deg_column = [0:N-K];
deg_row = [0:N];
lambda = zeros(1,N-K+1);
rho = zeros(1,N+1);
%--------------------------------------------------------------------------

%----------------Read the Weight Matrix of the Graph-------------------
fid = fopen(['/scratch/amir/ITW_Journal/Learn_Results/N_',num2str(N),'_K_',num2str(K),'/Weigh_matrix_alpha_',num2str(alpha0),'_beta_',num2str(beta0),'_theta_',num2str(theta0),'_index_',num2str(index_in),'.txt'], 'r');    
if (fid == -1)        
    error('No input matrix!');    
end
W = fscanf(fid, '%f',[N,N-K]);
W = W';
fclose(fid);
%--------------------------------------------------------------------------
        
%-------------------Compute Column Degree Distribution---------------------    
s = sum(sign(abs(W)));  
[lambda_temp,deg_column] = hist(s,deg_column);
lambda = lambda + lambda_temp;    
%--------------------------------------------------------------------------
        
%-------------------Compute Row Degree Distribution------------------------    
s = sum(sign(abs(W')));           
[rho_temp,deg_row] = hist(s,deg_row);   
rho = rho + rho_temp;       
%--------------------------------------------------------------------------

%-------------------------NORMALIZE DISTRIBUTIONS--------------------------
lambda = lambda/N;
rho = rho/(N-K);
%--------------------------------------------------------------------------
