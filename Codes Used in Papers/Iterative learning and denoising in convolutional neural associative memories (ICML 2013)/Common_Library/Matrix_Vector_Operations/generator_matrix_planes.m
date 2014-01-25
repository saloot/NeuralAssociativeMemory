%==========================================================================
%********************FUNCTION: generator_matrix_planes**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of columns in the output matrix
% N_loc: The number of columns in the random blocks containing non-zero elements of the matrix
% K_loc: The number of rows in the random blocks containing non-zero elements of the matrix
% deg_column_G_local: The number of non-zero elements in each column of the random blocks
% deg_row_G_local: The number of non-zero elements in each row of the random blocks
% delta_N: The amount of horizontal shift in each iteration
% delta_K: The amount of vertical shift in each iteration
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% G_lobal: the output matrix
% K: the number of rows in the output matrix
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function is part of a procedure to produces a random binary (0/1) 
% matrix which is used to generate random vectors such that when
% transformed into a 2D array, each group of neighboring "pixels" form a
% subspace. The shape of the matrix is described in our Allerton 2012
% draft. Briefly speaking, the function first generates a random block and
% then the whole matrix is constructed by shifting this block horizontally
% and vertically by the ammounts delta_N and delta_K, resepectively. So the
% non-zero elements in the whole matrix form a staircase shape around the
% main diagonal of the matrix. 
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================

function [G_global,K] = generator_matrix_planes(N,N_loc,K_loc,deg_column_G_local,deg_row_G_local,delta_N,delta_K)

%============================INITIALIZATIONS===============================
itr_no = 1+(N - N_loc)/delta_N;                                         % The number of times the horizontal and vertical shifting should occur
K = K_loc + delta_K*(itr_no-1);                                         % This is the number of rows in the output matrix
G_global = zeros(K,N);                                                  % This is the output matrix

%-----------------Check the Validity of the Input Parameters---------------
if (abs(itr_no-round(itr_no))>0)
    error('Invalid input parameters!');
end
%--------------------------------------------------------------------------

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                 % Include the library of common functions in the search path
%--------------------------------------------------------------------------

%==========================================================================


%==============================MAIN LOOP===================================
% G_local = bipartite(N_loc,K_loc,deg_column_G_local,deg_row_G_local);              
for itr = 1:itr_no    
    G_local = bipartite(N_loc,K_loc,deg_column_G_local,deg_row_G_local);                                  % Generate the random block
    G_global(1+(itr-1)*delta_K:K_loc+(itr-1)*delta_K,1+(itr-1)*delta_N:N_loc+(itr-1)*delta_N) = G_local;  % Put the block in the proper position in the whole matrix
end
%==========================================================================