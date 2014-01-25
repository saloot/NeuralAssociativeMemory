%==========================================================================
%********************FUNCTION: generator_matrix_planes**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N_horiz: The number of rows in the 2D data
% N_vert: The number of columns in the 2D data
% N_loc_horiz: The number of columns in the random blocks containing non-zero elements of the matrix
% N_loc_vert: The number of rows in the random blocks containing non-zero elements of the matrix
% K_loc: The number of rows in the random blocks containing non-zero elements of the matrix
% deg_column_G_local: The number of non-zero elements in each column of the random blocks
% deg_row_G_local: The number of non-zero elements in each row of the random blocks
% delta_N: The amount of horizontal shift within planes in each iteration
% delta_K: The amount of vertical shift within planes in each iteration
% delta: The amount of horizontal shift in each iteration of the global construction
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% G_tot: the output matrix
% K_tot: the number of rows in the output matrix
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function produces a random binary (0/1) 
% matrix which is used to generate random vectors such that when
% transformed into a 2D array, each group of neighboring "pixels" form a
% subspace. The shape of the matrix is described in our Allerton 2012
% draft. This function corresponds to building the whole matrix G from
% sub-matrices \hat{G}.
% Briefly speaking, the function first generates a random staircase block and
% then the whole matrix is constructed by shifting this block horizontally
% and vertically by the ammounts delta and K_global, resepectively. So the
% non-zero elements in the whole matrix form another staircase shape around 
% the main diagonal of the matrix. 
%--------------------------------------------------------------------------


%==========================================================================
%==========================================================================


function [G_tot,K_tot] = generator_matrix_vert(N_horiz,N_vert,N_loc_horiz,N_loc_vert,K_loc,deg_column_G_local,deg_row_G_local,delta_N,delta_K,delta)

%============================INITIALIZATIONS===============================
itr_no = 1+(N_vert - N_loc_vert)/delta;                                      % The number of times shifting should occur to construct the matrix


%-----------------Check the Validity of the Input Parameters---------------
if (abs(itr_no-round(itr_no))>0)
    error('Invalid input parameters!');
end
%--------------------------------------------------------------------------

%-------------------------Add Necessary Libraries--------------------------
addpath(genpath('/home1/amir/cluster/Common_Library'));                      % Include the library of common functions in the search path
%--------------------------------------------------------------------------

[~,K_global] = generator_matrix_planes(N_horiz,N_loc_horiz,K_loc,deg_column_G_local,deg_row_G_local,delta_N,delta_K);
K_tot = K_global *itr_no;                                                    % The number of rows in the output matrix
G_tot = zeros(K_tot,N_vert*N_horiz);                                         % The output matirx
%==========================================================================

%===============================MAIN LOOP==================================
for itr = 1:itr_no
    start_ind = 1+(itr-1)*delta*N_horiz;                                     % The starting position of the next plane (block)
   
    %----Generate a Random Plane and Position It Properly in the Matrix----
    for j = 1:N_loc_vert
        [G_global,K_global] = generator_matrix_planes(N_horiz,N_loc_horiz,K_loc,deg_column_G_local,deg_row_G_local,delta_N,delta_K);
        G_tot(1+(itr-1)*K_global:itr*K_global,start_ind+(j-1)*N_horiz:start_ind+j*N_horiz-1) = G_global;
    end
    %----------------------------------------------------------------------
    
end
%==========================================================================