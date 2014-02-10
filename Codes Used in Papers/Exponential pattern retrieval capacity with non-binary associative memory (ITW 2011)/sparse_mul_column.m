%==========================================================================
%*****************FUNCTION: sparse_mul_column(A,x,m)***********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% A: The matrix in its sparse representation. 
% x: The vector which we would like to multiply it with A.
% m: The number of rows in original matrix.
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% c: The result of multiplication A_1*x, where A_1 is the original representation of A.
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The input to this function is a matrix A in its sparse representation. 
% In other words, in each column of A we have only the index of rows with 
% non-zero elements in the corresponding row and column. 
% The function then gets a vector and returns the result of multiplying the
% given matrix with vector from the right side.
%--------------------------------------------------------------------------

%-------------------------------FURTHER NOTES------------------------------
% It is assumed that the non-zero elements in the given matrix are only
% 1's. This could be extended to the case for other non-zero values as well
% but A must have a new representation in which not only the index of the
% row with the non-zero element is stored, but also its value is specified.
% -------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function c = sparse_mul_column(A,x,m)

%-----------------------------Initialization-------------------------------
[dc,n] = size(A);
[k,l] = size(x);
c = zeros(m,1);
non_zero_index = [];
%--------------------------------------------------------------------------

%---------------Check the Vailidty of the Input Parameters-----------------
if (n ~=k)
    error('Invalid input!');
end
%--------------------------------------------------------------------------

%---------------------Gather Non-zero positions in x-----------------------
for i = 1:length(x)
    if (x(i)~=0)
        non_zero_index = [non_zero_index,i];
    end
end
%--------------------------------------------------------------------------

%----------------Calculate the Column Multiplication-----------------------
for i = 1:length(non_zero_index)
    temp = extract_column(A,non_zero_index(i),m);
    c = c + x(non_zero_index(i))*temp;
end
%--------------------------------------------------------------------------