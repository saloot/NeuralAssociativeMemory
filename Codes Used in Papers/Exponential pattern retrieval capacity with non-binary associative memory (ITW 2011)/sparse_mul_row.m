%==========================================================================
%*******************FUNCTION: sparse_mul_row(A,x,n)************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% A: The matrix in its sparse representation. 
% x: The vector which we would like to multiply it with A.
% n: The number of columns in original matrix.
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% c: The result of multiplication x*A_1, where A_1 is the original representation of A.
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The input to this function is a matrix A in its sparse representation. 
% In other words, in each row of A we have only the index of columns with 
% non-zero elements in the corresponding row and column. 
% The function then gets a vector and returns the result of multiplying the
% given matrix with vector from the left side.
%--------------------------------------------------------------------------

%-------------------------------FURTHER NOTES------------------------------
% It is assumed that the non-zero elements in the given matrix are only
% 1's. This could be extended to the case for other non-zero values as well
% but A must have a new representation in which not only the index of the
% row with the non-zero element is stored, but also its value is specified.
% -------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function c = sparse_mul_row(A,x,n)

%-----------------------------Initialization-------------------------------
[m,dr] = size(A);
[k,l] = size(x);
c = zeros(1,n);
non_zero_index = [];
%--------------------------------------------------------------------------

%---------------Check the Vailidty of the Input Parameters-----------------
if (m ~=l)
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
    temp = extract_row(A,non_zero_index(i),n);
    if (norm(size(c) - size(temp))~=0)
        error('Something is wrong!');
    end
    c = c + x(non_zero_index(i))*temp;
end
%--------------------------------------------------------------------------