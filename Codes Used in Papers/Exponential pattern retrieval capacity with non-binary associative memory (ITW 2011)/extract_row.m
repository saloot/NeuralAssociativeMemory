%==========================================================================
%*******************FUNCTION: extract_row(A,index,n)***********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% A: The matrix in its sparse representation. 
% index: The index of the row which we want to retrieve
% n: The number of columns in original matrix.
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% r: The desired row of the matrix in its non-sparse (traditional) format.
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The input to this function is a matrix A in its sparse representation. 
% In other words, in each column of A we have only the index of columns with 
% non-zero elements in the corresponding row and column. 
% The function then gets the index of a row and returns this row in 
% its traditional representation, i.e. the one in which the non-zero
% elements are shown in their relative positions and other elements are
% equal to zero. 
%--------------------------------------------------------------------------

%-------------------------------FURTHER NOTES------------------------------
% It is assumed that the non-zero elements in the given matrix are only
% 1's. This could be extended to the case for other non-zero values as well
% but A must have a new representation in which not only the index of the
% row with the non-zero element is stored, but also its value is specified.
% -------------------------------------------------------------------------
%==========================================================================
%==========================================================================


function r = extract_row(A,index,n)


%-----------------------------Initialization-------------------------------
temp1 = A(index,:);
r = zeros(1,n);
%--------------------------------------------------------------------------

%---------------------------Extract the Column-----------------------------
for i = 1:length(temp1)
    a = A(index,i);
    if (a ~=0)              % Ensure a is not the all zero row.
        r(a) = 1;           % Here we have assumed the non-zero elements were 1's. This could be extended to the case for other non-zero numbers as well.
    end
end
%--------------------------------------------------------------------------



