%==========================================================================
%*****************FUNCTION: extract_column(A,index,m)**********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% A: The matrix in its sparse representation. 
% index: The index of the column which we want to retrieve
% m: The number of rows in original matrix.
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% c: The desired column of the matrix in its non-sparse (traditional) format.
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The input to this function is a matrix A in its sparse representation. 
% In other words, in each column of A we have only the index of rows with 
% non-zero elements in the corresponding row and column. 
% The function then gets the index of a column and returns this column in 
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


function c = extract_column(A,index,m)


%-----------------------------Initialization-------------------------------
temp1 = A(:,index);             % Store the column 
c = zeros(m,1);
%--------------------------------------------------------------------------

%---------------------------Extract the Column-----------------------------
for i = 1:length(temp1)
    a = A(i,index);
    c(a) = 1;                   % Here we have assumed the non-zero elements were 1's. This could be extended to the case for other non-zero numbers as well.
end
%--------------------------------------------------------------------------