%==========================================================================
%***********FUNCTION: random_matrix_sparse(m,n,deg_column,deg_row)*********
%==========================================================================

%--------------------------------INPUTS------------------------------------
% m: The number of rows in the desired matrix
% n: The number of columns in desired matrix
% deg_column: The average column degree of each column
% deg_row: The average row degree of each row
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% A_column: The sparse bipartite matrix in its column format
% A_row: The sparse bipartite matrix in its row format
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the size and specifications of a regular bipartite 
% graphs and retunrs the adjacency matrix of such a graph in two 
% representations: in the column representation, only the position of the 
% non-zero elements in each column is stored. In the row representation, 
% the position of the non-zero elements in each row is stored.
%--------------------------------------------------------------------------

%-------------------------------FURTHER NOTES------------------------------
% It is assumed that the non-zero elements in the given matrix are only
% 1's. This could be extended to the case for other non-zero values as well
% -------------------------------------------------------------------------
%==========================================================================
%==========================================================================



function [A_column,A_row] = random_matrix_sparse(m,n,deg_column,deg_row)

%---------------Check the Validity of the Input Parameters-----------------
if (m*deg_row ~= n*deg_column)
    error('Invalid input parameters');
end
%--------------------------------------------------------------------------

%-----------------------------Initialization-------------------------------
A_column = zeros(deg_column,n);         % Initialize the column representation
A_row = zeros(m,deg_row);               % Initialize the row representation
E = n*deg_column;                       % Compute the number of edges
%--------------------------------------------------------------------------


%-------------------Generate a Random Bipartite Graph----------------------
edge = randperm(E);                     % Permute the check nodes sockets
for i = 1:E 
    var_node = 1+floor((i-1)/deg_column);       % Find the corresponding variable node to the current socket
    var_socket = 1+mod(i,deg_column);           % Find the corresponding socket number in the variable node.
    check_node = 1+floor((edge(i)-1)/deg_row);  % Find the corresponding check node to the current socket
    check_socket = 1+mod(edge(i),deg_row);      % Find the corresponding socket number in the check node.
    A_column(var_socket,var_node) = check_node; % Connect the variable socket to the check socket
    A_row(check_node,check_socket) = var_node;  % Connect the variable socket to the check socket
end
%--------------------------------------------------------------------------


%-----------------------Remove Similar Columns-----------------------------
for j = 1:n
    for i = j+1:n
        col1 = extract_column(A_column,i,m);    % Extract the first columns
        col2 = extract_column(A_column,j,m);    % Extract the second columns
        
        if (col1 == col2)                       % Check for columns with equal values.
            
            %-------------Shift all the columns one to the left------------ 
            for k = i:n-1                       
                A_column(:,k) = A_column(:,k+1);
            end
            %--------------------------------------------------------------
            
            %--------Update the positions in the row representation--------
            for ii = 1:m
                for jj = 1:deg_row
                    if (A_row(ii,jj) >i)
                        A_row(ii,jj) = A_row(ii,jj)-1;
                    end
                    if (A_row(ii,jj) ==i)
                        A_row(ii,jj) = 0;
                    end
                end
            end
            %--------------------------------------------------------------
                    
            n = n-1;                            % Update the number of columns.
        end
    end
end
%--------------------------------------------------------------------------

%-----------------Remove the Redundant Columns completely------------------
AA = A_column(:,1:n);
A_column = AA;
%--------------------------------------------------------------------------