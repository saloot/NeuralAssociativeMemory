%==========================================================================
%*********************FUNCTION: bipartite(n,m,l,r)*************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% n: The number of variable (left) nodes.
% m: The number of check (right) nodes.
% l: The degree of each variable (left) node.
% r: The degree of each check (right) node.
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% H: The adjacency matrix of a randomly generated bipartite graph with given parameters
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function gets the specifying parameters of a bipartite graph 
% and generates the adjacency matrix of a random bipartite graph with those
% specifications. 
% The function works by considering sockets attached to both variable and
% check nodes. A permutation of right sockets are then attached to the left
% sockets in order to construct the graph. 
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================



function H = bipartite(n,m,l,r)

%---------------Check the Vailidty of the Input Parameters-----------------
if (n*l ~= m*r)
    error('Invalid input parameters!');
end
%--------------------------------------------------------------------------


%-----------------------------Initialization-------------------------------
H = zeros(m,n);                 % Initialize the matrix.
E = n*l;                        % E is the number of edges.
check_sockets = zeros(E,2);     % Initialize the sockets
%--------------------------------------------------------------------------    

%------------------Assign Right Sockets to Check Nodes---------------------
for i=1:m
    for j = 1:r
        check_sockets((i-1)*r+j,1) = (i-1)*r+j;
        check_sockets((i-1)*r+j,2) = i;
    end
end
%--------------------------------------------------------------------------

%----------Connect Left Sockets to a Permutation of the Right Ones---------
perm = randperm(E);             % Permute the right sockets.

for i = 1:n
    for k = 1:l
        j = perm((i-1)*l+k);    % Connect right socket j to left node i.
        H(check_sockets(j,2),i) = 1;    % Identify the check node corresponding to right socket j and connect to variable node i.
    end
end
%--------------------------------------------------------------------------
    
    