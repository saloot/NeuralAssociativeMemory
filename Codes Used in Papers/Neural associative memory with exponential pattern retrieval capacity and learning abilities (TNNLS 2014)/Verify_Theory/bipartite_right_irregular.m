%==========================================================================
%************FUNCTION: bipartite_right_irregular(n,m,d,lambda)*************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% n: The number of variable (left) nodes.
% m: The number of check (right) nodes.
% d: The degree distribution of variable nodes.
% lambda: The distribution of degrees for variable nodes.
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% H: The adjacency matrix of a randomly generated bipartite graph with given parameters
% neighborsize: The total number of check nodes connected to the variable nodes in each round of the construction
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function gets the specifying parameters of a bipartite graph 
% and generates the adjacency matrix of a random bipartite graph with those
% specifications. 
% The function works by picking a variable node in each iteration and
% determining its degree according to the given degree distribution. Based
% on the determined degree, denoted by d, the picked variable node is 
% connected to d constraint nodes uniformly at random.
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================



function [H,neighborsize] = bipartite_right_irregular(n,m,d,lambda)


%-----------------------------Initialization-------------------------------
H = zeros(m,n);                 % Initialize the matrix.
neighborsize = zeros(1,n);      % The size of the neighborhood in each iteration
%--------------------------------------------------------------------------    

%------------------Assign Right Sockets to Check Nodes---------------------
for i=1:n

    %-------------Determine the Degree of the Random Variable--------------
    p = rand;
    F = lambda(1);
    
    for k = 2:length(lambda)
        if (p <= F)           
            break;
        else
            F = F + lambda(k);
        end
    end    
    %----------------------------------------------------------------------
        
    %----------Connect the Variable Nodes to the Constraint Nodes----------
    e = randperm(m);    
    for j = 1:d(k-1)
        H(e(j),i) = 1;
    end    
    HH = H';
    neighborsize(i) = sum(sign(sum(HH)));
    %----------------------------------------------------------------------
    
end
%--------------------------------------------------------------------------