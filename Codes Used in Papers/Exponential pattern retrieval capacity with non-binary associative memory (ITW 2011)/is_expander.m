%==========================================================================
%*********************FUNCTION: is_expander(H,beta)************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% H: The adjacency matrix of a bipartite graph
% beta: The expansion coefficient
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% S: The size of the neighborhood for which expansion prperty holds. 
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the adjacency matrix of an expander graph and an
% expansion coefficient. Then it checks various neghborhood sizes of
% variable nodes and see if the expansion property holds. This process is
% repeated with increasing neghborhood sizes until at some point the
% expansion property breaks. The size of this neighborhood minus one is
% returned as the output.
%--------------------------------------------------------------------------

%-------------------------------FURTHER NOTES------------------------------
% The method used here is pure brute force and for that it takes a long
% time to run. More intelligent ways will result in faster checks. 
% -------------------------------------------------------------------------
%==========================================================================
%==========================================================================

function S = is_expander(H,beta)


%-----------------------------Initialization-------------------------------
[m,n] = size(H);                    % Extract the size of the matrix. 
max_itr = 1000000;                  
%--------------------------------------------------------------------------

for S = 2:n                         % Repeat for increasing neighborhood sizes of left nodes.
    expander_flag = 0;              % This flag will be triggered when the expansion property is not hold.
    
    for itr = 1:max_itr        
        
        %---Count the Number of Outgoing Edges for a Random Neighborhood---
        p = randperm(n);            % Consider a permutation of left nodes
        summ = zeros(m,1);  
        
        optimal_neighborhood = 0;
        for i = 1:S
            summ = summ + abs(H(:,p(i)));   
            optimal_neighborhood = optimal_neighborhood + sum(abs(H(:,p(i))));  % Count the number of outgoing edges
        end
        %------------------------------------------------------------------
        
        neighbor_count = sum(sign(summ));   % Count the actual number of neighbors
        
        
        %-----------------Check the Expansion Property---------------------
        if (neighbor_count < beta * optimal_neighborhood)
            expander_flag = 1;
            break;
        end
        %------------------------------------------------------------------
    end
    
    
    if (expander_flag == 1)         % If the expansion property is violated, then terminate. 
        break;
    end
end
S = S-1;                % return the neighborhood size for which the expansion property holds. 
