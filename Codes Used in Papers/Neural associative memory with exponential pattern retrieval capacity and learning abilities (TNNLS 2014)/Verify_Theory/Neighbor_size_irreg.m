%==========================================================================
%**************FUNCTION: Neighbor_size_irreg(n,m,d,lambda)*****************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% n: The number of variable (left) nodes.
% m: The number of check (right) nodes.
% d: The degree distribution of variable nodes.
% lambda: The distribution of degrees for variable nodes.
% -------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% neighborsize_simul: The average neighborhood size obtained from simulations
% neighborsize_theory: The average neighborhood size derived theoretically
% -------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% The following function verifies the formula derived in the progress
% report of 17-30 March 2012 for the average number of constraint nodes
% connected to a set of variable nodes when constructing a bipartite graph
% with a given degree. 
%--------------------------------------------------------------------------
%==========================================================================
%==========================================================================

function [neighborsize_simul,neighborsize_theory,MSE] = Neighbor_size_irreg(n,m,d,lambda)

%----------------Check the Vailidity of the Input Parameters---------------
if (norm(sum(lambda) - 1)>0.001)
    error('Probability distribution must add up to 1');
end
if (length(lambda) ~= length(d))
    error('Invalid input!');
end
%--------------------------------------------------------------------------

%-----------------------------Initialization-------------------------------
deg_ave = sum(d.*lambda);               % The average degree of variable nodes
no_of_ensembles = 200;                  % The number of randomly generated irregular bipartite graphs.
neighborsize_theory = [deg_ave];
neighborsize_simul = zeros(1,n);
a=clock;                                % Initialize the seed for random number generation with the clock value.
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a))); 
%--------------------------------------------------------------------------

%---------------Calculate the Simulated Neighborhood Size------------------
for k = 1:no_of_ensembles
    [~,neigh] = bipartite_right_irregular(n,m,d,lambda);
    neighborsize_simul = neighborsize_simul + neigh;
end
neighborsize_simul = neighborsize_simul/k;
%--------------------------------------------------------------------------

%---------------Calculate the Theretical Neighborhood Size-----------------
for t = 1:n-1
      neighborsize_theory = [neighborsize_theory,m*(1-((1-deg_ave/m)^(t+1)))];
end
%--------------------------------------------------------------------------

%--------------------------Display the Results-----------------------------
MSE = norm(neighborsize_simul - neighborsize_theory)                % Calculate the MSE between theory and practice
figure
plot (neighborsize_theory,'g');
hold on
plot(neighborsize_simul,'b')
legend('Theory','Simulation');
%--------------------------------------------------------------------------