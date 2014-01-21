%==========================================================================
%*****************FUNCTION: P_e_clustered_theory_single********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% lambda_in: The input pattern nodes degree distribution from a cluster point of view in a clustered neural associative memory
% rho_in: The input cluster degree distribution from a cluster point of view in a clustered neural associative memory
% error_bits: Number of initial noisy nodes. It can be a vector
% try_max: The maximum number of iterations in outside-cluster recall algorithm
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe: The theoretical recall probability of error for each initial number of erroneous nodes. It can be a vector (one element for each error_bits entry).
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates the theoretical probability of recall error in a
% clustered neural associative memory assuming each cluster is capable of
% correcting only one bit of error.
% The calculation method is in fact very similar to finding the threshold
% for erasure codes.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================

function [Pe,Pb] = P_e_clustered_theory_single(N,K,L,lambda_in,rho_in,error_bits,try_max)

%-----------------------------Initialization-------------------------------
Pe = [];
Pb = [];
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           
%--------------------------------------------------------------------------

%----------------------------Main Loop-------------------------------------
for k = 1:length(error_bits)   
    x = error_bits(k)/N_tot;                                    % This is the input probability of error
    x0 = x;
    end_flag = 0;
    try_itr = 0;
    while (end_flag == 0)
        try_itr = try_itr + 1;
        x = x0*lambda_poly(1-rho_poly(1-x,rho_in),lambda_in);
        if (( x==0) || (x==1)||(try_itr>try_max))
            end_flag = 1;
        end
    end
    x = max(x,0);
    Pe = [Pe,1-(1-x)^N];
    Pb = [Pb,x];
end
%--------------------------------------------------------------------------
% 111