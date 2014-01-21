%==========================================================================
%*****************FUNCTION: P_e_clustered_theory_single********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% cluster_size: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% no_clusters: The number of clusters if we have multiple levels (no_clusters = 1 for single level)
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

function [Pe,Pe_ave,Pb,Pb_ave] = P_e_clustered_theory_main_v2(cluster_size,N_tot,lambda_in,rho_in,error_bits,try_max,pe_max,pe_ave)

%-----------------------------Initialization-------------------------------
Pe = [];
Pe_ave = [];
Pb = [];
Pb_ave = [];
%--------------------------------------------------------------------------


%----------------------------Main Loop-------------------------------------
for k = 1:length(error_bits)   
    x = error_bits(k)/N_tot;                                    % This is the input probability of error
    x_ave = x;
    x0 = x;
    end_flag = 0;
    try_itr = 0;
    while (end_flag == 0)
        try_itr = try_itr + 1;
%         P_e_av = 0;
%         for j = 1:length(deg_column)
%             P_e_av = 
        x = x0*lambda_poly(1-(1-pe_max)*rho_poly(1-x,rho_in),lambda_in);
        x_ave = x0*lambda_poly(1-(1-pe_max)*rho_poly(1-x,rho_in)+rho_poly(1-x,rho_in)*pe_ave,lambda_in);
        if (( x==0) || (x==1)||(try_itr>try_max))
            end_flag = 1;
        end
    end
    Pe = [Pe,1-(1-x)^cluster_size];
    Pe_ave = [Pe_ave,1-(1-x_ave)^cluster_size];
    Pb = [Pb,x];
    Pb_ave = [Pb_ave,x_ave];
end
%--------------------------------------------------------------------------