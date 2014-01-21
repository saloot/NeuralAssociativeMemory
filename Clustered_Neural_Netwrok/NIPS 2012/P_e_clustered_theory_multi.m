%==========================================================================
%*****************FUNCTION: P_e_clustered_theory_multi*********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% L: The number of clusters if we have multiple levels (L = 1 for single level)
% lambda_in: The input pattern nodes degree distribution from a cluster point of view in a clustered neural associative memory
% deg_column: The input pattern nodes degrees from a cluster point of view in a clustered neural associative memory
% rho_in: The input cluster degree distribution from a cluster point of view in a clustered neural associative memory
% deg_row: The input cluster degrees from a cluster point of view in a clustered neural associative memory
% error_bits: Number of initial noisy nodes. It can be a vector
% e: The number of initial errors each cluster can correct with probability very close to 1
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% Pe: The theoretical recall probability of error for each initial number of erroneous nodes. It can be a vector (one element for each error_bits entry).
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function calculates the theoretical probability of recall error in a
% clustered neural associative memory assuming each cluster is capable of
% correcting e initial errors.
% The calculation method is in fact very similar to finding the threshold
% for erasure codes.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function [Pe] = P_e_clustered_theory_multi(N,K,L,lambda_in,deg_column,rho_in,deg_row,error_bits,e,gamma_BFS,try_max,lambda_within,deg_column_within)

%-----------------------------Initialization-------------------------------
deg_ave = sum(deg_column_within.*lambda_within);
Pe = [];
load(['/scratch/amir/Clustered_Neural/Initialization_Files/N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'/clustered_neural_parameters_v1_N_',num2str(N),'_K_',num2str(K),'_L_',num2str(L),'_index_',num2str(1),'.mat']);           
%--------------------------------------------------------------------------

%---------------------------------Main Loop--------------------------------
for k = 1:length(error_bits)   
    q = error_bits(k)/N_tot;
    q0 = q;
    
    end_flag = 0;
    try_itr = 0;
    while (end_flag == 0)
        try_itr = try_itr + 1;
        
        %---------------------Calculate pai(t)-----------------------------
        pai = 0;        
        for j = 1:length(deg_row) 
            d = deg_row(j);
            
            P1 = 0;
            for i = 0:e
                if ( i==0)
                    P1 = P1 + nchoosek(d,i)*(q^i)*( (1-q)^(d-i));
                elseif (i ==1)
                    P1 = P1 + nchoosek(d,i)*(q^i)*( (1-q)^(d-i)) * (1-lambda_poly(deg_ave/(N-K),lambda_within));
                else
                    error('Temp error!');
                end
            end
            P1 = 1-P1;
            
            pai = pai + rho_in(j)*P1;
        end
        %------------------------------------------------------------------
        
        %-------------------Calculate q(t+1)-------------------------------
        q_temp = 0;
        for j = 1:length(deg_column) 
            d = deg_column(j);
            
            P2 = 0;
            for i = ceil(gamma_BFS*d):d
                P2 = P2 + nchoosek(d,i) * (pai^i) *( (1-pai)^(d-i) );
            end
            q_temp = q_temp + lambda_in(j)*P2;        
        end                    
        q = q_temp*q0;
        %------------------------------------------------------------------    
        
        %------------------Check the Stopping Condition--------------------
        if (( q==0) || (q==1)||(try_itr>try_max))
            end_flag = 1;
        end
        %------------------------------------------------------------------
        
    end
    Pe = [Pe,1-(1-q)^N];
end
%--------------------------------------------------------------------------