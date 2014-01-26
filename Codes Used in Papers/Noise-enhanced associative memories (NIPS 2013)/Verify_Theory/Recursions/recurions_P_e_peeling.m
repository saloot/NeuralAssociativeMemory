%==========================================================================
%********************FUNCTION: recurions_P_e_peeling***********************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% av_del_edge: Average number of edges deleted in each round of the peeling algorithm
% lambda: The degree distribution for the pattern nodes of the contracted graph (inter-cluster viewpoint)
% rho: The degree distribution for the cluster nodes of the contracted graph (inter-cluster viewpoint)
% p_e_range: the range of the external probability of +1/-1 noise
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% z: The lower bound on the error correction threshold of the peeling algorithm
% f: The value of the peeling polynomial for each instance of the external error probability (in "p_e_range")
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function numerically solves the recursion corresponding to the 
% genearlized peeling algorithm explained in our March 2013 Progress
% Report.
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function [f,z]=recurions_P_e_peeling(av_del_edge,lambda,rho,p_e_range)


%=============================INITIALIZATION===============================
x = [0:.01:1];
f = [];                                             % The BER bound for each value of the external probability of error
%==========================================================================

%==================COMPUTE THE RECURSION AND FIND BER======================
for j = 1:length(p_e_range)    
    p_e = p_e_range(j);
%     f = [f;x-1+rho_poly(rho,1-p_e*(lambda_poly(lambda,x).^(av_del_edge)))];
%     f = [f;x-1+rho_poly_v2(rho,1-p_e*(lambda_poly_bound_v2(lambda,ones(1,length(lambda)),1-x).^(av_del_edge)))];    
    f = [f;x-1+rho_poly(rho,1-p_e*(lambda_poly_bound_v2(lambda,ones(1,length(lambda)),1-x).^(av_del_edge)))];    
    if (sum(f(j,:)<0)>0)
        z = p_e_range(j-1);
        break;
    end
%     plot(x,f(j,:))
end
%==========================================================================
