%==========================================================================
%**********************FUNCTION: recurions_P_e_v2**************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% P_c: The array containing the probabilities of correcting one (and possibly more) error(s) by a cluster on average
% lambda: The degree distribution for the pattern nodes of the contracted graph (inter-cluster viewpoint)
% rho: The degree distribution for the cluster nodes of the contracted graph (inter-cluster viewpoint)
% p_e_range: the range of the external probability of +1/-1 noise
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% f: The upper bound on BER of the inter-cluster error correcing algorithm
%--------------------------------------------------------------------------


%--------------------------FUNCTION DESCRIPTION----------------------------
% This function numerically solves the recursion corresponding to the 
% theoretical upper bound on the probability of error in a clustered neural 
% networks built out of faulty neurons. The bound is a looser version of 
% the one we explained in our ISIT 2013 paper. However, it is theoretically
% more sound. 
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


function f = recurions_P_e(P_c,lambda,rho,p_e_range,version)


%=============================INITIALIZATION===============================
max_itr = 1000;                                                             % The maximum number of times the recursion is executed
f = zeros(1,length(p_e_range));                                             % The BER bound for each value of the external probability of error
%==========================================================================

%==================COMPUTE THE RECURSION AND FIND BER======================
for j = 1:length(p_e_range)
    p_e = p_e_range(j);
    z = p_e;
    for itr = 1:max_itr
        
        z_old = z;
        switch version
            case 1                
                z = p_e*(lambda_poly(lambda,1-P_c(1)*rho_poly(rho,1-z)));
            case 2
                z = p_e*(lambda_poly_bound_v2(lambda,P_c,rho_poly_v2(rho,1-z)));        
            case 3                                        
                z = p_e*(lambda_poly_bound_v3(lambda,P_c,rho,z));
            case 4
                z = p_e*lambda_poly_bound_v4(lambda,P_c,rho,z);                
            otherwise
                error('Invalid version!')
        end
        
        if (norm(z-z_old)/z_old < 0.0001)
            break;
        end        
        if (( norm(z-1)<.000001) || (z<1e-8))
            break;
        end
    end
    f(j) = z;
end
%==========================================================================
