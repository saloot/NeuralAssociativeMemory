% If we consider probability of correcting two errors
function f = lambda_poly_bound_v3(lambda,P_correct_cluster,P_correct_cluster2,rho,x)

P_c = sort(P_correct_cluster,'ascend');
P_c_2 = sort(P_correct_cluster2,'ascend');
f = 0;
for i = 1:length(lambda)
    if (lambda(i) == 0)
        continue
    end
    temp = 1;
    for j = 1:i
        temp = temp*(1-P_c(j)*rho_poly(rho,1-x)-P_c_2(j)*x*rho_poly_deriv(rho,1-x));
    end
    
    f = f + lambda(i)*temp;
%     f = temp;
%     break
end
