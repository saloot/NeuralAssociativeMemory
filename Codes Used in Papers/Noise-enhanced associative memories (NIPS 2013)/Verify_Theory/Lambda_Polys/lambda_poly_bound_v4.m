% If we consider probability of correcting two errors
function f = lambda_poly_bound_v4(lambda,P_c,rho,x)

P_c_val = zeros(1,4);
for l = 1:length(P_c_val)    
    P_c_val(l) = min(P_c(l,:));
end

f = 0;
for i = 1:length(lambda)
    if (lambda(i) == 0)
        continue
    end
    
    temp = 1-P_c_val(1)*rho_poly_v2(rho,1-x)...
         - P_c_val(2) * x  * rho_poly_deriv(rho,1-x)...
         - P_c_val(3) * x^2* (1/2) * rho_poly_deriv2(rho,1-x)...
         - P_c_val(4) * x^3* (1/6) * rho_poly_deriv3(rho,1-x);
    
    
    f = f + lambda(i)*temp.^i;
%     f = temp;
%     break
end

