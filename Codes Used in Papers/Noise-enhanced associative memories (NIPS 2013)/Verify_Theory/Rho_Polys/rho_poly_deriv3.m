function f = rho_poly_deriv3(rho,x)

f = 0;
for i = 1:length(rho)
    f = f + (i-1)*(i-2)*(i-3)*rho(i)*(x.^(i-4));
end
