function f = rho_poly_deriv2(rho,x)

f = 0;
for i = length(rho):-1:1
    if (rho(i) == 0)
        continue;
    end
    f = f + (i-1)*(i-2)*rho(i)*(x.^(i-3));    
end
