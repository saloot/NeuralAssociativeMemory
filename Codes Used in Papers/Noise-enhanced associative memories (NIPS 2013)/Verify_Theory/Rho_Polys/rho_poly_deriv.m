function f = rho_poly_deriv(rho,x)

f = 0;
for i = length(rho):-1:1
    if (rho(i) == 0)
        continue
    end
    f = f + (i-1)*rho(i)*(x.^(i-2));
end
