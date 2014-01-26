function f = rho_poly(rho,x)

f = 0;
for i = 1:length(rho)
    f = f + rho(i)*(x.^(i-1));
end
