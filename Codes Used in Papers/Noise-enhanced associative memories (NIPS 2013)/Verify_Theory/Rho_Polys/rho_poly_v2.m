function f = rho_poly_v2(rho,x)

f = 0;
for i = length(rho):-1:1
    if (rho(i) == 0)
        continue
    end
    f = (x.^(i-1));
    break
end
