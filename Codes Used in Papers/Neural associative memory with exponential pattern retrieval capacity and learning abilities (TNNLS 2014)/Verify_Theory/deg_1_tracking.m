

x = [0:.01:1];
p1_0 = (1-x);
f_rho = rho_poly(rho,x);
f1 = x.*lambda_poly(lambda,1-x.*f_rho);
f2 = 0.5*p1_0.*lambda_poly(lambda,sqrt(ones(1,length(x))-f_rho.^2));
plot(x,f1-f2)