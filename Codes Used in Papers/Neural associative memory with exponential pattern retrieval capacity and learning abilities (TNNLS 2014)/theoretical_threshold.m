function x_thr = theoretical_threshold(n,m,lambda,deg_column,rho,deg_row)

x_thr = 0;
x = [0:.01:1];
p1_0 = (1-x);
f_rho = rho_poly(rho,x);
f1 = x.*lambda_poly(lambda,1-x.*f_rho);



% deg_ave = sum(d.*lambda);

P2_tot = zeros(1,length(x));                
for j = 1:length(deg_column)    
    P2 = zeros(1,length(x));                            
    dp = deg_column(j);            
    alpha = rho_poly(rho,x);    
    for i = ceil(0.5*dp):dp    
        P2 = P2 + nchoosek(dp,i) * ((1+alpha).^i) .* ((1-alpha).^(dp-i)) ;          
    end    
    P2_tot = P2_tot + lambda(j)*P2*(.5)^dp;    
end
g = p1_0.*P2_tot;



a = ((f1-g)<0);
if (sum(a) == 0)
    x_thr = 1;
    display('It is not possible to find a threshold');
else
    for i = length(a)-1:-1:1
        if (a(i)==0)
            x_thr = x(i+1);
            break;
        end
    end
end