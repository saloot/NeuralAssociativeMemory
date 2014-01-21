lambda_in = [0 0 1];
rho_in = [0 0 0 0 0 1];
try_max = 100;
for x0 = 0:0.001:.9
    x_ave = x0;
    x_ave2 = x0;
    end_flag = 0;
    try_itr = 0;
    x = x0;
    while (end_flag == 0)
        try_itr = try_itr + 1;
%         P_e_av = 0;
%         for j = 1:length(deg_column)
%             P_e_av = 
        x = x0*lambda_poly(1-rho_poly(1-x,rho_in)-x*rho_poly_der(1-x,rho_in)+rho_poly(1-x,rho_in)*pe_max+pe_2_max*x*rho_poly_der(1-x,rho_in),lambda_in);
%         x_ave = x0*lambda_poly(1-rho_poly(1-x,rho_in)-x*rho_poly_der(1-x,rho_in)+rho_poly(1-x,rho_in)*pe_ave+pe_2_ave*x*rho_poly_der(1-x,rho_in),lambda_in);
        x_ave = x0*lambda_poly_v2(1-(1-pe_ave)*rho_poly(1-x_ave,rho_in),lambda_in);
        x_ave2 = x0*lambda_poly_v2(1-rho_poly(1-x_ave2,rho_in),lambda_in);
        if (( x_ave2==0) || (x_ave2==1)||(try_itr>try_max))
            end_flag = 1;
        end
    end
    if (x_ave2 > .001)
            break;
    end
end
