% Calculating the potential function based on the results of the paper 
% "A Simple Proof of Threshold Saturation for Coupled Scalar Recursions"

%============================INITIALIZATIONS===============================

%--------------------Initialize the Degree Distributions-------------------
lambda =[0.0011,0.0032,0.0043,0.0722,0,0.0054,0,0.0841,0.0032,0,0,0.098,0,0,0,0.7284];
% lambda = [0 1];
rho = [zeros(1,63),1];
%--------------------------------------------------------------------------

%--------------------Reverse the Order of Degree Distributions-------------
lambda_r = [fliplr(lambda),0];
rho_r = [fliplr(rho)];
%--------------------------------------------------------------------------

%-------------------------Other Initializations----------------------------
x = [0:.001:1];
e = 2;                                          % e is the number of errors the system can correct
p_e_range = [0:.001:1];
Omega_min = zeros(1,length(p_e_range));
Delta_E = zeros(1,length(p_e_range));
U = zeros(length(p_e_range),length(x));
p_e_star_flag = 0;
p_e_dagger_flag = 0;
p_e_dagger = 1;
p_e_star = 1;
%--------------------------------------------------------------------------

%==========================================================================


%=========================CALCULATE f(x) and g(x)==========================
if (e == 2)
    g = 1-polyval(rho_r,1-x)-x.*polyval(polyder(rho_r),1-x);
    g_0 = 1-polyval(rho_r,1);
elseif (e == 1)
    g = 1-polyval(rho_r,1-x);
    g_0 = 1-polyval(rho_r,1);
else
    error('Unsupported e!');
end
%==========================================================================


for ij = 1:length(p_e_range)
    p_e = p_e_range(ij);
    %========================CALCULATE epsilon_s^*=========================
    %--------------According to Definition 4 of the Paper------------------
    if (p_e_dagger_flag == 0)
        U_der = x-p_e*polyval(lambda_r,g);
        if (min(U_der) < 0)        
            p_e_dagger = p_e-p_e_range(2);
            p_e_dagger_ind = ij;
            p_e_dagger_flag = 1;
        end
    end
    %======================================================================


    %========================CALCULATE u_epsilon===========================
    %---------------According to Formula (3) of the Pap--------------------
    U_der = x-p_e*polyval(lambda_r,g);
    for i = 1:length(x)    
        if (U_der(i) < 0)        
            u_epsilon = x(i-1);
            u_e_ind = i;
            break
        end
    end
    %======================================================================


    %==========================CALCULATE K_f_g=============================
    %-------------------According to Lemma 5 of the Paper------------------
    if (e == 2)
    % f = p_e * lambda(x)
    %  g = 1-rho(1-x) + x*rho_der(1-x)
        f_deriv = p_e*polyval(polyder(lambda_r),x);
        g_deriv = x.*polyval(polyder(polyder(rho_r)),1-x);
        g_deriv_deriv = polyval(polyder(polyder(rho_r)),1-x) - x.*polyval(polyder(polyder(polyder(rho_r))),1-x);

        norm_inf_f_deriv = max(abs(f_deriv));
        norm_inf_g_deriv = max(abs(g_deriv));
        norm_inf_g_deriv_deriv = max(abs(g_deriv_deriv));

        K_f_g = norm_inf_g_deriv + (norm_inf_g_deriv^2 * norm_inf_f_deriv) + norm_inf_g_deriv_deriv;
    elseif (e == 1)
        f_deriv = p_e*polyval(polyder(lambda_r),x);
        g_deriv = polyval(polyder(rho_r),1-x);
        g_deriv_deriv = -polyval(polyder(polyder(rho_r)),1-x);

        norm_inf_f_deriv = max(abs(f_deriv));
        norm_inf_g_deriv = max(abs(g_deriv));
        norm_inf_g_deriv_deriv = max(abs(g_deriv_deriv));

        K_f_g = norm_inf_g_deriv + (norm_inf_g_deriv^2 * norm_inf_f_deriv) + norm_inf_g_deriv_deriv;
    end
    %======================================================================

    
    %=====================CALCULATE THE POTENTIAL==========================
    if (e == 2)            
        G_int = x + polyval(polyint(rho_r),1-x) + (x.*polyval(rho_r,1-x) + polyval(polyint(rho_r),1-x))-(0 + polyval(polyint(rho_r),1) + (0.*polyval(rho_r,1) + polyval(polyint(rho_r),1)));            
        F_int = p_e*polyval(polyint(lambda_r),g)-p_e*polyval(polyint(lambda_r),g_0);        
    elseif (e == 1)            
        G_int = x + polyval(polyint(rho_r),1-x) -(0 + polyval(polyint(rho_r),1) );                   
        F_int = p_e*polyval(polyint(lambda_r),g)-p_e*polyval(polyint(lambda_r),g_0);
    end    
    U(ij,:) = x.*g - G_int - F_int;
    %======================================================================
        
    if (p_e > p_e_dagger)
        %========================CALCULATE Delta_E=========================        
        U_short = U(ij,u_e_ind:end);
        Delta_E(ij) = min(U_short);
        Omega_min(ij)= K_f_g/Delta_E(ij);
        %==================================================================

        %======================CALCULATE epsilon^*=========================
        %------------According to Definition 4 of the Paper----------------
        if (p_e_star_flag == 0)                
            if (min(U(ij,:)) < 0)        
                p_e_star = p_e-p_e_range(2);
                p_e_star_flag = 1;
%                 break;
            end
        end   
        %==================================================================
    end
end

plot(x,U(ij,:))
% plot(p_e_range(1:ij-1), Delta_E(1:ij-1));
