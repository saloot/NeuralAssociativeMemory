function G = erdos_reny_inhibitory(m,n_exc,n_inh,p)

%-------------------------Vectorized Version-------------------------------
G = (rand(m,n_exc)<p);
G = [G,-(rand(m,n_inh)<p)];
%--------------------------------------------------------------------------

if (m==n_exc+n_inh)
    for i = 1:min(m,n_exc+n_inh)
        G(i,i) = 0;
    end
end