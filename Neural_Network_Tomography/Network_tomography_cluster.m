%---------------------In-loop Initializations--------------------------
    theta = 0.05;   % update threshold
   
    p = 0.3;                        % connection probability
    n = 200;                        % number of neurons
    T = parameter_range(itr);                       % number of sample recordings
    tau = 500000;                        % This is the rate according to which membrane potential drops
    tau2 = 1*tau;
    acc_ave = 0;                    % average number of correctly infered edges
    err_ave = 0;                    % average number of mistakenly infered edges + missed edges
    
    q = .65*(1-1/tau)*theta/p;                    % stimulus firing probability    
    avg_no_edge = round(p*n);       % estimated number of edges in the graph
    no_edges_fire = ceil(theta*n);  % estimated number of edges required for a neuron to fire
    no_averaging_itrs = 30;         % number of times we perform each simulation for the sake of averging