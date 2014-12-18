%==========================================================================
%***********************FUNCTION: neural_journal_v1******************************
%==========================================================================

%--------------------------------INPUTS------------------------------------
% N: The number of pattern nodes in the graph
% K: The dimension of the subspace of the pattern nodes
% N_const_per_job: Number of constraints learned by each submitted job to the cluster
% alpha0: The step size in the learning algorithm
% beta0: The sparsity penalty coefficient in the learning algorithm
% theta0: The sparsity threshold in the learning algorithm
% error_bits: Number of initial noisy nodes. It can be a vector
% max_noise_amp: The maximum amplitude of noise during the recall process
% no_simulated_instances: The number of patterns considered in the recall phase.
% simul_instances_per_job: Number of simulated instances in the recall phase per submitted job to the cluster
% gamma_BFO: The update threshold in the original bit-flipping recall algorithm 
% gamma_BFS: The update threshold in the simplified (yet better) bit-flipping recall algorithm 
%--------------------------------------------------------------------------

%--------------------------------OUTPUTS-----------------------------------
% None
%--------------------------------------------------------------------------

%--------------------------FUNCTION DESCRIPTION----------------------------
% This function gets the specification of a neural network and performs the
% learning and recall phases of a bipartite neural associative memory. More 
% specifically, the function first learns a pre-specified number of 
% constraints. Then the function goes directly to the recall phase. 
% Each step is performed by dividing the jobs in parallel and submitting
% them to the cluster. 
% The learning algorithm is based on the approach mentioned in our journal
% paper. The recall process considers three methods: The winner-take-all,
% the original bit flipping (mentioned in our ITW 2011 paper) and the
% simplified bit-flipping (proposed in ISIT 2012).
%--------------------------------------------------------------------------

%==========================================================================
%==========================================================================


%%
%=============================INITIALIZATION===============================
if (~exist('initialization_done_by_compare','var'))    
    N = 300;                                % N is the number of neurons in network.
    K = 150;                                % K is the number of message bits.
    z_max = 1;                              % This is the maximum value of message bits.
    z_min = 0;                              % This is the minimum value of message bits.
    index_max = 50;                         % This is the maximum number of random scenarios generated for simulation
    random_flag = 1;                        % Determines if patterns are drawn from a subspace or generated randomly
    no_of_patterns = 200;                   % The number of patterns in the dataset
    no_simulated_instances = 500;         % The number of patterns that are going to be denoised during the recall phase
    max_noise_amp = 1;                      % Maximum value of integer-valued noise added to each bit  
    err_bits_range = [0:10];                % The number of bits that will be corrupted initially for the recall phase
    addpath(genpath('../../Common_Library'));                          % Include the library of common functions in the search path
    plot_flag = 1;                          % Plot the graphs    
    a=clock;                                % Initialize the seed for random number generation with the clock value.
    RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*a)));    
end
if (random_flag)
    mkdir('../../Recall_Results/Lee',['N_',num2str(N),'_Random'])
else
    mkdir('../../Recall_Results/Lee',['N_',num2str(N),'_K_',num2str(K)])
end
initialization_done_by_Lee = 1;
%==========================================================================       


%%
for random_flag = 0:1
%==========================PREPARE THE DATASET==============================
    dataset_learn = zeros(no_of_patterns,N);
    
    for l = 1:index_max
        display(['Doing the phase for graph ',num2str(l),' of the ensemble'])
        %------------------Pick a Pattern at Random--------------------
        if (random_flag)
            for mu = 1:no_of_patterns    
                dataset_learn(mu,:) = [randi(2,1,N)-1];
            end
        else
            fid = fopen(['../../Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(l),'.mat'], 'r');                        % The path towards the dataset
            if (fid > -1 )
                fclose(fid);
                load(['../../Initialization_Files/N_',num2str(N),'_K_',num2str(K),'/neural_journal_train_set_N_',num2str(N),'_K_',num2str(K),'_index_',num2str(l),'.mat']);
            else    
                error('I can not find the learning dataset!');    
            end
            for mu = 1:no_of_patterns 
                index = 1+1 + floor(rand*2^K);
                temp = dec2bin((index),K);
                message = zeros(1,K);           % Generate the message from the index                                                         
                for j = 1:K                                    
                    message(j) = (z_max-z_min)*(temp(j) - 48)+z_min;
                end
                dataset_learn(mu,:) = [message*G];
            end 
        end
        S = max(max(abs(dataset_learn)));
        phi0 = 2*pi/(S+1);

%==========================================================================

%%
%============================LEARNING PHASE================================
        Sigma = exp(1i*phi0*dataset_learn)';
        W = Sigma*pinv(Sigma);
%==========================================================================

%%
%==============================RECALL PHASE================================
        for iki = 1:length(err_bits_range)
            pattern_error_count = 0;
            bit_error_count = 0;
            err_bits = err_bits_range(iki);            
        
            for iji = 1:no_simulated_instances    
                mu = 1 + floor(rand*no_of_patterns);
                p = dataset_learn(mu,:);     
        
                nois = zeros(1,N);
                pp = 1+floor((N-1)*rand(1,err_bits));                
                for h = 1:err_bits        
                    nois(pp(h)) =(1+floor((max_noise_amp-1)*rand))*((-1)^randi(2));            
                end        
                x = p + nois;
                x = exp(1i*x*phi0);
        
                for itr = 1:max_itr_recall        
                    ind = 1 + floor(rand*N);
                    h = W*x.';            
                    x_old = x;
                    x = csign(h*exp(1i*phi0/2),phi0);            
                    x = x.';            
                    if ( norm(abs(x-x_old))/norm(abs(x)) < 1e-5)
                        break
                    end
                end
            
                if ( sum((abs(x-exp(p*1i*phi0))>1e-3)) > 0)        
                    pattern_error_count = pattern_error_count + 1;
                    bit_error_count = bit_error_count + sum((abs(x-exp(p*1i*phi0))>1e-3));
                end    
            end 
    
            %=========================STORE THE RESULTS============================
            BER = bit_error_count/N/iji;
            PER = pattern_error_count/iji;

            if (random_flag)
                fid = fopen(['../../Recall_Results/Lee/N_',num2str(N),...
                    '_Random/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'a+');  
            else
                fid = fopen(['../../Recall_Results/Lee/N_',num2str(N),'_K_',num2str(K),...
                    '/Recall_results_Capacity_',num2str(no_of_patterns),'.txt'], 'a+');  
            end
    
            if (fid > -1)
                fprintf(fid, 'e \t %d \t per \t %f \t ber \t %f \t',err_bits,PER,BER);
                fprintf(fid,'\n');
                fclose(fid);
            else
                error('Can not store the results');
            end
            %======================================================================
        end
    end
%==================================================================================

%===============================READ THE RESULTS===================================
    run Read_Lee_Results
%==================================================================================    
end





%================================PLOT THE RESULTS==================================
if plot_flag
    figure;
    plot(sort(processed_error_bits_Lee_non_random),sort(processed_BER_Lee_non_random),'b-.','LineWidth',2,'Color','red');
    hold on        
    plot(sort(processed_error_bits_Lee_random),sort(processed_BER_Lee_random),'b-*','LineWidth',2,'Color','blue');
    legend('Subspace patterns', 'Random patterns')    
    xlhand = get(gca,'xlabel');            
    ylhand = get(gca,'ylabel');            
    set(xlhand,'string','e','fontsize',30)            
    set(ylhand,'string','Final BER','fontsize',30)
    title('Lee''s work')
end
%==================================================================================

