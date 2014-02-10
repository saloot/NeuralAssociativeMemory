clc
%clear all;

%=============================INITIALIZATION===============================
N = 12;                     % N is the number of neurons in network.
k = 4;                      % k is the number of message bits.
M = 2^k;                    % M is the number of patterns.
deg_ave = 2;                % This is the average degree of the LDGM code.
beta = 4.5;                  % beta is the coefficient of stochastic networks in the sigmoid function
theta = 0;                  % theta is the neural threshold. 
max_itr = 1;               % This is the maximum number of iteration for the evolution of network for a given input pattern
max_pattern_itr = 100;      % Maximum number of patterns initiated 
max_no_of_codes = 50;       % Maximum number of codes over the ensemble determined by deg_ave
error_count = 0;            % Counts the number of errors
nois_prob = 0.01;           % This is the BSC cross-over probability
sjs = 3;

Y = zeros(M,N);             % Y is the matrix of patterns;
x = zeros(1,N);             % x is the vector of states;
%==========================================================================


%====================PICK PATTERNS AND THE INITIAL STATE===================


for code_loop = 1:max_no_of_codes
    temp_err_count = 0;
    average_distance = 0;
    %--------------------------Build the Generator Matrix----------------------
    G = zeros(k,N);
    for j = 1:N
        one_positions = 1+round((k-1)*rand(deg_ave,1));         % Randomly pick deg_ave positions in each column and make them equal to 1
        for i = 1:deg_ave        
           G(one_positions(i),j) = 1;
        end
    end
    %--------------------------------------------------------------------------

    %--------------------------Construct the Codewords-------------------------
    for i = 0:M-1
        temp = dec2bin(i,k);
        message = zeros(1,k);
        for j = 1:k
            message(j) = temp(j) - 48;                          % Mapping from ASCII to digit
        end
        Y(i+1,:) = mod(message*G,2);
    end
    %--------------------------------------------------------------------------

    %--------------------------Calculate Minimum Distance----------------------
    min_dist = 10*N;
    for m = 1:M
        for j = m+1:M
            average_distance = average_distance + sum(abs(Y(m,:)-Y(j,:)));
            if ( sum(abs(Y(m,:)-Y(j,:))) < min_dist)
                min_dist = sum(abs(Y(m,:)-Y(j,:)));
                ind = [m,j];
            end
        end
    end
    Y = -2*Y + ones(M,N);                                   % Change from 0/1 to +1/-1
    average_distance = average_distance/M;
    %--------------------------------------------------------------------------


    %=========================DETERMINE THE WEIGHTS============================
    A = zeros(M*N,N^2);
    for m = 1:M                       
        for i = 1:N
            for j = 1:N
                A((m-1)*N+i,N*(i-1)+j) = Y(m,i)*Y(m,j);
                if (j == sjs)
                    A((m-1)*N+i,N*(i-1)+j) = -Y(m,i)*Y(m,j);
                end
            end
        end
    end
    
    

    b = zeros(M*N,1);
%     for m = 1:M
%         b((m-1)*N+1:m*N) = Y(m,:)';
%     end
    b = -10*ones(M*N,1)/(2*beta);

    options = optimset('MaxIter',1000,'Display','off');
    
    f=zeros(1,N^2);
    W = linprog(2*f,-A,b,[], [], [], [], [],options);         % Weights are transformed from a matrix form to a row vector form...
     W = W/max(W);                                           % ...by writing each row of the matrix after the previous one in a vector. 
    W_mat = zeros(N,N);
    for i = 1:N
        W_mat(i,:) = W((i-1)*N+1:i*N);
    end
%     max(abs(W_mat - eye(N)))
%     pause
    %======================================================================

   
    %====================SIMULATE THE NETWORK BEHAVIOUR====================
    for pattern_loop = 1:max_pattern_itr    % Pick max_papattern_itr patterns randomly as the input pattern and simulate.
        
        %-------------------Initialize with the First Pattern--------------
        p = round((M-1)*rand)+1;
        x = Y(p,:);                 % Initialize the network with the first pattern.
        x0 = x;
        %------------------------------------------------------------------
        
        %----------------------Add Some BSC Noise--------------------------
        nois = ones(1,N);
        for j = 1:N
            pp = rand;
            if (pp < nois_prob)
                nois(j) = -1;
            end
        end
        %x = x .*nois;
        
        x(sjs) = -x(sjs);
        %------------------------------------------------------------------

        %------------------Update Neural States Serially-------------------
        x_last = zeros(1,N);
        itr = 0;
        while (sum(abs(x_last-x)) ~=0)
            itr = itr+1;
            x_last = x;
            for i = 1:N             % Update the state of each neuron
               summ = 0; 
               %------------------Calculate the Weighted Input Sum---------
               for j = 1:N
                  summ = summ + W((i-1)*N+j)*x(j); 
               end
               %-----------------------------------------------------------
       
               %-------Compare to the Threshold and Maje the Decision------
               x(i) = sign(summ);
               if (summ == 0)
                   x(i) = 1;
               end
               q = rand;
               if (q <= 1/(1+exp(-2*beta*summ)))        % Stochastic network
                   x(i) = 1;
               else
                   x(i) = -1;
               end
               %-----------------------------------------------------------
            end
        end
        
        %------------------------------------------------------------------

        %------------------------Look for Errors---------------------------
        
        
        if (norm(x - x0)~=0)
            temp_err_count = temp_err_count + 1;
            
        end       
        %------------------------------------------------------------------        

    end
    %======================================================================

    display(['Error rate for code ', num2str(code_loop),' = ',num2str(temp_err_count/max_pattern_itr)]);
    display(['Minimum distance = ', num2str(min_dist),' and average distance = ',num2str(average_distance)]);
    display('  ');
    if ((temp_err_count/max_pattern_itr >= .4)||(temp_err_count/max_pattern_itr <0))
        W_mat(:,3)
    pause
    end
    error_count = error_count + temp_err_count;    
end
display('------------------------------------------------------------------');
display(['Error rate = ',num2str(error_count/(max_pattern_itr*max_no_of_codes))]);
display('   ');
display(['Probability of error = ',num2str(nois_prob),' and beta = ',num2str(beta)]);
display('   ');
display(['Number of errors = ',num2str(error_count)]);
display('   ');
display(['Out of ',num2str(max_pattern_itr),' iterations over ',num2str(max_no_of_codes),' codes.']);
display('------------------------------------------------------------------');