%this is the 2nd one 
%in marboto be golde 5 hastaesh
clc
%clear all;


%=============================INITIALIZATION===============================
N = 2^7-1;                     % N is the number of neurons in network.
k = 8;                       % k is the number of message bits.
M = N;

theta = 0;                   % theta is the neural threshold.
max_itr = 1;                 % This is the maximum number of iteration for the simulation.
deg_ave = 3;                 % This is the average degree of the LDGM code.error_count = 0;             % This counts the number of errors over iterations.

beta0 = 1000;                   % beta determines the degree of stochasticity of the Hopfield network
error_count = zeros(1,21);             % This tracks the number of errors.
bit_error_count = zeros(1,21);

error_count_stochastic = zeros(1,21);             % This tracks the number of errors.
bit_error_count_stochastic = zeros(1,21);

W = zeros(N,N);              % W is the weight matrix.
Y = zeros(M,N);              % Y is the matrix of patterns;
x = zeros(1,N);              % x is the vector of states;
%==========================================================================



%====================PICK PATTERNS AND THE INITIAL STATE===================

m = 7;
N = 2^m -1;
l = 3;
M = N;
alph = gf([2],m);
% G = gold(m,l);
load gold;
% load gold_q_9_deterministic
Y(1:N,:) = G(1:N,:);
% Y(N+1:2*N,:) = G(N+2:2*N+1,:);

%-------------------Calculate Correlation of Patterns--------------
%[corr,pd,ax] = correlation_patterns(Y);
%------------------------------------------------------------------


%=========================DETERMINE THE OMEGA'S========================
Omega = zeros(1,M);

% Omega(M) = 1;       %%%%%%%%%%%%%%
for mu = 0:N-1
    aa = gftrace(alph^(-mu),m);
    Omega(mu+1) = (-1)^(double(aa.x));
    bb = gftrace(alph^(-mu-1),m);
    Omega(N+mu+1) = (-1)^(double(bb.x));
end
%======================================================================



%=======================DETERMINE THE WEIGHTS==========================
for i = 1:N
    for j = 1:N
        summ = 0;
        for m = 1:M
            summ = summ + Omega(m)*Y(m,i)*Y(m,j);
        end
        W(i,j) = summ;
    end
end
W = W/M;

%======================================================================
error_itr = 0;
for err_bits = 0:5:20                % This is the number of erroneous bits
    error_itr = error_itr + 1;
%=========INTIALIZE THE NETWORK WITH A PATTERN AND ITERATE=============
error_ind = [];
for kk = 1:5000
    
    %------------------------------Add Noise---------------------------
    nois = ones(1,N);
    %         for nn = 1:N
    %             pp = rand;
    %             if (pp<nois_prob)
    %                 nois(nn) = -1;
    %             end
    %         end
    pp = 1+round((N-1)*rand(1,err_bits));
    for h = 1:err_bits
        nois(pp(h)) = -1;       
    end
    %------------------------------------------------------------------
    
    
    p = round((M-1)*rand)+1;
    x = Y(p,:).*nois; % Initialize the network with the first pattern.
    x0 = Y(p,:);
    x_last = zeros(1,N);
    
    %--------Update Neurons' States Synchronously and Iteratively------
    itr = 0;
    exit_flag = 0;
    111
    pause
    while (exit_flag == 0)
        itr = itr + 1;
        
        %------------------Calculate the Weighted Input Sum-----------------
        x_temp = W*x';
        
%         beta = max(5,beta0/(10*itr));
        beta = beta0;
        for iii = 1:N
            pp = rand;
            if (pp > 1/(exp(-2*beta*x_temp(iii))))
                x(iii) = -1;
            else
                x(iii) = 1;
            end
        end
        
        
        
        % Remain in the loop until convergence or until we exceed 1000
        % iterations:
        if ( ((mod(itr,2) == 0) && (norm(x-x_last) == 0)&(itr>0))||(itr == 1000))
            exit_flag = 1;
%             itr
%             pause
%             clc
        end
        
        
        if ((mod(itr,2) == 0)&(exit_flag == 0))           
            x_last = x;            
%         else
%             11111
%         (norm(x-x0*Omega(p))^2)/4
        end
%        pause    
%         
    end
    %------------------------------------------------------------------    
    if (norm(x - x0)~=0)
        error_count(error_itr) = error_count(error_itr) + 1;
        error_ind = [error_ind p];
        bit_error_count(error_itr) =  bit_error_count(error_itr)+ (norm(x-x0)^2)/4;       
    end
    
    if (mod(kk,100) == 0)
        display(['Iteration:',num2str(kk), ' Error count so far: ',num2str(error_count(error_itr))]);
    end
    
end
%==========================================================================


display(' ');
display(['Pattern Error Rate ',num2str(error_count(error_itr)/kk)]);
display(['Bit Error Rate ',num2str(bit_error_count(error_itr)/kk/N)]);
end


