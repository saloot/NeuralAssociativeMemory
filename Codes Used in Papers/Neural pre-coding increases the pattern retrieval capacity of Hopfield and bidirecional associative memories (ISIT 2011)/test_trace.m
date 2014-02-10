clc

%=============================INITIALIZATION===============================
N = 2^7-1;                     % N is the number of neurons in network.
k = 8;                       % k is the number of message bits.
M = N;



Y = zeros(M,N);              % Y is the matrix of patterns;

%==========================================================================



%====================PICK PATTERNS AND THE INITIAL STATE===================

m = 7;
N = 2^m -1;
l = 3;
M = N;
alph = gf([2],m);


Tr_value = 2*ones(N,M);
for i = 0:N-1
    for mu = 0:N-1
        if (gftrace((alph^i)*alph^mu+(alph^(-mu)),m) == gf([0],m))
            Tr_value(i+1,mu+1) = 0;
        else
            Tr_value(i+1,mu+1) = 1;
        end
    end
end
% G = gold(m,l);
% load gold;
% % load gold_q_9_deterministic
% Y(1:N,:) = G(1:N,:);
% % Y(N+1:2*N,:) = G(N+2:2*N+1,:);
% 
% %-------------------Calculate Correlation of Patterns--------------
% %[corr,pd,ax] = correlation_patterns(Y);
% %------------------------------------------------------------------
% 
% 
% %=========================DETERMINE THE OMEGA'S========================
% Omega = zeros(1,M);
% 
% % Omega(M) = 1;       %%%%%%%%%%%%%%
% for mu = 0:N-1
%     aa = gftrace(alph^(-mu),m);
%     Omega(mu+1) = (-1)^(double(aa.x));
%     bb = gftrace(alph^(-mu-1),m);
%     Omega(N+mu+1) = (-1)^(double(bb.x));
% end
% %======================================================================
% 
% 
% 
% %=======================DETERMINE THE WEIGHTS==========================
% for i = 1:N
%     for j = 1:N
%         summ = 0;
%         for m = 1:M
%             summ = summ + Omega(m)*Y(m,i)*Y(m,j);
%         end
%         W(i,j) = summ;
%     end
% end
% W = W/M;