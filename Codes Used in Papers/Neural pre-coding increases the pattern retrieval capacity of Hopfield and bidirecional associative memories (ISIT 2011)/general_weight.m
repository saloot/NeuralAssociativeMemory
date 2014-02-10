clc
clear all;
error_count = 0;

%=============================INITIALIZATION===============================
N = 100;                     % N is the number of neurons in network.
M = 1*N;
theta = 0;                  % theta is the neural threshold. 
max_itr = 40000;              % This is the maximum number of iteration for the simulation.
mu = -.1;
sigma = .8;
epsilon0 = 0;


W = zeros(N,N);             % W is the weight matrix.
Y = zeros(M,N);             % Y is the matrix of patterns;
x = zeros(1,N);             % x is the vector of states;
%==========================================================================


%====================PICK PATTERNS AND THE INITIAL STATE===================
Y = randint(M,N);
Y = -2*Y + ones(M,N);

min_distance = N;
for m = 1:M
    for j = m+1:M
        dist = sum(abs(Y(m,:)-Y(j,:)));
        if (dist < min_distance)
            min_distance = dist;
            ind = [m,j];
        end
    end
end

min_distance   
%--------------------------------------------------------------------------


%==========================================================================



%===========================DETERMINE THE WEIGHTS==========================
for i = 1:N
    for j = 1:N
        p = rand;
        if (p < epsilon0)
            W(i,j) = 0;
        else
            W(i,j) = mu + sigma*randn;
        end        

    end
end
        


for i = 1:N                 
    W(i,i) = 0;             % Remove the self loops.
end
%==========================================================================


for kk = 1:100
% Y = 2*randint(M,N)-ones(M,N); % Pick M random N-bit patterns. 
p = round((M-1)*rand)+1;
x = Y(p,:);                 % Initialize the network with the first pattern.
x0 = x;
%======================UPDATE NEURONS' STATE ITERATIVELY===================
for r = 1:max_itr           
    
        
       %------------------Calculate the Weighted Input Sum-----------------
       x_temp = W*x';
       %-------------------------------------------------------------------
       
       %----------Compare to the Threshold and Maje the Decision-----------
       x = sign(x_temp'-theta*ones(1,length(x)));   
       for jj = 1:N
           if (x(jj) == 0)
               pp = rand;
               if (pp <= -1)
                   x(jj) = 1;
               else
                   x(jj) = -1;
               end
                 
           end
       end
       %-------------------------------------------------------------------
    
end
%==========================================================================

if (norm(x - x0)~=0)
    clc
%     display('Error');
     (norm(x-x0))^2
error_count = error_count + 1
kk

end
kk

end
display(error_count);
