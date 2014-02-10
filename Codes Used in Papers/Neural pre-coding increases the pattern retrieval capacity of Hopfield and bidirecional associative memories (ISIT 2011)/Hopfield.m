clc
%clear all;


%=============================INITIALIZATION===============================

N = 2^7-1;                     % N is the number of neurons in network.
k = 8;                       % k is the number of message bits.
M = (N+1)/2;

theta = 0;                   % theta is the neural threshold. 
max_itr = 2;                 % This is the maximum number of iteration for the simulation.
deg_ave = 3;                 % This is the average degree of the LDGM code.error_count = 0;             % This counts the number of errors over iterations.
err_bits = 5;                % This is the number of erroneous bits
beta0 = 10;                   % beta determines the degree of stochasticness of the Hopfield network
error_count = 0;             % This tracks the number of errors.
ensemble_flag = 2;           % 0 means random patterns are used. 1 is for codewords. 2 is for a given set of patterns
convergence_flag = 1;        % 0 means you repeat for max_itr times. 1 means you wait untill convergence or do 1000 iteration

W = zeros(N,N);              % W is the weight matrix.
Y = zeros(M,N);              % Y is the matrix of patterns;
x = zeros(1,N);              % x is the vector of states;
%==========================================================================


if (ensemble_flag == 1)
    no_codes_generated = 100;% This is the number different random codes generated by the progra
else
    no_codes_generated = 1;
end
%====================PICK PATTERNS AND THE INITIAL STATE===================



for out_loop = 1:no_codes_generated
    
    if (ensemble_flag == 0)
        Y = randint(M,N);
        Y = -2*Y+ones(M,N);    
    elseif (ensemble_flag == 1)
        M = 2^k;                    
        %----------------------Build the Generator Matrix------------------
        G = zeros(k,N);
        for j = 1:N
            for i = 1:k
                p = rand;
                if (p<deg_ave/k)
                    G(i,j) = 1;
                end
            end
        end
        %------------------------------------------------------------------

        %------------------Construct the Codewords-------------------------       
        for i = 0:M-1  
            temp = dec2bin(i,k);
            message = zeros(1,k);
            for j = 1:k
                message(j) = temp(j) - 48;
            end            
            Y(i+1,:) = mod(message*G,2);    
        end
        %------------------------------------------------------------------
        Y = -2*Y + ones(M,N);
        YY = Y';

%         %-----------------------Flip One Bit Randomly----------------------
%         for m = 1:M
%             prand = 1+round((N-1)*rand);
%             Y(m,prand) = -1*Y(m,prand);
%         end
%         %------------------------------------------------------------------

%         %--------------------------Pick the First M Elements-------------
%         M_orig = M;
%         M = .2*N;
%         Y = Y(1:M,:);
%         %----------------------------------------------------------------

    else
        m = 7;
        N = 2^m -1;
        alph = gf([2],m);
        %run gold
       
        Y = G(1:M,:);

        %-------------------Calculate Correlation of Patterns--------------
        %[corr,pd,ax] = correlation_patterns(Y);
        %------------------------------------------------------------------
    end
    

    %------------------------Calculate Minimum Distance--------------------
    min_distance = 10*N;
    for m = 1:M
        for j = m+1:M        
            dist = sum(abs(Y(m,:)-Y(j,:)));
            if ((dist < min_distance))
                min_distance = dist;
                ind = [m,j];
            end                
        end            
    end
    %----------------------------------------------------------------------
    
    
   


    %======================================================================


    %=========================DETERMINE THE OMEGA'S========================
    Omega = zeros(1,M);
    alph = gf([2],7);

    %Omega(M) = 1;
    for mu = 0:M-1
        aa = gftrace(alph^(-mu),7);      
        
        if ((-1)^(double(aa.x))==1)
            Omega(mu+1) = 1;
        else
            Omega(mu+1) = -1;
        end
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
    W = W/N;        

%   for i = 1:N                 
%       W(i,i) = 0;             % Remove the self loops.
%   end
%     A = zeros(N,N^2);
%     B = [];
%     for m = 1:M
%         for i = 1:N
%             for j = 1:N
%                 A(i,(i-1)*N+j) = Y(m,i)*Y(m,j);
%             end
%         end
%         B = [B;A];
%     end
%     %b = -.0001*ones(M*N,1);
%     
%     AA = zeros(N,N^2);
%     for i = 1:N
%         for j = 1:N
%             if (j == 1)
%                 A(i,(i-1)*N+j) = -1;
%             else
%                 A(i,(i-1)*N+j) = 1;
%             end
%         end
%     end
%     B = [B;AA];
%     
%     b = -.0001*ones(M*N+N,1);    
%     
%     ww = linprog(zeros(N^2,1),-B,b);

    %======================================================================

    %=========INTIALIZE THE NETWORK WITH A PATTERN AND ITERATE=============
    error_ind = [];
    for kk = 1:12000   
        
        %------------------------------Add Noise---------------------------
        nois = ones(1,N);
        pp = 1+round((N-1)*rand(1,err_bits));        
        for h = 1:err_bits
            nois(pp(h)) = -1;       
        end
        %------------------------------------------------------------------
    
        
        p = round((M-1)*rand)+1;   
        x = Y(p,:).*nois;                 % Initialize the network with the first pattern.
        x0 = Y(p,:);
        x_last = zeros(1,N);
        %--------Update Neurons' States Synchronously and Iteratively------
        itr = 0;
        exit_flag = 0;
        while (exit_flag == 0)                                        
            itr = itr + 1;
       
            %------------------Calculate the Weighted Input Sum-----------------
             x_temp = W*x';    


             beta = beta0*itr;
             x_last_last = x_last;
             x_last = x;
             
            for iii = 1:N
                pp = rand;
                if (pp > 1/(1+exp(-2*beta*(x_temp(iii))-theta)))
                    x(iii) = -1;
                else
                    x(iii) = 1;
                end
            end
            
            
%           Omega(p)
%              norm(x-x0)^2
%              norm(x+x0)^2
%               pause
            
            if (convergence_flag == 0)      %If 0, do max_itr iterations and quite
                if (itr >=max_itr)
                    exit_flag = 1;                    
                end
            else                            % Else, remain in the loop until convergence or exceeding 1000 iterations.
                if (((norm(x-x_last_last) == 0)&(mod(itr,2)==0)&(itr>40))|(itr > 1000))
                    exit_flag = 1;        
                    
                end
            end
          
        end
        %------------------------------------------------------------------
        
        
            
       
       if ((norm(x - x0)~=0))
           error_count = error_count + 1;          
           error_ind = [error_ind p];            

       end

       if (mod(kk,100) == 0)
           display(['Iteration:',num2str(kk), ' Error count so far: ',num2str(error_count)]); 
       end
       
    end
   %==========================================================================

   

end

display(['Number of errors out of ',num2str(kk),' iterations over ',num2str(out_loop),' random codes:']);
display(error_count);
