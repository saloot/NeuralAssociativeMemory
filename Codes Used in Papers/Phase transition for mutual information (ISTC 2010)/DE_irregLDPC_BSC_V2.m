%==========================================================================
%=========================DENSITY EVOLUTION FOR BSC========================
%==========================================================================
%--------------------Function Name: DE_regLDPC_BSC()-----------------------
%----------------------------------Inputs----------------------------------
% p: BSC transition probability
% dv: Variable node degree distribution. The first row contains the degrees...
% ...and the second row contains coefficients.
% dc: Check node degree distribution. The first row contains the degrees...
% ...and the second row contains coefficients.
% del: Quantization step
% N: Number of quantization levels
% max_itr: Maximum number of density evolution iterations
%--------------------------------------------------------------------------

%----------------------------------Output----------------------------------
% Perr: Probability of error in each iteration of density evolution.
%--------------------------------------------------------------------------

%--------------------------Function Descritption---------------------------
% The function DE_regLDPC_BSC(p,dv,dc) 
%--------------------------------------------------------------------------


function Perr = DE_irregLDPC_BSC(p,dv,dc,del,N,max_itr)

%------------------------Check the Validity of Input-----------------------
if ((p <0)|(p>1)|(del <=0)|(N<=0)|(numel(dv)==0)|(numel(dc)==0))
    error('Error in the input of the density evolution function!');
end
if ((abs(sum(dc(2,:))-1)>1e-6) | (abs(sum(dv(2,:))-1)>1e-6))
    error('Error in input degree distribution');
end
%--------------------------------------------------------------------------

%------------------------------Initialiation-------------------------------
LLRvals = del*(-N:N);               % LLR density parameters

P0 = zeros(1,2*N+1);
[val1 pos1] = min(abs(LLRvals - ( -log((1-p)/p) )));
[val2 pos2] = min(abs(LLRvals - (  log((1-p)/p) ))); 
P0(pos1) = p;
P0(pos2) = 1-p;
P0 = [P0 0];    % Add the mass point at +infinity
Pvc = P0;
Perr = zeros(1,max_itr);
itr = 0;
end_flag = 0;
Pvc = P0;
%--------------------------------------------------------------------------


%-------------------------Computhe the Index Matrix------------------------
T = tgen(N,del);
%--------------------------------------------------------------------------


%---------------------------Do Density Evolution---------------------------
while (end_flag == 0),
    itr = itr+1
    
    %-------------Determine Check to Variable Node Messages----------------
    Pcv_final = zeros(1,2*N+2);
    for i = 1:length(dc(1,:))
    
        Pcv = conv_chk(Pvc,Pvc,del,T);        % Need dc >= 3
        if (dc(1,i)>=3)
            for cnt1 = 1:dc(1,i)-3,
                Pcv = conv_chk(Pcv,Pvc,del,T);
            end
        end
        Pcv_final = Pcv_final + dc(2,i)*Pcv;
    end
    Pcv = Pcv_final;
    %----------------------------------------------------------------------
    
    %-------------Determine Variable to Check Node Messages----------------
    Pvc_final = zeros(1,2*N+2);
    for i = 1:length(dv(1,:))
        Pvc = conv_var(P0,Pcv);        % Need dv >= 2
        if (dv(1,i)>=2)
            for cnt1 = 1:dv(1,i)-2,
                Pvc = conv_var(Pvc,Pcv);
            end
        end
        Pvc_final = Pvc_final + dv(2,i)*Pvc;
    end
    Pvc = Pvc_final;
    %----------------------------------------------------------------------
    
    
    %------------------Calculate the Probability of Error------------------
    Perr(itr) = sum(Pvc(1:N)) + 0.5*Pvc(N+1);
    Perr(itr)
    %-------------Determine Check to Variable Node Messages----------------
    
    if ( (itr>max_itr))
        end_flag = 1;
    elseif  (itr >=2)
        if (abs(Perr(itr)-Perr(itr-1)) <1e-8 )
            end_flag = 1;
        end
    end
end
Perr = Perr(1:itr);
%--------------------------------------------------------------------------

return;
%==========================================================================
%==========================================================================
%==========================================================================



%==========================================================================
%=========================VARIABLE NODE CONVOLUTION========================
%==========================================================================
%----------------------Function Name: conv_var(a,b)------------------------
%----------------------------------Inputs----------------------------------
% a: One of the input densities
% b: The other one of the input densities
%--------------------------------------------------------------------------

%----------------------------------Output----------------------------------
% c: The resulted density of convoluting two input densities. 
%--------------------------------------------------------------------------

%--------------------------Function Descritption---------------------------
% The function conv_var(a,b) gets two probability density functions and 
% returns their convolution using Fast Fourier Transforms (FFTs) according
% to the method described in section B.2 of the book "Modern Coding Theory"
% by Rudiger Urbanke.
%--------------------------------------------------------------------------

function c = conv_var(a,b)

%------------------------Check the Validity of Input-----------------------
if ((abs(sum(a)-1)>.0000000001)|(abs(sum(b)-1)>.0000000001))
    error('Error in variable update function input!');
end
%--------------------------------------------------------------------------

%----------------------------Initilization---------------------------------
M = length(a) - 1;      % -> 2N+1 in the book
a1 = a(1:M);            % Only consider the mass points in the interval [-N,N] and not infinity.
b1 = b(1:M);            % Only consider the mass points in the interval [-N,N] and not infinity.
%--------------------------------------------------------------------------


%-----------------------Compute the Convolution Using FFT------------------
c1 = ifft(fft(a1,2*M-1).*fft(b1,2*M-1));    % of length 2M-1 = 4N+1
negsum = sum(c1(1:(M-1)/2));                % sum of probs of elements correspoding to -2N to -N-1
possum = sum(c1(2*M-((M-1)/2):2*M-1));      % sum of probs of elements correspoding to N+1 to 2N
c = c1((M-1)/2+1:2*M-((M-1)/2)-1);          % elements corresponding to -N to N (length M)
c(1) = c(1) + negsum;
c = [c possum+( a(M+1)+b(M+1)-a(M+1)*b(M+1) )]; % length M+1
%--------------------------------------------------------------------------

%-----------------------Check the Validity of the Output-------------------
if (sum(c)~=1)    
    
    %----------------------Check for Numerical Issues----------------------
    if (abs(sum(c)-1)>1e-6)
        error('Error in the convoluted density of the varibale update function');
    else                                    %Take care of small numerical issues
        c = max(c,zeros(1,M+1));
        c(M+1) = c(M+1) + 1-sum(c);
    end
    %----------------------------------------------------------------------
    

end
%--------------------------------------------------------------------------

return;
%==========================================================================
%==========================================================================
%==========================================================================



%==========================================================================
%===========================CHECK NODE CONVOLUTION=========================
%==========================================================================
%--------------------Function Name: conv_chk(a,b,del,T)--------------------
%----------------------------------Inputs----------------------------------
% a: One of the input densities
% b: The other one of the input densities
% del: Quantization step of probability density functions
% T: A matrix that stores the index k such that T(i,j)=k if Q(i,j) = k. See
% section B.3 of the book "Modern Coding Theory" by Rudiger Urbanke.
%--------------------------------------------------------------------------

%----------------------------------Output----------------------------------
% c: The resulted density of convoluting two input densities. 
%--------------------------------------------------------------------------

%--------------------------Function Descritption---------------------------
% The function conv_chk(a,b,del,T) gets two probability density functions and 
% returns their convolution using the table method according to the approach 
% described in section B.2 of the book "Modern Coding Theory" by Rudiger Urbanke.
% Here the matrix T stores the index k such that T(i,j)=k if Q(i,j) = k.
%--------------------------------------------------------------------------

function c = conv_chk(a,b,del,T)

%------------------------Check the Validity of Input-----------------------
if ((abs(sum(a)-1)>.000001)|(abs(sum(b)-1)>.000001))
    error('Error in Check update function input!');    
end
if (del <=0)
    error('Error in the step size in the Check update function input!');    
end
%--------------------------------------------------------------------------

%------------------------------Initialization------------------------------
N = (length(a)-2)/2;
temp = a(N+1:2*N+1) + [0 a(N:-1:1)];   % length N+1
aplus = [temp a(2*N+2)];                 % including infinity at mass points
temp = [0 a(N+2:2*N+1)-a(N:-1:1)];    % length N+1
aminus = [temp a(2*N+2)];
temp = b(N+1:2*N+1) + [0 b(N:-1:1)];   % length N+1
bplus = [temp b(2*N+2)];
temp = [0 b(N+2:2*N+1)-b(N:-1:1)];    % length N+1
bminus = [temp b(2*N+2)];
cplus = [zeros(1,N+1) a(2*N+2)*b(2*N+2)];
cminus = cplus;
%--------------------------------------------------------------------------


%-------------------------Compute the Convolution--------------------------
for i =0:N+1    
    for j = 0:N+1
        if ((j<=N)|(i<=N))              % Ignore the mass at infinity at this point            
            cplus(T(i+1,j+1)+1) = cplus(T(i+1,j+1)+1) + aplus(i+1)*bplus(j+1);
            cminus(T(i+1,j+1)+1) = cminus(T(i+1,j+1)+1) + aminus(i+1)*bminus(j+1);                               
        end
    end
end

c = [0.5*(cplus(N+1:-1:2)-cminus(N+1:-1:2)) cplus(1) 0.5*(cplus(2:N+2)+cminus(2:N+2))]; % Have extended the same recovery mechanism to c_infinity also, CHECK!
%--------------------------------------------------------------------------




%-----------------------Check the Validity of the Output-------------------
if (sum(c)~=1)    
    
    %----------------------Check for Numerical Issues----------------------
    if (abs(sum(c)-1)>1e-6)
        error('Error in the convoluted density of the check update function');
    else                                    %Take care of small numerical issues
        c = max(c,zeros(1,length(c)));
        c(2*N+2) = c(2*N+2) + 1-sum(c);
    end
    %----------------------------------------------------------------------
    

end
%--------------------------------------------------------------------------


return;
%==========================================================================
%==========================================================================
%==========================================================================



%==========================================================================
%===========================INDEX MATRIX GENERATION========================
%==========================================================================
%----------------------Function Name: tgen(N,del)--------------------------
%----------------------------------Inputs----------------------------------
% N: number of quantization levels
% del: Quantization step of probability density functions
%--------------------------------------------------------------------------

%----------------------------------Output----------------------------------
% T: A matrix that stores the index k such that T(i,j)=k if Q(i,j) = k. See
% section B.3 of the book "Modern Coding Theory" by Rudiger Urbanke.
%--------------------------------------------------------------------------

%--------------------------Function Descritption---------------------------
% The funtion tgen(N,del) gets the number of quantization levels and the
% step between two consecutive levels and returns the matrix T such that
% T(i,j) = k if Q(i,j) = k. For more information about the deifinitions see
% section B.3 of the book "Modern Coding Theory" by Rudiger Urbanke.
%--------------------------------------------------------------------------

function T = tgen(N,del)
%------------------------Check the Validity of Input-----------------------
if ((del <=0)|(N<=0))
    error('Error in the index matrix generation function input!');    
end
%--------------------------------------------------------------------------

%-------------------------------Initialization-----------------------------
T = zeros(N+2,N+2);
%--------------------------------------------------------------------------


%--------------------------Determine the Index Matrix----------------------
for i = 0:N+1
    for j = 0:N+1
        if ((i == N+1)&(j ~=N+1))           % Takes care of Q(infintiy,j)
            Q = j;
        elseif ((j == N+1)&(i ~=N+1))       % Takes care of Q(i,infintiy)
            Q = i;
        elseif ((j<=N) &(i<=N))             % Takes care of other cases 
            Q = floor( (2*atanh(tanh(0.5*i*del)*tanh(0.5*j*del))/del) + 0.5 );
        else                                % Takes care of Q(infintiy,infinity)
            Q = -1;
        end
        
        flag = 0;
        for k = 0:N
            
            %----------------Check for Numerical Issues in Q---------------
            if (isreal(Q) == 0),
                error('Numerical trouble in the index matrix generation!');
            end
            %--------------------------------------------------------------
            
            if (Q == k)                                
                T(i+1,j+1) = k;  
                flag = 1;
                break;
            end
            
        end
        
        %--------------Check the Validity of the Output--------------------
        if ((flag == 0)&(Q~=-1))
            error('Error in index matrix generation output!');                        
        end
        %------------------------------------------------------------------
    end
end
if (sum(sum(isinf(Q)))>0),
    error('Infinite values for Q owing to finite precision effects, decrease del or N!');
end
return;
%==========================================================================
%==========================================================================
%==========================================================================