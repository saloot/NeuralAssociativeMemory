function montanari_bound_ldpc()

% DE parameters
del = .25;
Ndash = 100;
dv = [3;1];
dc = [5;1];
max_itr = 1000;

ldash = 3;
%[U PU V PV N] = loaddensity();
%load llrfile;
phival = [];
for p = 0.11:.01:.14,
    
    
    [Perr, PV, PU] = DE_irregLDPC_BSC(p,dv,dc,del,Ndash,max_itr);
    
    N = length(PU);
    
    U = [del*(-Ndash:Ndash) Inf];
    V = U;
    
    temp = phicalc(p,ldash,U,PU,V,PV,N);
    phival = [phival temp]
%     if (temp > 0),
%         break;
%     end
end
phival
% fprintf('phival = %f',phival);
% fprintf('BSC threshold is %f',p);
return;

function sum1 = phicalc(p,ldash,U,PU,V,PV,N)
sum1 = 0;
fprintf('Total iterations = %d',N^ldash);
for i = 0:N^ldash-1,
    str = dec2bigbase(i,N,ldash);
    prod0 = 1;
    prod1 = 1;
    prob = 1;
    for j = 1:ldash,
        prod0 = prod0*func0(U(str(j)+1));
        prod1 = prod1*func1(U(str(j)+1));
        prob = prob*PU(str(j)+1);
    end
    sum1 = sum1 + prob*(.5*log2(prod0 + p/(1-p)*prod1)+.5*log2(prod0 + (1-p)/p*prod1) );
    if (i < N^2),
        str = dec2bigbase(i,N,2);
        ind1 = str(1)+1;
        ind2 = str(2)+1;
        prob = PU(ind1)*PV(ind2);
        sum1 = sum1 - ldash*prob*log2( func0(U(ind1))*func0(V(ind2)) + func1(U(ind1))*func1(V(ind2)) );
    end
    
    if (mod(i,100000) == 0),
        i
    end
end
return;
    
        
function f = func0(u)
f = (1+exp(-2*u))^(-1);
return;

function f = func1(u)
f = ((1+exp(-2*u))^(-1))*exp(-2*u);
return;


