function montanari_bound_ldpc_BEC()

Nsimul = 4000000;    % # of times monte-carlo sampling is performed in evaluating Montanari's bound

% DE parameters
del = .25;
Ndash = 100;
dv = [3;1];
dc = [6;1];
max_itr = 1000;

ldash = 3;
kdash = 6;

phival = [];
for e = 0.45:.05:.60,
    
    cd('/Users/raj/Documents/Matlab work/Density Evolution');
    [Perr, PV, PU] = DE_irregLDPC_BEC(e,dv,dc,del,Ndash,max_itr);
    cd('/Users/raj/Documents/Matlab work/Modern codes')
    N = length(PU);
    U = 0.5*[del*(-Ndash:Ndash) Inf];   % The 0.5 up-front is to compensate for Montanari's strange definition of LLRs
    V = U;
    
    [PU IX] = sort(PU,'descend');
    U = U(IX);
    [PV IX] = sort(PV,'descend');
    V = V(IX);
    
    temp = phicalc(e,ldash,kdash,U,PU,V,PV,Nsimul);
    phival = [phival temp]
%     if (temp > 0),
%         break;
%     end
end
phival
% fprintf('phival = %f',phival);
% fprintf('BSC threshold is %f',p);
return;

function sum1 = phicalc(e,ldash,kdash,U,PU,V,PV,Nsimul)

neglect_val = 1e-6;

sum1 = 0;
indices1 = find(PU > neglect_val);
indices2 = find(PV > neglect_val);
indices = union(indices1,indices2);
N = length(indices);

fprintf('\nTotal number of significant points in the pdf = %d',N);
fprintf('\nIterations needed for an exact result = %d',N^max(ldash,kdash));
fprintf('\nNumber of Monte-Carlo iterations = %d',Nsimul);
%pause;

for cnt = 1:min(Nsimul,N^max(ldash,kdash)),
    
    if (cnt <= N^ldash),
        % str = randi([1 N],1,ldash);
        str = dec2bigbase(cnt-1,N,ldash) + 1;
        prod0 = 1;
        prod1 = 1;
        prob = 1;
        for j = 1:ldash,
            prod0 = prod0*func0(U(indices(str(j))));
            prod1 = prod1*func1(U(indices(str(j))));
            prob = prob*PU(indices(str(j)));
        end
        sum1 = sum1 + (1-e)*prob*log2(prod0) + e*prob*log2(prod0 + prod1);
    end
    
    if (cnt <= N^2),        
        %str = randi([1 N],1,2);
        str = dec2bigbase(cnt-1,N,2) + 1;
        ind1 = indices(str(1));
        ind2 = indices(str(2));
        prob = PU(ind1)*PV(ind2);
        sum1 = sum1 - ldash*prob*log2( func0(U(ind1))*func0(V(ind2)) + func1(U(ind1))*func1(V(ind2)) );
    end
    
    if (cnt <= N^kdash),
        % str = randi([1 N],1,kdash);
        str = dec2bigbase(cnt-1,N,kdash) + 1;
        prob = 1;
        for j = 1:kdash,
            prob = prob*PV(indices(str(j)));
        end
        sum2 = 0;
        for cnt2 = 0:2^kdash-1,
            str2 = dec2bigbase(cnt2,2,kdash);
            if (mod(sum(str2),2) == 0),
                prod2 = 1;
                for cnt3 = 1:kdash,
                    if (str2(cnt3) == 0),
                        prod2 = prod2*func0(V(indices(str(cnt3))));
                    else
                        prod2 = prod2*func1(V(indices(str(cnt3))));
                    end
                end
                sum2 = sum2 + prod2;
            end
        end
        sum1 = sum1 + ldash/kdash*prob*log2(sum2);
    end
    
    if (mod(cnt,1000) == 0),
        cnt
        sum1
    end
end
return;
    
        
function f = func0(u)
f = (1+exp(-2*u))^(-1);
return;

function f = func1(u)
f = ((1+exp(-2*u))^(-1))*exp(-2*u);
return;