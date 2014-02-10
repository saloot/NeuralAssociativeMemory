function montanari_bound_ldpc_BSC()

Nsimul = 40000;    % # of times monte-carlo sampling is performed in evaluating Montanari's bound

% DE parameters
del = .25;
Ndash = 100;
dv = [3;1];
dc = [5;1];
max_itr = 1000;

ldash = 3;
kdash = 5;

phival = [];
for p = 0.13:.01:.15,
    
    
    [Perr, PV, PU] = DE_irregLDPC_BSC(p,dv,dc,del,Ndash,max_itr);    
    N = length(PU);
    U = 0.5*[del*(-Ndash:Ndash) Inf];   % The 0.5 up-front is to compensate for Montanari's strange definition of LLRs
    V = U;
    
%     [PU IX] = sort(PU,'descend');   % This sorts the PDFs so that the "significant" values come up-front
%     U = U(IX);
%     [PV IX] = sort(PV,'descend');
%     V = V(IX);
    
    temp = phicalc(p,ldash,kdash,U,PU,V,PV,Nsimul);
    phival = [phival temp];
%     if (temp > 0),
%         break;
%     end
end
phival
% fprintf('phival = %f',phival);
% fprintf('BSC threshold is %f',p);
return;

function sum1 = phicalc(p,ldash,kdash,U,PU,V,PV,Nsimul)

neglect_val = -1;

sum1 = 0;
indices1 = find(PU > neglect_val);  % This part of the code neglects PDF values that are smaller than neglect_val
indices2 = find(PV > neglect_val);
indices = union(indices1,indices2);
N = length(indices);
if (N == 1),    % The code below doesn't work if only one point is chosen, so I include one more in this case
    indices = [indices 2];
    N = 2;
end

fprintf('\nTotal number of significant points in the pdf = %d',N);
fprintf('\nIterations needed for an exact result = %d',N^max(ldash,kdash));
fprintf('\nMaximum number of Monte-Carlo iterations = %d',Nsimul);
fprintf('\nFraction of probability mass in U, V being taken into account = %f',(sum(PV(indices))+sum(PV(indices)))/2);
fprintf('Press any key to continue...');
% pause;


    
    for cnt = 1:N^2,        % First term of (6.2) in Montanari's paper
        %str = randi([1 N],1,2);
        str = dec2bigbase(cnt-1,N,2) + 1;
        ind1 = indices(str(1));
        ind2 = indices(str(2));
        prob = PU(ind1)*PV(ind2);
        sum1 = sum1 - ldash*prob*log2( func0(U(ind1))*func0(V(ind2)) + func1(U(ind1))*func1(V(ind2)) );
    end

    for cnt = 1:Nsimul
%     if (cnt <= N^ldash),     % Second term of (6.2) in Montanari's paper
        % str = randi([1 N],1,ldash);
%         str = dec2bigbase(cnt-1,N,ldash) + 1;

        index_a = zeros(1,ldash);
        for iii = 1:ldash
            rando = rand;
            ppp = 0;
            kkk = 1;        
            while (rando>ppp)
                ppp = ppp + PU(kkk);
                kkk = kkk +1;
            end
            index_a(iii) = kkk-1;
        end
        prod0 = 1;
        prod1 = 1;
        prob = 1;
        for j = 1:ldash,
            prod0 = prod0*func0(U(index_a(j)));
            prod1 = prod1*func1(U(index_a(j)));
            prob = prob*PU(index_a(j));
        end
        sum1 = sum1 + (1-p)*prob*log2(prod0 + p/(1-p)*prod1) + p*prob*log2(prod0 + (1-p)/p*prod1);
%     end
    
    
%     if (cnt <= N^kdash),    % Third term of (6.2) in Montanari's paper
        % str = randi([1 N],1,kdash);
%         str = dec2bigbase(cnt-1,N,kdash) + 1;
        index_b = zeros(1,kdash);
        for iii = 1:kdash
            rando = rand;
            ppp = 0;
            kkk = 1;        
            while (rando>ppp)
                ppp = ppp + PV(kkk);
                kkk = kkk +1;
            end
            index_b(iii) = kkk-1;
        end
        prob = 1;
        for j = 1:kdash,
            prob = prob*PV(index_b(j));
        end
        sum2 = 0;
        for cnt2 = 0:2^kdash-1,
            str2 = dec2bigbase(cnt2,2,kdash);
            if (mod(sum(str2),2) == 0),
                prod2 = 1;
                for cnt3 = 1:kdash,
                    if (str2(cnt3) == 0),
                        prod2 = prod2*func0(V(index_b(cnt3)));
                    else
                        prod2 = prod2*func1(V(index_b(cnt3)));
                    end
                end
                sum2 = sum2 + prod2;
            end
        end
        sum1 = sum1 + ldash/kdash*prob*log2(sum2);
%     end
    
    if (mod(cnt,1000) == 0),
        cnt
        sum1
        % indices(str(1))
    end
end
return;
    
        
function f = func0(u)
f = (1+exp(-2*u))^(-1);
return;

function f = func1(u)
f = ((1+exp(-2*u))^(-1))*exp(-2*u);
return;