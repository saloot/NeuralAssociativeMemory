function [G] = kasami(m)

% Code generates (2^m-1)(2^(m/2)) sequences of length 2^m-1, whose maximum
% magnitude of correlation is sqrt(2^(m/2))+1 (number of cyclically 
% distinct sequences is 2^(m/2)). Sequences are over the {+1,-1} alphabet.
% -------------------------
%
% Initialization:
% m is any even positive integer.
% -------------------------

if (mod(m,2)~=0),
    error('m has to be even!');
end

k = m/2;
d = 2^k+1;
alph = gf([2],m);
G = [];

for i = 0:2^m-2,
    
    b = alph^i;
    if (b^(2^k) == b),
        g = [];
        for t = 1:2^m-1,
            
            sumgf = gf([0],m);
            sumgf = sumgf + gftrace(alph^t,m);
            for i2 = 0:k-1,
                sumgf = sumgf + (b*alph^(d*t))^(2^i2);
            end
            g = [g sumgf];
            
        end
        G = [G; g];
    end
    
end

g = [];
for t = 1:2^m-1,
    g = [g gftrace(alph^t,m)];
end
G = [G; g];


G1 = G;
n = 2^m-1;
for i = 1:n-1,
    G = [G; circshift(G1,[0,i])];
end

G = (-1).^double(G.x);

% corrchk = G*(G') - n*eye(n*(2^(m/2)));
% max(max(abs(corrchk)))
return;