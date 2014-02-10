% Code generates (2^m-1)(2^m+1) sequences of length 2^m-1, whose maximum 
% magnitude of correlation is 1+sqrt(2^(m+1)). Sequences are over the {+1,-1}
% alphabet, and are known to be correlation optimal for binary sequences 
% families of this particular size and length
% -------------------------


% Initialization
m = 9; % m has to be odd
l = 2; % Choose any l such that gcd(l,m)=1
% -------------------------


if (gcd(l,m)~=1),
    error('GCD of l and m is not equal to one!');
end

d = 2^l+1;
alph = gf([2],m);
G = [];

for i = 0:2^m-2,
    i
    g = [];
    for t = 1:2^m-1,
        g = [g gftrace((alph^i)*alph^t+(alph^(d*t)),m)];
    end
    G = [G; g];
end

% g = [];
% for t = 1:2^m-1,
%     g = [g gftrace(alph^(d*t),m)];
% end
% G = [G; g];
% 
% g = [];
% for t = 1:2^m-1,
%     g = [g gftrace(alph^t,m)];
% end
% G = [G; g];
% 
% G1 = G;
% n = 2^m-1;
% for i = 1:n-1,
%     G = [G; circshift(G1,[0,i])];
% end

G = (-1).^double(G.x);

% corrchk = G*(G') - n*eye(n*(2^m+1));
% max(max(abs(corrchk)))