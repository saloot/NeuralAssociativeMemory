function [tr] = gftrace(x,m)

tr = gf([0],m);

for i = 0:m-1,
    tr = tr + x^(2^i);
end

%arr = 2.^[0:m-1];
%tr = sum(x.^arr);

return;