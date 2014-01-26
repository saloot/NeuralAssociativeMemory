% this piece of code plots parts of the binomial sum which helps us find
% suitable bounds on the value of the sum based on the number of terms in
% it. 

warning off
U = 100;
a = 10;
summ = zeros(1,U+1);

for j = 0:U

    L = U-j;
    temp = 0;


    for i = L:U
        temp = temp + nchoosek(U,i) * (a^i);% * ((1-a)^(U-i));
    end

    summ(j+1) = temp;
end

figure
L = [U:-1:0];
plot(L,log(summ))
   