x = [0.001:.05:1];
y = x;
entr = zeros(length(x),length(x));
cst = zeros(length(x),length(x));

for i = 1:length(x)
    for j = 1:length(y)
        entr(i,j) = x(i)*log(x(i))-y(j)*log(y(j));
        cst(i,j) = exp(y(j)-x(i));
    end
end
entr = entr/log(2);

surf(x,y,entr-cst);
grid on
axis square

xlabel('x');
ylabel('y');

