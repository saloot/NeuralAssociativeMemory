function f = lambda_poly(lambda,x)

f = 0;
for i = 1:length(lambda)
    f = f + lambda(i)*(x.^(i));
end
