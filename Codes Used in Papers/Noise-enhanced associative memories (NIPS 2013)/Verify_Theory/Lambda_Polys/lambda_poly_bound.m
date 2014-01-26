function f = lambda_poly_bound(lambda,P_correct_cluster,x)

P_c = sort(P_correct_cluster,'ascend');

f = 0;
for i = 1:length(lambda)
    if (lambda(i) == 0)
        continue
    end
    temp = 1;
    for j = 1:i
        temp = temp*(1-P_c(j)*x);
    end
    
    f = f + lambda(i)*temp;
end
