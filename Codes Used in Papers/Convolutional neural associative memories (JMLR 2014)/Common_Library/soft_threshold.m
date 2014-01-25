function W_out = soft_threshold(W_in,theta)
W_out = zeros(length(W_in),1);
for i = 1:length(W_in)
    if (abs(W_in(i)) > theta)
        W_out(i) = W_in(i);%-theta;
    end
end
        