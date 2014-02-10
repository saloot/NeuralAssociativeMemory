sigma = 1;
x = [-1:.01:1];
f1 = x.*(1-(tanh(sigma*x.^2).^2));
f2 = x;
plot(x,f1,'r')
hold on
sigma = 10;
f1 = x.*(1-(tanh(sigma*x.^2).^2));
plot(x,f1,'g')
sigma = 100;
f1 = x.*(1-(tanh(sigma*x.^2).^2));
plot(x,f1,'b')
plot(x,f2,'k')
legend('w_i(1-tanh^2(w_i^2))','w_i(1-tanh^2(10 w_i^2))','w_i(1-tanh^2(100 w_i^2))','w_i')