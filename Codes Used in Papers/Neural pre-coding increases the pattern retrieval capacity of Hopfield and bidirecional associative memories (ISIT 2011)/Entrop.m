figure;
hold on;
for x1 = 0:.001:1
    for x2 = 0:.001:1
        if (x1+x2 < 1)
            H = -x1*log(x1) - x2*log(x2) - (1-x1-x2)*log(1-x1-x2);
            plot3(x1,x2,H);
        end
    end
end
