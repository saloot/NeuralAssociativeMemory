function Pe = p_1_2_j(x,epsilon,gamma,dj)


if (x<=(gamma-epsilon)*dj)
    temp1 = 0;
elseif ( x>=(gamma+epsilon)*dj ) 
    temp1 = 1;
else
    if (epsilon > 0)
        temp1 = (epsilon-gamma+x/dj)/(2*epsilon);
    else
        temp1 = 0;
    end
end

if (x>=(epsilon-gamma)*dj)
    temp2 = 0;
elseif ( x<=-(epsilon+gamma)*dj)
    temp2 = 1;
else
    if (epsilon > 0)
        temp2  = (epsilon-gamma-x/dj)/(2*epsilon);
    else
        temp2 = 0;
    end
end

Pe = temp1 + temp2;
    