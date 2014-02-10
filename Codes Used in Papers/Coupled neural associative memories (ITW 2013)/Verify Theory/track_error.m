addpath(genpath('/home/amir/cluster/Common_Library'));
n = 4096
e = 2;
p_e_total = [];
p_e_range = [0:.01:.31];
for iji = 1:length(p_e_range)
    p_e = p_e_range(iji);
    
L_horiz = 29;
Omega = 2;
try_max = 10*L_horiz;
pi = zeros(1,L_horiz);
z = p_e *ones(1,L_horiz);
d = 64;
end_flag = 0 ;
try_itr = 0;
zz = [];
while (end_flag == 0)
for l = 1:L_horiz
    z_bar = 0;
    av_count = 0;
    for i = -Omega:Omega
        if ((l-i>0)&&(l-i<L_horiz))
            z_bar = z_bar+z(l-i);
            av_count= av_count + 1;
        end
    end
    z_bar = z_bar/av_count;
    if (e == 1)
        pi(l) = 1-((1-z_bar)^(d-1));
    elseif (e == 2)
        pi(l) = 1-((1-z_bar)^(d-1))-z_bar*(d-1)*(1-z_bar)^(d-2);
    else
        error('e non-supported');
    end
        
end
for l = 1:L_horiz
    
    pi_bar = 0;
    av_count = 0;
    for i = -Omega:Omega
        if ((l-i>0)&&(l-i<L_horiz))
            pi_bar = pi_bar+pi(l-i);
            av_count= av_count + 1;
        end
    end
    pi_bar = pi_bar/av_count;
    z(l) = p_e *lambda_poly_v2(pi_bar,lambda(l,:));
end


try_itr = try_itr + 1;
if ((try_itr > try_max)||(norm(z) == 0))
    end_flag = 1;
end
zz = [zz,max(mean(z),0)];
% z
111;
end
% plot(zz)
mean(z)
p_e_total = [p_e_total,1-(1-max(mean(z),0))^n];
end
% plot(1-(1-zz).^64)

hold on

plot(p_e_range,p_e_total,'r');
