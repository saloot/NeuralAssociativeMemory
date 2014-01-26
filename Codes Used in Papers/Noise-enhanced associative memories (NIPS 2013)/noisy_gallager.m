function [BER,itr] = noisy_gallager(n,m,d_v,d_c,upsilon,nu,varphi,itr_max,p_e,no_tested_instances)

success_flag = 0;
message_var_to_check = zeros(m,n);
message_check_to_var = zeros(m,n);



BER = 0;

H = bipartite(n,m,d_v,d_c);

for k = 1:no_tested_instances

    channel_output = zeros(1,n);
nois = rand(1,n);
nois = (nois  < p_e);    
channel_output = mod(channel_output + nois,2);
for i = 1:m
    for j = 1:n
        if (H(i,j) > 0)
            message_var_to_check(i,j) = channel_output(j);
        end
    end
end



for itr = 1:itr_max
    
    %===================Check-to-Variable Nodes Messages===================
    for i = 1:m
        for j = 1:n
            if (H(i,j) > 0)
                message_check_to_var(i,j) = abs(mod(sum(message_var_to_check(i,:))-message_var_to_check(i,j),2));
                
                %-------------Add the Effect of Internal Noise-------------
                p = rand;
                if (p < nu)
                    message_check_to_var(i,j) = mod(1+message_check_to_var(i,j),2);
                end
                %----------------------------------------------------------
                
            end
        end
    end
    %======================================================================
    
    %===================Check-to-Variable Nodes Messages===================
    for j = 1:n
        for i = 1:m
            if (H(i,j) > 0)
                normalized_feedback = sum(message_check_to_var(:,j))-message_check_to_var(i,j);
                update_decision = normalized_feedback/sum(abs(sign(H(:,j)))) + upsilon*(rand-0.5);
                
                if (update_decision > varphi )
                    message_var_to_check(i,j) = 1;
                elseif (update_decision < 1 - varphi )
                    message_var_to_check(i,j) = 0;
                else
                    message_var_to_check(i,j) = channel_output(j);
                end
            end
        end
    end
    if (norm(message_var_to_check) == 0)
        success_flag = 1;
        break;
    end
    
    %======================================================================
    
end
BER = BER + sum(sign(sum(message_var_to_check)));
end

BER = BER/n/k;
