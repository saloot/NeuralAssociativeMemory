clc

%G = kasami(8);

%==========================FIND THE CORRELATIONS===========================
S = G(1:500,:);                             % S is the set of M*N patterns that we would like to see the correlation
[M,N] = size(S);                            % find the size of S
[corr,pd,ax] = correlation_patterns(S);     % The function gives us the correlation, its pdf (in pd) and the axis for the pdf (in ax)

% ---------------Find the Most Tilted Pattern Toward Right-----------------
% Among all patterns, the one which is most shifted toward right is the
% worst on in terms of error probability. Here, we find such a pattern to
% compare error probabilities.
nonzero_ind = -1000;
M_ind = 1;
for m = 1:M    
    for j = length(pd(1,:)):-1:1
        if (pd(m,j) ~= 0)
            if (j > nonzero_ind)
                nonzero_ind = j;
                M_ind = m;
            end
        end
    end
end
%--------------------------------------------------------------------------        
pd_ave = sum(pd)/M;                 % Find the average pdf over all M patterns
%==========================================================================
    

%===================COMPARE WITH RANDOM SEQUENCES==========================
[M,N] = size(S);        
rn = randint(M,N);              % Generate binary random sequences
rn = -2*rn + ones(M,N);         % Transform 0/1 into +1/-1
[corr_rnd,pd_rnd,ax_rnd] = correlation_patterns(rn);

% ---------------Find the Most Tilted Pattern Toward Right-----------------
% Among all random patterns, the one which is most shifted toward right is the
% worst on in terms of error probability. Here, we find such a pattern to
% compare error probabilities.
nonzero_ind_rnd = -100000;
M_ind_rnd = 1;
for m = 1:M    
    for j = length(pd_rnd(1,:)):-1:1    
        if (pd_rnd(m,j) ~= 0)            
            if (j > nonzero_ind_rnd)                         
                nonzero_ind_rnd = j;
                M_ind_rnd = m;
            end
        end
    end
end
%--------------------------------------------------------------------------        
pd_rnd_ave = sum(pd_rnd)/M;     % Find the average pdf over all M random patterns    
%==========================================================================


%===========================PLOT THE RESULTS===============================
figure;
plot(ax,pd_ave,'b')                     % Plot average pdf for the patterns
hold on
plot(ax,pd(M_ind,:),'g')                % Plot the most tilted pattern toward right
plot(ax_rnd,pd_rnd(M_ind_rnd,:),'k')    % Plot the average pdf for random patterns
plot(ax_rnd,pd_rnd_ave,'r')             % Plot the most tilted random pattern toward right
title('Comparison of correlation for Carefully Chosen vs. Completely Random Binary Patterns');
xlabel('Correlation');
ylabel('pdf');
legend('Average pdf for all patterns','pdf of the most tilted pattern','Average pdf for all random patterns','pdf of the most tilted random pattern');
%==========================================================================


%========================Calculate Probability of Error====================
th = 1;                                 % th is the integral threshold

%------Calculate the Probability of Error for Average Pattern pdf----------
Pe_pd_ave = 0;                          % This is the probability of error for average pdf
for m = 1:M
    for j = 1:length(pd(1,:))
        if ((ax(j) > th)&(pd(m,j)~=0))
            Pe_pd_ave = Pe_pd_ave + 1;
            break
        end
    end
end
Pe_pd_ave = Pe_pd_ave/(M);
%--------------------------------------------------------------------------


%-----Calculate the Probability of Error for the Worst Pattern's pdf-------
Pe_pd_tilted = 0;                       % This is the probability of error for the worst pattern's pdf
for j = 1:length(pd(1,:))
    if ((ax(j) > th)&(pd(M_ind,j)~=0))
        Pe_pd_tilted = Pe_pd_tilted + 1;%pd(M_ind,j);
        break
    end
end
%Pe_pd_tilted = Pe_pd_tilted/sum(pd(M_ind,:));
%--------------------------------------------------------------------------


%----Calculate the Probability of Error for Average Random Pattern pdf-----
Pe_rnd_ave = 0;                         % This is the probability of error for average pdf of random patterns
for m = 1:M
    for j = 1:length(pd_rnd_ave)
        if ((ax_rnd(j) > th)&(pd_rnd(m,j)~=0))
            Pe_rnd_ave = Pe_rnd_ave + 1;%pd_rnd(m,j);
            break
        end
    end
end
Pe_rnd_ave = Pe_rnd_ave/(M);
%--------------------------------------------------------------------------

%---Calculate the Probability of Error for the Worst Random Pattern's pdf--
Pe_rnd_tilted = 0;                      % This is the probability of error for the worst random pattern's pdf
for j = 1:length(pd_rnd_ave)
    if ((ax_rnd(j) > th)&(pd_rnd(M_ind_rnd,j)~=0))
        Pe_rnd_tilted = Pe_rnd_tilted + 1;%pd_rnd(M_ind_rnd,j);
        break
    end
end
%Pe_rnd_tilted = Pe_rnd_tilted/sum(pd_rnd_ave);
%--------------------------------------------------------------------------

corr_ave = sum(sum(corr))/(M*N);

%----------------------------Display Results-------------------------------
display(['Number of patterns = ',num2str(M)]);
display(['Length of patterns = ',num2str(N)]);
display('-------------------------------------------------')
display(['Average error probability for carefully chosen patterns = ']);
display(num2str(Pe_pd_ave));
display(['Error probability for the most tilted carefully chosen patterns = ']);
display(num2str(Pe_pd_tilted));
display('-------------------------------------------------')
display(['Average error probability for random patterns = ']);
display(num2str(Pe_rnd_ave));
display(['Error probability for the most tilted random patterns = ']);
display(num2str(Pe_rnd_tilted));
display('-------------------------------------------------')
%--------------------------------------------------------------------------

%==========================================================================