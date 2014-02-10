%=================Function: correlation_patterns===========================

%================================INPUT=====================================
% S: is the set of M binary patterns with length N stored in an M*N matrix.
%==========================================================================

%===============================OUTPUT=====================================
% coor: is the correlation term between each pattern and all other patterns for each bit.
% pd: is the pdf of correlations for each pattern
% ax: is the range of the pdf
%==========================================================================

%==========================================================================

%=============================DESCRIPTION==================================
% This function gets a set of patterns as input and calculates the
% correlation term between each pattern and all other ones based on the
% equation 2.13 of the book "Introduction to the Theory of Neural
% Computation" by Hertz et al. (the blue book!)
%==========================================================================


function [corr,pd,ax] = correlation_patterns(S)
[M,N] = size(S);                        

%---------------------------Compute the Correlation------------------------
corr = zeros(M,N);
for nu = 1:M
    for i = 1:N
        summ = 0;
        for mu = 1:M
            if (mu ~= nu)
                summ = summ + S(mu,i)*(S(mu,:)*(S(nu,:))');
            end
        end
        corr(nu,i) = -S(nu,i)*summ/N;
    end
end
%--------------------------------------------------------------------------


%-------------------------Calculate the pdf--------------------------------
pd_length = 100;                % Divide the correlation range into bins (the number of bins = pd_length)
step_size = ((max(max(corr))-min(min(corr))))/pd_length;    % Calculate the step size
pd = zeros(M,pd_length);        % Initalize the pdf

for m = 1:M
    for i = 1:pd_length
        left = min(min(corr)) + (i-1)*step_size;
        right = min(min(corr)) + (i)*step_size;
        for j = 1:N
            if ((corr(m,j)>=left)&(corr(m,j)<right))
                pd(m,i) = pd(m,i)+1;
            end
        end
    end
end
if (step_size ~=0)              % Take care of a scpecial case when all patterns have the same correlation
    ax = [min(min(corr)):step_size:max(max(corr))];
else
    for m = 1:M
        pd(m,50) = 1;
    end
    ax = [min(min(corr))-50:1:min(min(corr))+49];
end

ax = ax(1:pd_length);
%--------------------------------------------------------------------------
