% This function generates a full rank random n*m binary matrix A
% in which each column has deg_column 1's. 
% m and n are updated in the end since identical columns are omitted
% during the process
function A = random_matrix(m,n,deg_column)
A = zeros(m,n);
for j = 1:n
    p = randperm(m);  % Generate a permutation of integers from 1 to n.
    for i = 1:deg_column
        A(p(i),j) = 1;      % Here we put integer elements in the generator matrix randomly.
    end    
end

for j = 1:n
    for i = j+1:n
        if (A(:,j) == A(:,i))   % Check for columns with equal values.
            temp_flag = 0;
            while (temp_flag == 0)
                temp_flag = 1;
                A(:,i) = zeros(m,1);    
                p = randperm(m);  % Generate a permutation of integers from 1 to n.
                for ii = 1:deg_column
                    A(p(ii),i) = 1;      % Here we put integer elements in the generator matrix randomly.
                end  
                for k = 1:j       % Remove duplicate columns by shifting the columns to the left.
                    if (A(:,i) == A(:,k))
                        temp_flag = 0;
                    end
                end
            end
            
        end
    end
end

% % AA = A(:,1:n);
% A = AA;