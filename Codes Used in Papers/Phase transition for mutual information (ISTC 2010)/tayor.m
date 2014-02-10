n = 5
for i = 1:n
    y=sprintf('%s', 'x', num2str(i));
    syms(y);
end

f = 0;
for i = 1:n
    y=sprintf('%s', 'x', num2str(i));
    f=f+eval(y)^2;
end
f = log(f);

% Find Taylor Coeffecients
for i=1:n
    eval(['x' num2str(i) '=1/n']); 
end



subs(diff(f,x1,4))
