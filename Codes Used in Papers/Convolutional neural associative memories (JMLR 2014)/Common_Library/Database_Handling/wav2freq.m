function V = wav2freq(WAV,N_out,threshold)

V = zeros(1,N_out);
z = fft(WAV)/length(WAV);
% z = z/sum(abs(z));
step_size = floor( (floor(length(z)/2)+1)/(N_out+1));
% ind = [1:step_size:floor(length(z)/2)+1-step_size];
% ind = [ind,floor(length(z)/2)+1];
ind = linspace(1,1+(N_out)*floor((floor(length(z)/2)+1)/(N_out)),N_out+1);
ind(N_out+1) = floor(length(z)/2)+1;
b = sum(abs(z))/length(z);

for i = 2:N_out+1
    deltax = (ind(i)-ind(i-1))/(floor(length(z)/2)+1);
%     s = sum(abs(z(ind(i-1):ind(i))))/deltax;
    s = max((abs(z(ind(i-1):ind(i)))));
    V(i-1) = s.*(s>=threshold);%*(floor(length(z)/2)+1); 
end
% V = V/max(V);
% V = V*0.15/threshold;

       