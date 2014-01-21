function WAV = freq2wav(V,N_out,duration)

z = zeros(1,duration);
step_size = floor((floor(length(z)/2)+1)/(N_out));
% z(1) = V(1);
for i = 1:N_out
    z((i-1)*step_size+1) = V(i)+.00025;
    z(duration-(i-1)*step_size-1) = V(i)+.00025;
end
z = z*duration;
WAV = ifft(z);

