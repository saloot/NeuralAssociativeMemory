function img_out = binary_to_image(binary_img_in,quantization_levels)

KK = (1+floor(log(quantization_levels)/log(2))); %ceil(log(quantization_levels)/log(2));

N = length(binary_img_in);
N = round(N/KK);
img_out = bi2de(reshape(binary_img_in,[N,KK]));
%pointer = 1;
%img_out = [];
%while(pointer < N)
%    temp_pixel = [];
%    for i = 0:KK-1
%        temp_pixel = [temp_pixel,(binary_img_in(pointer+i))];
%    end
%    img_out = [img_out,bi2de(temp_pixel)];
%    pointer = pointer + KK ;
%end
111;    