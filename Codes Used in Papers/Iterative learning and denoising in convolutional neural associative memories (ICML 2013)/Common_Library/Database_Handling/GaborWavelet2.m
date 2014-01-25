function phi = GaborWavelet2 (width, height, K, a,b, theta,P)
% Create the Gabor Wavelet Filter
% Author : Chai Zhi  
% e-mail : zh_chai@yahoo.cn



phi = zeros ( width , height );
w_0 = theta;
F_0 = 1;

for x = -width/2 + 1 : width/2
    
    for y = -height/2 + 1 : height/2
        
        x_r = x*cos(theta) + y*sin(theta);
        y_r = -x*sin(theta) + y*cos(theta);
        s = exp(1i*(2*pi*F_0*(x*cos(w_0)+y*sin(w_0))+P));
        phi(x+width/2,y+height/2) = K*exp(-pi*( (a*x_r)^2 + (b*y_r)^2)) *s;
%         phi(x+width/2,y+height/2) = ( f^2/(pi*gamma*eta) ) * exp( f^2*((x_r/gamma)^2-(y_r/eta)^2) ) *  exp( i*2*pi*f*x_r);
    
    end

end
