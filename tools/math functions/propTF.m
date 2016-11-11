function u2 = propTF(u1,L,lambda,z)
% propagation - transfer fucntion approach
% assume the same x and y side lengths and 
% uniform sampling
% u1 - source plaine field
% L - source and observation plane side length
% lambda - wavelength
% z - porpagation distance
% u2 - observation plane field

[M,~] = size(u1);           % get input field array size
dx = L/M;                   % sample interval
% k = 2*pi/lambda;            % wavenumber

fx = -1/(2*dx):1/L:1/(2*dx)-1/L;    % frequence coords
[FX,FY] = meshgrid(fx,fx);

H = exp(-1i*pi*lambda*z*(FX.^2+FY.^2));  % transfer function
H = fftshift(H);
U1 = fft2(fftshift(u1));
U2 = H.*U1;
u2 = ifftshift(ifft2(U2));
end


