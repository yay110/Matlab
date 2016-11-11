function out = jinc(x)
%
% jinc function
%
% J1(2*pi*x)/x;
% divided by zero fix
%
% located non-zero elements of x
mask = (x~=0);
% initialize ouput with pi (value for x = 0)
out = pi*ones(size(x));
% compute output values for all other x
out(mask) = besselj(1,2*pi*x(mask))./x(mask);
