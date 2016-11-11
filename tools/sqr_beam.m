% sqr_ beam propagtation example
%
L1 = 0.01;                   % side length
M = 5000;
dx1 = L1/M;
x1 = -L1/2:dx1:L1/2-dx1;
y1 = x1;

lambda = 0.5e-6;
k = 2*pi/lambda;
w = 0.002;                  % source half width (m)
z = 0.1;                   % propagation distance (m)
critical_dx = lambda*z/L1;

[X1,Y1] = meshgrid(x1,y1);
% u1 = rect(X1/(2*w)).*rect(Y1/(2*w));  % src field of rect 
u1 = exp(-(X1.^2+Y1.^2)/w^2);           % Gaussian profile of scr field
% u1 = rect(sqrt(X1.^2+Y1.^2)/(2*w));           % circular aperture
I1 = abs(u1.^2);                    % src irradiance

figure(1)
subplot(2,2,1);imagesc(x1,y1,I1);
axis image;axis xy;
xlabel('x (m)');ylabel('y (m)');
title('Irradiance at z = 0 m');

% deg = pi/180;
% alpha = 5.0e-5;     %rad
% theta = 135*deg;
% u1 = tilt(u1,L1,lambda,alpha,theta);

% ============== adding a focus term =============
zf = z;
u1 = focus(u1,L1,lambda,zf);

% ============== adding a tilt term  ===============
% subplot(2,2,3);imagesc(x1,y1,abs(u1.^2));
% axis image;axis xy;
% xlabel('x (m)');ylabel('y (m)');
% title('Irradiance at z = 0 m with tilt');


u2 = propTF(u1,L1,lambda,z);    %propagation by transfer func tion;
u3 = propIR(u1,L1,lambda,z);    %propagation by impluse response;

x2 = x1;                % obs coords
y2 = y1;            
I2 = abs(u2.^2);        % obs irradiance
I3 = abs(u3.^2);        % obs irradiance

subplot(2,2,2);imagesc(x2,y2,I2);
axis image; axis xy;
xlabel('x (m)');ylabel('y (m)');
title(['Irradiance at z = ', num2str(z), ' m by TF']);

subplot(2,2,3);plot(x1,I2(:,M/2+1));
xlabel('x (m)');ylabel('irradiance');
title('Profile of the beam');


subplot(2,2,4);imagesc(x2,y2,I3);
axis image; axis xy;
xlabel('x (m)');ylabel('y (m)');
title(['Irradiance at z = ', num2str(z), ' m by IR']);
