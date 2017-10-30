%% Model for LER in SNW-FET
clear all;
%% Change parameters here
radius = 5;         % NW radius
l = 10;             % NW length
A = 1/3;              % Roughness amplitude
lambda = 2;         % Correlation length
alpha = 1;          % Roughness coefficient
corr = 1;           % Correlation cofactor
Fs = 8;             % Sampling frequency
Fr = 64;            % Samples along circumference
%%
sp = -l/2:l/Fs:l/2-1/Fs;
theta = 0:2*pi/Fr:2*pi-1/Fr;
xm = [5*cos(theta)' 5*sin(theta)'];
z = (0:2*Fs-2)*l/(Fs-1)        % CHECK HERE!

%acf = A*gaussmf(x,[lambda^2,0]);
acf = A*exp((-sp.^2/(2*lambda^2))*alpha);
n = 2^nextpow2(length(sp));
facf = fft(acf,n);
facf = facf/n;
%f = Fs*(-n/2:n/2-1)/n;

%facf = abs(facf/n);
%facf1 = facf(1:n/2+1);
%facf1 (2:end-1) = 2*facf1(2:end-1);

%plot(f,facf);

%corr_mat = [1 0; corr sqrt(1-corr^2)];
 
noise1 = wgn(1,Fs,lambda^2);
% noise2 = wgn(1,Fs,lambda^2);
% pnoise = corr_mat*[noise1 ; noise2];
noise = wgn(Fs,Fr,10*log10(lambda^2));

fnoise = fft(noise)/n;

for i = 1:Fr
    a = conv(fnoise(:,i),facf);
    a = n*ifft(a);
    frx = real(a)*cos(theta(i)) + xm(i,1);
    fry = real(a)*sin(theta(i)) + xm(i,2);
    
    plot3(frx,fry,z,'.b');
    hold on
    axis([-2*radius 2*radius -2*radius 2*radius 0 2*l])
    
end



