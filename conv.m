%Conventional Beamformers
clc; clear; close all;

M      = 16;        % Number of Array Elements.
lambda = 1;         % Incoming Signal Wavelength in (m).
d      = lambda/2;  % Interelement Distance in (m).
phi_s  = -30;        % Steering Angle of Beamformer in (deg).

u_s = (d/lambda)*sin(phi_s*pi/180);  % Normalized Spatial Frequency of the signal of interest.

c = exp(-1i*2*pi*u_s*(0:M-1).')/sqrt(M);

angle = -90:0.1:90;
L = length(angle);
C1 = zeros(1,L);

for k=1:L
    u = (d/lambda)*sin(angle(k)*pi/180);
    v = exp(-1i*2*pi*u*(0:M-1).')/sqrt(M); % Azimuth Scanning Steering Vector.
    C1(k) = c'*v;
end

%figure('NumberTitle', 'off','Name','Figure 11.9','Position',[0 0 600 850]);
subplot(1,2,1);
% This plots the instantaneous power for every element (M waveforms).
plot(angle,10*log10(abs(C1).^2));
ylim([-70 5]);
xlim([-90 90]);
grid on;
title('Beampattern of Conentional Beamformer');
xlabel('Sample Number');
ylabel('Output Power (dB)');

subplot(1,2,2);
polardb(angle*pi/180,10*log10(abs(C1).^2),-60,'b');
title('Beampattern (polar plot) of Conventional Beamformer');
grid on;
