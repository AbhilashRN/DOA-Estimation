%MUSIC Simulation
% Cohernce Problem
clc
clear all
close all
format long %The data show that as long shaping scientific

for kk = 1:6
    
%Transmitter
doa=[-30 30]/180*pi; %Direction of arrival
%N=1*ceil(5*(2^(kk/1.3))/10);%Snapshots
N=200;
w=[pi/4 pi/4+0.001*kk-0.003]';%Frequency
M=8;%Number of array elements
P=length(w); %The number of signal
lambda=150;%Wavelength
d=lambda*0.5;%Element spacing

D=zeros(P,M); %To creat a matrix with P row and M column
for k=1:P
D(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]); %Assignment matrix
end
D=D';
xx=2*exp(j*(w*[1:N])); %Simulate message signal
x=D*xx; %Signal after beamforming at transmitter

%Channel
snr=0;%SNR
x=x+awgn(x,snr);%Insert Gaussian white noise

%Reciever
Sx=x*x'; %Data covarivance matrix
J=fliplr(eye(M)); %Exchange matrix


Sx=Sx+J*conj(Sx)*J;
%Sx=Sx+J*conj(Sx)*J;%IMPROVED MUSIC

[nn,vv]=eig(Sx); %Find the eigenvalues and eigenvectors of Sx
NN=nn(:,1:M-P); %Estimate noise subspace
%Creating array manifold vector V
theta=-90:0.5:90; %Peak search
for ii=1:length(theta)
V=zeros(1,length(M));%array manifold vector
for jj=0:M-1
V(1+jj)=exp(-j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
end

%Solution for MUSIC constrainted equations
PP=V*(NN*NN')*V';
Pmusic(ii)=abs(1/ PP);
end
Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function

subplot(2,3,kk)
plot(theta,Pmusic,'-k')
xlabel('angle \theta/degree')
ylabel('spectrum function P(\theta) /dB')
title(['dw= ',num2str(round(0.001*kk-0.003,6)),' '])
grid on

end
 %sgtitle('MUSIC-COHERENCE PROBLEM');
 sgtitle('IMPROVED MUSIC-COHERENCE');