%MUSIC Simulation
clc
clear all
close all
format long %The data show that as long shaping scientific
hm=zeros(10,3);
for h=1:10
%Transmitter
doa=[45 -30]/180*pi; %Direction of arrival
N=h;%Snapshots
w=[pi/4 pi/6]';%Frequency
M=8;%Number of array elements
P=length(w); %The number of signal
lambda=150;%Wavelength
d=lambda/2;%Element spacing

D=zeros(P,M); %To creat a matrix with P row and M column
for k=1:P
D(k,:)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]); %Assignment matrix
end
D=D';
xx=exp(j*(w*[1:N])); %Simulate message signal
x=D*xx; %Signal after beamforming at transmitter

%Channel
snr=0%SNR
hm(h,1)=N;
x=x+awgn(x,snr);%Insert Gaussian white noise

%Reciever
Sx=x*x'; %Data covarivance matrix
J=fliplr(eye(M)); %Exchange matrix
Sx=Sx+J*conj(Sx)*J;

[nn,vv]=eig(Sx); %Find the eigenvalues and eigenvectors of Sx

%NN=nn(:,1:M-P); %Estimate noise subspace
SS=nn(:,M-P+1:M); %Estimate signal subspace

%Solution for ESPRIT constrainted equations
phi=linsolve(SS(1:M-1,:),SS(2:M,:));
esp_doa=asin(angle(eig(phi))/(2*pi*d/lambda))*180/pi;
esp_doa=round(esp_doa);
esp_doa=sort(esp_doa)
hm(h,2)=min(esp_doa);
hm(h,3)=max(esp_doa);
end
