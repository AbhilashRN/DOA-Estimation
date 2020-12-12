%MVDR Simulation
clc
clear all
format long %The data show that as long shaping scientific

%Transmitter
doa=[-40 20; 70 -50]/180*pi; %Direction of arrival
N=100;%Snapshots
w=[1 pi/4]';%Frequency
M=8;%Number of array elements Y axis
NX=8;
P=length(w); %The number of signal
lambda=150;%Wavelength
d=lambda/2;%Element spacing
snr=0;%SNR
D=zeros(P,M,NX); %To creat a matrix with P row and M column
for k=1:P
    for l=0:M-1
        D(k,l+1,:)=exp(j*(pi*sin(doa(1,k))*cos(doa(2,k))*l+ pi*sin(doa(1,k))*sin(doa(2,k))*[0:NX-1])); %Assignment matrix
    end
end
D=reshape(D,P,M*NX);
xx=2*exp(j*(w*[1:N])); %Simulate message signal
x=D'*xx; %Signal after beamforming at transmitter

%Channel
x=x+awgn(x,snr);%Insert Gaussian white noise

%Reciever
Sx=x*x'; %Data covarivance matrix

%Creating array manifold vector V
theta=-90:1:90; %Peak search
phi = -90:1:90;
Pmvdr=zeros(length(theta),length(phi));
for ii=1:length(theta)
    for iii =  1:length(phi)
        V=zeros(length(M),length(NX));%array manifold vector
        for jj=-(M-1)/2 : 1 : (M-1)/2
            for jjj = -(NX-1)/2 : 1 : (NX-1)/2
            V(1+((M-1)/2)+jj,1+(NX-1)/2 +jjj)=exp(j*(jj*pi*sin(theta(ii)/180*pi)*cos(phi(iii)/180*pi)+...
                jjj*pi*sin(theta(ii)/180*pi)*sin(phi(iii)/180*pi)));
            Vnew(1+((M-1)/2)+jj,1+(NX-1)/2 +jjj)= exp(j*(jj*pi*sin(50/180*pi)*sin(120)/180*pi)+...
                jjj*pi*sin(50/180*pi)*cos(120/180*pi));
            end
        end
        V=reshape(V,1,M*NX);
        PP=(V/Sx)*(V)';
        Pmvdr(ii,iii)=abs(1/ PP);
    end

end
Pmvdr=10*log10(Pmvdr/max(max(Pmvdr))); %Spatial spectrum function
mesh(phi,theta,Pmvdr)
xlabel('azimuthal angle \phi/degree')
ylabel('elevation angle \theta/degree')
zlabel('spectrum function P(\theta,\phi) /dB')
title('DOA estimation based on planar MVDR algorithm ')


