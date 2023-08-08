%%
clc
clear all
close all

%% Define cross spectral density matrix S
wmax=100;
dw=0.05

w=[dw:dw:wmax];
nw=length(w);

S(1,1,1:nw)=1;
S(2,2,1:nw)=10./w;

figure(); hold on; grid on;
plot(w,squeeze(S(1,1,:)),'b');
plot(w,squeeze(S(2,2,:)),'r');
ylim([0 3]);
legend({'S_{11}' 'S_{22}'});

%% Simulate

sim=MCCholeskyFast(w,S,'domegasim',0.001,'Nsim',10);

t=sim{1};
x=sim{2};

% Plot one of the simulations
figure(); hold on; grid on;
plot(t,x(1,:),'b');
plot(t,x(2,:),'r');
legend({'x_{1}' 'x_{2}'});