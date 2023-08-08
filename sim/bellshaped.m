function x=BellShapedNoise(f_low,f_up,t,n)
%% Simulate noise from bell-shaped spectrum
%
% Inputs:
% f_low: lower frequency bound (5%)
% f_up: upper frequency bound (95%)
% t: time axis for simulation
% n: number of simulations
%
% Outputs:
% x: matrix with simulated time series as rows
%

%% Check input

if nargin==3
    n=1;
end

if f_low>f_up
   error('f_low must be lower than f_up'); 
end

%% Simulate

% Find maximum omega 
dt_target=diff(t(1:2));
T_target=t(end);

[domegasim,omega_max]=mc_freqaxis(dt_target,T_target);
f_max=omega_max/(2*pi);

% Simulation axis and spectrum
f_axis=[1e-3:1e-3:f_max];
omega_axis=f_axis*2*pi;

f_mean=(f_low+f_up)/2;
omega_mean=f_mean*2*pi;

% S(w)=exp(-0.5*(w-w_mean).^2/c^2);
omega_up=2*pi*f_up;
c=sqrt((-0.5*(omega_up-omega_mean).^2)/log(0.05));

S(1,1,:)=exp(-0.5*(w_axis-w_mean).^2/c^2);

% Simulate by MC
x_sim=mc_sim(omega_axis,S,t,n);

x_sim=x_sim./std(x_sim,0,2);