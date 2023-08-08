function [S omega S0]=posdefpsd(S0,omega0,varargin)

%% Modify spectrum so that the eigenvalues are positive

% Small value is added to the diagonal

% Inputs:
% S0: [n,n,N] original spectrum
% omega0: original frequency vector in rad/s

%% Inputs 

p=inputParser;
addParameter(p,'omegamax',[],@isnumeric) % Maximum desired omega, S will be appended with zeros
addParameter(p,'factor_add',1e-3,@isnumeric) % Modification factor (add to diagonal)
parse(p,varargin{:});

omegamax = p.Results.omegamax;
factor_add = p.Results.factor_add;

%%

%Make real diagonal
for i=1:size(S0,1)
    S0(i,i,:)=real(S0(i,i,:));
end

%Add to every diagonal vertically along frequency axis
for i=1:size(S0,1)
    S0(i,i,:)=S0(i,i,:)+max(S0(i,i,:))*factor_add;
end

%Make hermittian symmetric
for i=1:size(S0,3)
    S0(:,:,i)=(S0(:,:,i)+S0(:,:,i)')/2;
end

%Trail zeros at spectrum end
if ~isempty(omegamax)
    domega=omega0(2)-omega0(1);
    omega_add=[domega:domega:omegamax];
    omega=[omega0 omega0(end)+omega_add];
    S=interp1zbig(omega0,S0,omega,'linear',0);
else
   omega=omega0;
   S=S0;
end

end

function [ Mi ] = interp1zbig(z,M,zi,varargin)
%%interp1z Interpolates input data in similar fashion to interp1, but along z-axis (3. dimension)
%
% INPUT:    z:                  original z-axis
%           M:               	original matrix
%           zi:                 interpolated z-axis              
%                            
% OUTPUT:   Mi:                	interpolated matrix
%   
%
% Knut Andreas Kvaale (c) 2013
%

z=z(:);
zi=zi(:);

if nargin==3
    method='linear';
%     extrapolate='extrap';
    extrapolate=0;
elseif nargin==5
    method=varargin{1};
    extrapolate=varargin{2};
end



[Lx,Ly,Lz] = size(M);
Lzi=length(zi);

Mmod=reshape(M,[],Lz).'; 
Mmodi = interp1(z,Mmod,zi,method,extrapolate);
Mi=zeros(Lx,Ly,Lzi);

%calculate range
clear z_range
n_max=3e7;

nbins=ceil(Lzi*Lx*Ly/n_max);
z_range=rangebinequal(Lzi,nbins);

% disp(['Range bins: ' num2str(length(z_range))]);
% 
for k=1:length(z_range)
	range_k=z_range{1,k};
	Mi(1:Lx,1:Ly,range_k)=reshape(Mmodi(range_k,:).',Lx,Ly,length(range_k));
end

end

function range=rangebinequal(n,d)

%% Devide range [1:n] into d almost equal bins

% Inputs:
% n: range end
% d: number elements in bin

% Outputs:
% range: cell with range for each bin

%%

if n==0 | d==0
    error('n or d zero')
end

n_elements=ceil(n/d);

nbins=d;

if nbins==1;
    range{1,1}=1:n;
    return
end

for k=1:d-1
    range{1,k}=[1:n_elements]+(k-1)*n_elements;
end

for k=d
	n_low=range{1,k-1}(end)+1;
	n_up=n;
    range{1,k}=n_low:n_up; 
end

end
