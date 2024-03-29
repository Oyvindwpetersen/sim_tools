function varargout=MCCholeskyFast(omegaaxisinput,SS,varargin)
%% Description
% This function simulates realizations (time series) using the
% cross-spectral density matrix SS(Ndof, Ndof, Nomegaaxisinput) as basis. To improve the performance the
% factorized cross spectral density is interpolated to a denser frequency
% axis befor the time series are generated by means of ifft.
%
% Ole Oiseth 25.10.2014
% OWP 2023

%% Define parameters
p=inputParser;
addParameter(p,'Nsim',1,@isnumeric);                     %Number of time series
addParameter(p,'domegasim',0.001,@isnumeric);           %delta omega used in the simuilations
addParameter(p,'omegamax_chol',1e6,@isnumeric);                  %File Name. If no file name is given, the data will be exported as a cell array.
addParameter(p,'appendomega',[],@isnumeric);                  %File Name. If no file name is given, the data will be exported as a cell array.

parse(p,varargin{:})
Nsim = p.Results.Nsim;
domegasim = p.Results.domegasim;
omegamax_chol = p.Results.omegamax_chol;
appendomega = p.Results.appendomega;

%% Cholesky decomposition of the cross spectral density matrix

[s1,s2]=size(omegaaxisinput);  
if s1>s2; omegaaxisinput=omegaaxisinput.'; end

if ~isempty(appendomega)
    domegainput=omegaaxisinput(2)-omegaaxisinput(1);
    appendomegaaxis=[(omegaaxisinput(end)+domegainput):domegainput:appendomega];
    
    omegaaxisinput=[omegaaxisinput appendomegaaxis];
    SS=cat(3,SS,zeros(size(SS,1),size(SS,2),length(appendomegaaxis)));
end

GG=zeros(size(SS));

treshold_S=max(max(max(abs(SS))))*1e-16;

for n=1:length(omegaaxisinput)
    if  omegaaxisinput(n) > omegamax_chol | ~any(any(SS(:,:,n))) | max(max(abs(SS(:,:,n))))<treshold_S
    else
        GG(:,:,n)=chol(SS(:,:,n),'lower');
    end
end


%% Simulate time series
omegaaxissim=domegasim:domegasim:max(omegaaxisinput);
NFFT=2^nextpow2(2*length(omegaaxissim));
t=linspace(0,2*pi/domegasim,NFFT);
XX=cell(1,Nsim+1);
XX{1,1}=t;

c_precalc3=cell(1,size(GG,1));
for m=1:size(GG,1)
    c_precalc3{m}=zeros(m,length(omegaaxissim));
	for n=1:m    
	c_precalc3{m}(n,:)=interp1(omegaaxisinput,permute(GG(m,n,:),[1,3,2]),omegaaxissim,'linear',0);
	end
end

for z=1:Nsim
    phi=2*pi*rand(size(SS,1),length(omegaaxissim));
    x=zeros(size(SS,1),NFFT);
    for m=1:size(GG,1)
        c_all=c_precalc3{m}.*exp(1i.*phi(1:m,:));
        x(m,:)=sum(real(ifft(c_all,NFFT,2)),1)*NFFT*sqrt(2*domegasim);
    end
	XX{1,z+1}=x;
end

if nargout==1
    varargout{1}=XX;
elseif nargout==2
    varargout{1}=XX{1};
    varargout{2}=XX(2:end);
end

%% Simulate time series
% omegaaxissim=domegasim:domegasim:max(omegaaxisinput);
% NFFT=2^nextpow2(2*length(omegaaxissim));
% t=linspace(0,2*pi/domegasim,NFFT);
% XX=cell(1,Nsim+1);
% XX{1,1}=t;
% 
% c_precalc=cell(size(GG,1),size(GG,1));
% for m=1:size(GG,1)
% 	for n=1:m
% 	c_precalc{m,n}=interp1(omegaaxisinput,permute((GG(m,n,:)),[1,3,2]),omegaaxissim,'linear','extrap');
% 	end
% end
% 
% for z=1:Nsim
%     phi=2*pi*rand(size(SS,1),length(omegaaxissim));
%     x=zeros(size(SS,1),NFFT);
%     for m=1:size(GG,1)
%         for n=1:m
%             c=c_precalc{m,n}.*exp(1i.*phi(n,:));
%             x(m,:)=x(m,:)+real(ifft(c,NFFT))*NFFT*sqrt(2*domegasim);
%         end
%     end
%     XX{1,z+1}=x;
% %     if strcmpi(FileName,'No'); XX{1,z+1}=x; else save([FileName '_Nr_' num2str(z) '.mat'],'t','x'); end
% end
% 
% if nargout==1
%     varargout{1}=XX;
% elseif nargout==2
%     varargout{1}=XX{1};
%     varargout{2}=XX(2:end);
% end

%% Simulate time series
% omegaaxissim=domegasim:domegasim:max(omegaaxisinput);
% NFFT=2^nextpow2(2*length(omegaaxissim));
% t=linspace(0,2*pi/domegasim,NFFT);
% XX=cell(1,Nsim+1);
% if strcmpi(FileName,'No');  XX{1,1}=t; end
% for z=1:Nsim
%     phi=2*pi*rand(size(SS,1),length(omegaaxissim));
%     x=zeros(size(SS,1),NFFT);
%     for m=1:size(GG,1)
%         for n=1:m
%             c=interp1(omegaaxisinput,permute((GG(m,n,:)),[1,3,2]),omegaaxissim,'linear','extrap').*exp(1i.*phi(n,:));
%             x(m,:)=x(m,:)+real(ifft(c,NFFT))*NFFT*sqrt(2*domegasim);
%         end
%     end
%     if strcmpi(FileName,'No'); XX{1,z+1}=x; else save([FileName '_Nr_' num2str(z) '.mat'],'t','x'); end
% end

