clear;clc;

T=1700;
M=18.01;
x=0.0186;
P=1;
L=10;

%%
format long;
dat=load('1700-3-0.0186.txt');
lines=importhitran('1392-3.par');
nu=lines.transitionWavenumber;
s0=lines.lineIntensity;
g_a0=lines.airBroadenedWidth;
g_s0=lines.selfBroadenedWidth;
Edd=lines.lowerStateEnergy;
n=lines.temperatureDependence;
d0=lines.pressureShift;
fil1=fopen('q1.txt','r');   
partfun=fscanf(fil1,'%f %f', [2 Inf]); 
partfun=partfun';

d=d0.*(296/T).^0.0;
g_a=g_a0.*(296/T).^n;
g_s=g_s0.*(296/T).^0.75;
nu0=nu+P.*(1-x).*d;
S0=s0.*2.479371939e19*P*x;         % effect of pressure, length, mole fraction comes here
S=linestrength(S0,nu0,Edd,T,partfun);
delnu_g=0.5*7.162242257e-7.*nu0.*sqrt(T/M);
delnu_l=P.*(x.*g_s+(1-x).*g_a);

par0=[nu0 S delnu_g delnu_l]';

%%
[parmin,resnom,res,exitflag]= fit2voigt(dat,par0);
%%
a=voigt(dat(:,1),parmin);
S_max=max(parmin(2,:));
j=find(parmin(2,:)==S_max);
par_plot=parmin(:,j);
a_j=voigt(dat(:,1),par_plot);
phi_j=a_j/S_max;
intphi_j=trapz(dat(:,1),phi_j);
A_j=trapz(dat(:,1),a_j);
S_j=A_j/intphi_j;
T0=1500;
solT=fzero(@(T)(S_j-linestrength(S0(j),nu0(j),Edd(j),T,partfun)),T0);

%%
plot(dat(:,1),a*L);
xlabel('Wavenumber (cm^{-1})')
ylabel('Absorbance (entire region)')
figure()
plot(dat(:,1),a_j*L);
xlabel('Wavenumber (cm^{-1})')
ylabel('Absorbance (seperated line)')
figure()
subplot(2,1,1);
plot(dat(:,1),dat(:,2)*L,'b',dat(:,1),a*L,'r');
xlabel('Wavenumber (cm^{-1})')
ylabel('Absorbance (entire region)')
legend('data','fit')
subplot(2,1,2);
plot(dat(:,1),res,'k'); 
xlabel('Wavenumber (cm^{-1})')
ylabel('Residual')



function q=Q(T,partfun)
T=round(T);
T1=T-1;
T2=T+1;
Tind1=find(partfun(:,1)==T1,1);
Tind2=find(partfun(:,1)==T2,1);
q=interp1(partfun(Tind1:Tind2,1),partfun(Tind1:Tind2,2),T);
end

function S=linestrength(S0,nu0,Edd,T,partfun)
T0=296;
h=6.62607004e-34;
c=29979245800;
k=1.38064852e-23;
S=S0.*(T0/T.*(174.5813./Q(T,partfun)).*exp(-h*c.*Edd/k*(1/T-1/T0)).*(1-exp(-h*c.*nu0/(k*T)))./(1-exp(-h*c.*nu0/(k*T0))));
end