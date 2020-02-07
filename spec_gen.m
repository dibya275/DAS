clear;clc;

T=1700;
M=18.01;
x=0.0186;
P=1;
L=1;

%%
format long;
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

d=d0.*(296/T).^0.00;
g_a=g_a0.*(296/T).^n;
g_s=g_s0.*(296/T).^0.75;
nu0=nu+P.*(1-x).*d;
S0=s0.*2.479371939e19*P*x*L;
S=linestrength(S0,nu0,Edd,T,partfun);
delnu_g=0.5*7.162242257e-7.*nu0.*sqrt(T/M);
delnu_l=P.*(x.*g_s+(1-x).*g_a);

par0=[nu0 S delnu_g delnu_l]';


%%
dat=linspace(7184,7187,10000)';
fit=voigt(dat,par0);
dat1=load('1700-3-0.0186.txt');

%%
plot(dat,fit,'g',dat1(:,1),dat1(:,2),'r');


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
S=S0.*(T0/T.*(174.5813/Q(T,partfun)).*exp(-h*c.*Edd/k*(1/T-1/T0)).*(1-exp(-h*c.*nu0/(k*T)))./(1-exp(-h*c.*nu0/(k*T0))));
end



