%% generate voigt profile  using HITRAN data
clear all;
tic
%par0=load('par0.txt');
dat=load('spectraplot.dat');
dat1=load('hitran_ru.dat');
vrange = (min(dat(:,1))-0.5):0.001:(max(dat(:,1))+0.5);
hitran_par
%%
%[parmin,resnom,res,exitflag]= fit2voigt(dat,par0);
%%
%fit=voigt(dat(:,1),parmin);
fit=voigt(vrange',par0);
%%
%subplot(2,1,1);
hold on
plot(vrange,fit,'g',dat(:,1),dat(:,2),'b',dat1(:,1),dat1(:,2)*10,'k');
toc
% clear all
% legend('data','fit')
% xlim([600,1700]);ylim([-0.01, 0.25]);
% subplot(2,1,2);
% plot(dat(:,1),res,'k'); 
% legend('residual','location','north');
% xlim([600,1700]);ylim([-7, 7]*1e-3);