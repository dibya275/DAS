%% calculate Voigt parameters from HITRAN data
p = 1; %pr atm
chi = 0.00141;%mole fr of H2O
p_self = p*chi;
L = 10;%cm
NA = 6.022e23;%Avogadro number
c = 2.99792458e10;% light speed cm/s
K = 1.3806488e-16; %erg/K
M = 18.01; %gm/mol
T = 296;%K
load H2O_multiline.par

v0 = H2O_multiline(:,2)+H2O_multiline(:,4)*p; %cm-1
s = H2O_multiline(:,3);  %cm-1/(molec.cm-2) at 296K
s = s*1013250/(K*T);
s = s*p*chi*L;
s =s*2*sqrt(log(2)/pi);
al = H2O_multiline(:,5)*p_self+(p-p_self)*H2O_multiline(:,6); %Lorentzian lnewidth par
ag = H2O_multiline(:,2)/c*sqrt(2*NA*K*T*log(2)/M); %gaussian linewidth Ref: hitran.ru[HWHM]

par0 = [v0';s';ag';al']*1.000;

