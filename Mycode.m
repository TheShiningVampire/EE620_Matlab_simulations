%% Constants
q=1.6e-19;
k=8.314/(6.023*10^23);
T=300;
eps_0=8.85e-12; 
Vt=k*T/q; %thermal voltage

%% Substrate parameters

N_sub=-1e17*1e6;    % -ve=p type substrate, +ve=n type substrate, input doping density here
ni=1.5e16;          % number of instrinsic carriers in silicon substrate
Eg_sub=1.12*q;           %band gap of silicon, 1.1eV
eps_sub=12*eps_0; %permitivitty of silicon substrate
EA_sub=4.05*q;       %electron affinity for substrate
%L_sub=1e-6;
%calculate n and p for substrate
if(N_sub<0)%p type substrate(NMOS)
    p=abs(N_sub);
    n=ni^2/p;
elseif(N_sub>0)%n type substrate(PMOS)
    n=N_sub;
    p=ni^2/n;
else %intrinsic semiconductor
    n=ni;
    p=pi;
end   

%% oxide parameters
eps_ox=3.9*eps_0;  %permitivitty of SiO2
t_ox=10e-9;        %input the oxide thickness required here
Ec_offset=3.160*q; %conduction band offset
Eg_ox=9*q;
c_ox=eps_ox/t_ox;

%% Metal Parameters;
phi_m = (4.05)*q;  %work function of metal
phi_b = -sign(N_sub)*Vt*q*log(abs(N_sub)/ni);%Ei-Ef for semiconductor
phi_sub=EA_sub+Eg_sub/2+phi_b;  %work function for substrate
V_fb=(phi_m-phi_sub)/q ;  %flatband voltage
V_g=1;             %gate voltage
%% Solver
dx=1e-9;
NR_max=10;
tol1=1e-6;

%% Input
psi_s=0.7941;
pre=(2*eps_sub*q*Vt*abs(N_sub))^0.5;

for i = 1:10
psi0=psi_s;   
Q_s = pre*(exp(-psi_s/Vt)+psi_s/Vt-1+ni^2/(abs(N_sub))^2*(exp(psi_s/Vt)-psi_s/Vt-1))^0.5;
f= -V_g+V_fb-Q_s/c_ox+psi_s;
df=-1/c_ox*(pre)^2*0.5*(-1/Vt*exp(-psi_s/Vt)+1/Vt+ni^2/abs(N_sub)^2*(1/Vt*exp(psi_s/Vt)-1/Vt))/Q_s + 1;
psi_s=psi_s-f/df;
if(abs(psi_s-psi0)<1e-6)
   break;
end

end
 