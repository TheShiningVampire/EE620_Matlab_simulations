%% input parameters
phi_m =(4.05) ;      %work function of metal
V_g=-3;            %gate voltage
t_ox=10e-9;          %oxide thickness 
N_sub=-1e17*1e6;     % -ve=p type substrate(NMOS), +ve=n type substrate(PMOS), input doping density here

%% Constants
q=1.6e-19;
k=8.314/(6.023*10^23);
T=300;
eps_0=8.85e-12; 
Vt=k*T/q;         %thermal voltage

%% Substrate parameters
ni=1.5e16;          % number of instrinsic carriers in silicon substrate
Eg_sub=1.12*q;      %band gap of silicon, 1.12eV
eps_sub=12*eps_0;   %permitivitty of silicon substrate
EA_sub=4.05*q;      %electron affinity for substrate
L_sub=1e-6;
k_sub=12;
%calculate n and p for substrate
if(N_sub<0)%p type substrate(NMOS)
    p_sub=abs(N_sub);
    n_sub=ni^2/p_sub;
elseif(N_sub>0)%n type substrate(PMOS)
    n_sub=N_sub;
    p_sub=ni^2/n_sub;
else %intrinsic semiconductor
    n_sub=ni;
    p_sub=pi;
end   

%% oxide parameters
eps_ox=4*eps_0;  %permitivitty of SiO2
Ec_offset=3.1*q; %conduction band offset
Eg_ox=9*q;         %band gap of oxide
c_ox=eps_ox/t_ox;
k_ox=4;

%% Wf,flatband values;

phi_b = -sign(N_sub)*Vt*log(abs(N_sub)/ni); %Ei-Ef/q for semiconductor
phi_s=EA_sub/q+Eg_sub/(2*q)+phi_b;%work function for substrate
V_fb=(phi_m-phi_s) ;  %flatband voltage
tm = 10e-9; 
%% Solver
dx=1e-9;
I_max=50;
tol1=1e-6;
tol2=1e-7;
% %% Calculating psi_s from Vg
psi_s=-sign(V_g-V_fb)*0.01;             %guess for psi_s
pre=(2*eps_sub*q*abs(N_sub))^0.5;

func = @(psi_s) -V_g + V_fb -(1/c_ox -sign(V_g-V_fb)*pre(exp(-psi_s/Vt)+psi_s/Vt-1+(ni^2/(abs(N_sub))^2)*(exp(psi_s/Vt)-psi_s/Vt-1))^0.5) + psi_s;
x0 = psi_s;
psi_s = fzero(func,x0);

%% I TRIED TO DEBUG MY NEWTON RAPHSON BUT WAS NOT ABLE TO FIND THE PROBLEM
%% So I used built in function fzero to complete the assignment. The following block is my newton raphson

% for i = 1:50
% psi0=psi_s;   
% Q_s = sign(psi_s)*pre*(Vt*exp(-psi_s/Vt)+psi_s-Vt+(ni^2/(abs(N_sub))^2)*(Vt*exp(psi_s/Vt)-psi_s-Vt))^0.5;
% f= -V_g+V_fb+Q_s/c_ox+psi_s;
% df=1/c_ox*(pre)^2*0.5*(-1*exp(-psi_s/Vt)+1+(ni^2/abs(N_sub)^2)*(1*exp(psi_s/Vt)-1))/Q_s + 1;
% psi_s=psi0-f/df;
% if(abs(psi_s-psi0)<tol2)
%    break;
% end
% 
% end
 
%psi_s = phi_b;                  % Surface Potential

%% Mesh
x = 0:dx:L_sub;                 % semiconductor mesh, dx seperation, from 0 to L_sub
N = length(x);                  % number of points= length of array x

%% Newton raphson for calculating potential inside substrate
A = (-2*diag(ones(1, N),0) + 1*diag(ones(1,N-1),1) + 1*diag(ones(1,N-1),-1)); %constructing A matrix
A(1, :) = 0; 
A(1,1) = 1;
A(N, :) = 0; 
A(N,N) = 1;

b = zeros(N,1); %Constructing B matrix
b(1)= psi_s;
b(N) = 0;

psi = zeros(N,1);  % guess for psi

%% Calculating psi matrix using newton raphson method, reference taken from class tutorial

for i=1:I_max
    p = p_sub*exp(-psi/Vt);
    n = n_sub*exp(psi/Vt);
    rho=(p-n+N_sub)*q;
    b= -rho/eps_sub*dx^2;
    b(1) = psi_s;
    b(N) = 0;        
    f2 = A*psi - b;
    
 
    drho= (-1/Vt*p+1/Vt*n)*q;
    db=-drho/eps_sub*dx^2;
    db(1)=0; db(N)=0;           %since b(1) and b(end) are constants
    Jacobian= A - diag(db);
    dV=-Jacobian\f2;
    
   
    psi = psi + dV;
    
end

%% Addition of Oxide

E_s = (psi(1) - psi(2))/dx;          % E = -dpsi/dx, which can be approximated by finite difference
E_ox = eps_sub/eps_ox*E_s;           % Epsilon*field same across interface (D=epsilon*E)    
xox = -t_ox:dx:0;                    
psiox = psi(1) - E_ox*xox;           
Vox = E_ox*t_ox;                     % voltage drop for the oxide 

% %% Plots
% figure(1); plot([xox x]*1e9, [psiox'; psi]);
% xlabel('Position (nm)'); ylabel('Potential (V)');
% xlim([-10, 150]);
% figure(2);
% plot([xox(1:end-1) x(1:end-1)]*1e9, [-diff(psiox)'/dx;-diff(psi)/dx]*1e-8);
% xlabel('Position (nm)'); ylabel('Electric Field (MV/m)');
% xlim([-10, 150]);

%% Band Diagram

Ei = -q*psi;                        % Ei in bulk is taken as zero reference
Ec = Ei + Eg_sub/2;
Ev = Ei - Eg_sub/2;
Ef = (Ei(end) - q*phi_b) * ones(size(Ei)); % Because phi_b = (Ei - Ef)/q

Ecox = -q*psiox + Eg_sub/2 + Ec_offset;
Evox = Ecox - Eg_ox;

V_g = V_fb + Vox + psi(1);
Em = Ef(end) - q*V_g;                % Because q*Vg = Em - Ef
xm = -t_ox - tm:dx:-t_ox;             % Metal mesh
Em = Em*ones(size(xm)); 
figure(4);
plot([xm xox x]*1e9, [Em'; Ecox'; Ec]/q, 'r', [xm xox x]*1e9, [Em';Evox'; Ev]/q, 'b', x*1e9, Ei,'k--',x*1e9,Ef/q,'k');
xlim([-20, 150]);
axes = gca;
axes.LineWidth=1;axes.FontSize=18;axes.FontWeight='bold';axes.Box='on';
xlabel('Position (nm)'); ylabel('Energy (eV)');
lines = axes.Children;
set(lines, 'LineWidth', 1.5);
leg = legend('E_C', 'E_V', 'E_i', 'E_F');
leg.Orientation = 'horizontal';