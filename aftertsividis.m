%% Input parameters
phi_m =4.05 ;      %work function of metal
V_g=1.5;            %gate voltage
t_ox=10e-9;          %oxide thickness 
N_sub=1e17*1e6;      %-ve=p type substrate(NMOS), +ve=n type substrate(PMOS), input doping density here


%% Constants
q=1.6e-19;
k=8.314/(6.023*10^23);
T=300;
eps_0=8.85e-12; 
Vt=k*T/q;          %thermal voltage

%% Substrate parameters

ni=1.5e16;          % number of instrinsic carriers in silicon substrate
Eg_sub=1.12*q;      %band gap of silicon, 1.12eV
eps_sub=12*eps_0;   %permitivitty of silicon substrate
EA_sub=4.05*q;      %electron affinity for substrate
L_sub=1e-6;
k_sub=12;

%calculating the value of n and p for substrate
if(N_sub<0)                %p type substrate(NMOS)
    p_sub=abs(N_sub);
    n_sub=ni^2/p_sub;
elseif(N_sub>0)            %n type substrate(PMOS)
    n_sub=N_sub;
    p_sub=ni^2/n_sub;
else                       %intrinsic semiconductor
    n_sub=ni;
    p_sub=pi;
end   

%% Oxide parameters

eps_ox=4*eps_0;    %permitivitty of SiO2
Ec_offset=3.1*q;   %conduction band offset
Eg_ox=9*q;         %band gap of oxide
c_ox=eps_ox/t_ox;
k_ox=4;

%% Workfunction,flatband values;

phi_b = -sign(N_sub)*Vt*log(abs(N_sub)/ni);   %Ei-Ef/q for semiconductor
phi_s = EA_sub/q+Eg_sub/(2*q)+phi_b;          %work function for substrate
V_fb=(phi_m-phi_s) ;                          %flatband voltage
tm = 10e-9;                                   %Gate thickness, used for polysilicon, not of use in metalgate

%% Solver parameters
dx=1e-9;
I_max=50;
tol1=1e-6;

% Sign of substrate charges in the current mode, since we have signum(psi_s) in formula of charge

if(V_g>=V_fb)
    signq = -1;       %Depletion for NMOS, Accumulation for PMOS
else                  % i.e. (V_g<V_fb)
    signq = +1;       %Accumulation for NMOS, Depletion for PMOS
end

%Same charge sign for depletion and inversion

if sign(N_sub) <= 0
    F = @(s) V_fb-(signq*((2*eps_sub*k*T*abs(N_sub))^0.5*((exp(-q*s/(k*T))+q*s/(k*T)-1)+(ni^2/(abs(N_sub))^2)*(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g;
    psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve
    
else                              % sign(N_sub) > 0
    F = @(s) V_fb-(signq*((2*eps_sub*k*T*abs(N_sub))^0.5*((ni^2/(abs(N_sub))^2)*(exp(-q*s/(k*T))+q*s/(k*T)-1)+(exp(+q*s/(k*T))-q*s/(k*T)-1))^0.5)/c_ox)+s-V_g;
    psi_s = fsolve(F,-signq*0.1); % Initial guess, sign taken due to fsolve
end

%% Newton raphson for calculating potential inside substrate,as explained in tutorial

x = 0:dx:L_sub;                 % semiconductor mesh, dx seperation, from 0 to L_sub
N = length(x);                  % number of points= length of array x

A = (-2*diag(ones(1, N),0) + 1*diag(ones(1,N-1),1) + 1*diag(ones(1,N-1),-1)); %constructing A matrix
A(1, :) = 0; 
A(1,1) = 1;
A(N, :) = 0; 
A(N,N) = 1;

b = zeros(N,1); 
b(1)= psi_s;
b(N) = 0;

psi = zeros(N,1);  % guess for psi for newton raphson

%% Calculating psi matrix using newton raphson method, reference taken from class tutorial

for i=1:I_max
    p = p_sub*exp(-psi/Vt);
    n = n_sub*exp(psi/Vt);
    rho=(p-n+N_sub)*q;
    b= -(rho/eps_sub)*dx^2;
    b(1) = psi_s;
    b(N) = 0;        
    f2 = A*psi - b;
    
 
    drho= (-p/Vt-n/Vt)*q;
    db=-drho/eps_sub*dx^2;
    db(1)=0; db(N)=0;           %since b(1) and b(end) are constants as per our boundary conditions
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

%% Band Diagram

%silicon substrate

Ei = -q*psi;             % Ei in bulk is taken as zero reference
Ec = Ei + Eg_sub/2;
Ev = Ei - Eg_sub/2;
Ef = (Ei(end) - q*phi_b) * ones(size(Ei));

%oxide

Ecox = -q*psiox + Eg_sub/2 + Ec_offset;
Evox = Ecox - Eg_ox;

V_g = V_fb + Vox + psi(1);
%metal               
xm = -t_ox - tm:dx:-t_ox;             
Em = (Ef(end) - q*V_g)*ones(size(xm)); 

figure(1);
plot([xm xox x]*1e9, [Em'; Ecox'; Ec]/q, 'r', [xm xox x]*1e9, [Em';Evox'; Ev]/q, 'b', x*1e9, Ei/q,'k--',x*1e9,Ef/q,'k');
xlim([-20, 200]);
xlabel('Position (nm)'); 
ylabel('Energy (eV)');
legend('E_C', 'E_V', 'E_i', 'E_F');


% carrier concentration with distance
% Charge Density
figure(2)
semilogy(x*1e9, p*1e-6, 'r', x*1e9, n*1e-6, 'b');
xlabel('Position (nm)');
ylabel('Concentration (cm^{-3})');
xlim([0, 200]);
legend('p (Holes)', 'n (Electrons)');
 






