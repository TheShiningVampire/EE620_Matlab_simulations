
%% Parameters

m0 = 9.1e-31;
h = 6.6e-34;
q = 1.6e-19;
kB = 1.38e-23;
hbar = h/(2*pi);

% Inputs
meff = 0.3*m0;
tox = 5e-9;
Eox = 10*1e6/1e-2;
T = 300;
phit = 3;


mox = meff;
pref = 4*pi*meff*q/h^3;
pref1 = -2/hbar*(2*mox)^0.5;
kBT = kB*T;

%% 

close all;
V = Eox * tox;

Ef1 = q*V;
Ef2 = 0;

f1 = @(E) 1./(1 + exp((E-Ef1)/kBT));
f2 = @(E) 1./(1 + exp((E-Ef2)/kBT));

Ec = @(x) Ef1 + q*phit - q*Eox*x;

dJs = [];
Ns = [];
TCs = [];

Emin = Ef2 - 0.5*q;
Emax = Ef1 + q*phit;
Einf = q*10;
dEx = q*0.1;
Exs = Emin:dEx:Emax;

for Ex = Exs
    E = @(Ep) (Ex + Ep);
    fp = @(Ep) (f1(E(Ep)) - f2(E(Ep)));
    N = integral(fp , 0 , Einf);
    Ns = [Ns N];
    
    x1 = (Ef1 + q*phit - Ex)/(q*Eox);
    TCx = @(x) ( Ec(x) - Ex ).^0.5;
    
    TC = exp(pref1 * integral(TCx , 0, x1 ) );
    TCs = [TCs TC];
    
    dJ = pref*N*TC;
    dJs = [dJs dJ];
    
end

J = sum(dJs)*dEx;

plot(Exs, dJs );



