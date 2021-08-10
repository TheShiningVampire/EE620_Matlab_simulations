
%% Parameters

m0 = 9.1e-31;
h = 6.6e-34;
q = 1.6e-19;
kB = 1.38e-23;
hbar = h/(2*pi);
eps_0 = 8.85e-12;
k_si = 12; 
eps_si = k_si*eps_0;
k_ox = 4;
eps_ox = k_ox*eps_0;
chi_si = 4.05;
Eg_ox = 9*q;

% Inputs
m_si_per = 0.98*m0;
m_si_par = 0.19*m0;
eta = 2;
tox = 2.5e-9;
Eox = 10*1e6/1e-2;
T = 300;
phit = 0;


mox = 0.38*m0;
pref = 4*pi*m_si_per*q/h^3;
pref1 = -2/hbar*(2*mox)^0.5;
kBT = kB*T;


%% On my way

Foxs = 0:0.001:8;
Foxs = Foxs;

J1s = [];

J2s = [];

T2s = [];
T1s = [];

for Fox = Foxs
    disp(Fox);
    Q = eta*q*m_si_par/(pi*hbar^2);
    
    f = 0.6*( 2*q/( 3*pi*hbar*q*m_si_per )^(1/3) ) * ( eps_ox*Fox/eps_si )^(2/3);
    
    E_si_per = 0.6*(( 3*pi*hbar*q*m_si_per )^(2/3) /(2*m_si_per) ) * ( eps_ox*Fox/eps_si )^(2/3);
    v_si_per = @(E_si_per) (2*E_si_per / m_si_per)^0.5;
    
    qphi_cat = @(E_si_par) q*chi_si -(E_si_per + E_si_par) ; 
    qphi_an = @(E_si_par) q*chi_si -(E_si_per + E_si_par)-q*Fox*tox ;
    
    gamma = @(Eox) Eox*(1- Eox/Eg_ox);
    gamma_prime = @(Eox) 1 - 2*Eox /Eg_ox;
    
     
     fac = @(Eox) 2*gamma_prime(Eox) * gamma(Eox)^0.5 + (Eg_ox)^.5 * asin(gamma_prime(Eox));
     Twkb =@(E_si_par) exp( Eg_ox*(2*mox)^0.5/(4*hbar*q*Fox) * (fac(qphi_cat(E_si_par)) - fac(qphi_an(E_si_par)) )  );
     
     vox = @(Eox) (2*gamma(Eox)/mox )^0.5/(gamma_prime(Eox));
     
     TR = @(E_si_par) 4*v_si_per(E_si_per)*vox(qphi_cat(E_si_par)) * 4*v_si_per(E_si_per + q*Fox*tox)*vox(qphi_an(E_si_par))/ ( (   (v_si_per(E_si_per))^2 + (vox(qphi_cat(E_si_par)))^2   ) * (  (v_si_per(E_si_per + q*Fox*tox ))^2  + (vox(qphi_an(E_si_par)))^2  ) )   ;                       
     
     
     T1 = Twkb(0.5*pi*hbar^2*eps_ox*Fox/(q*m_si_par * eta) );
     
     T1s = [T1s T1];
     
     T2 = Twkb(0.5*pi*hbar^2*eps_ox*Fox/(q*m_si_par * eta) ) * TR(E_si_per);
     
     T2s = [T2s T2];
     
     J1 = Q*f*T1;
     J1s = [J1s J1];
    
     
     
     J2 = Q*f*T2;
     J2s = [J2s J1];
    
end

hold on
plot(Foxs , log10(J1s));
plot(Foxs , log10(J2s));

