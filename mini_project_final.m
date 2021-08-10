
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
m_si_per = 0.32*m0;
m_si_par = 0.25*m0;
eta = 4;
tox = 2.2e-9;
Eox = 10*1e6/1e-2;
T = 300;
phit = 0;


mox = 0.38*m0;
pref = 4*pi*m_si_per*q/h^3;
pref1 = -2/hbar*(2*mox)^0.5;
kBT = kB*T;


x = 0:1e-10:1e-9;

%% On my way

Foxs = 0:0.001:7;
Foxs = Foxs * 1e8;

Js_with_reflection = [];
Js_without_reflection = [];


for F_ox = Foxs
    disp(F_ox);
    % Using the relations mentioned in the refernce paper
    Q = F_ox * eps_ox; 
    
    f = 0.6 * ((2*q)/((3 * pi * hbar * q * m_si_per)^(1/3)))*...
                     ((eps_ox/eps_si) * F_ox)^(2/3);
    
    E_si_per = 0.6 * (((3 * pi * hbar * q * m_si_per)^(2/3))/(2 * m_si_per))...
                   * ((eps_ox/eps_si) * F_ox)^(2/3);
    Ef = E_si_per + (pi * hbar^2 * eps_ox * F_ox)/(q * m_si_par * eta);
    E_si_par = 0.5 * (Ef - E_si_per);
    
    
    v_si_per_cat = sqrt((2 * E_si_per)/m_si_per);
    v_si_per_an = sqrt((2 * E_si_per + q * F_ox * tox)/m_si_per);
    
    Eox_prime = - q * (F_ox);
    Eox = Eox_prime * x;
    
    gamma = Eox .* ( 1 - (Eox./Eg_ox));
    gamma_prime = 1 - 2*(Eox./Eg_ox);
    
    vox = (1./gamma_prime).*sqrt((2*gamma)/(mox)); 
    
    q_phi_cat = Ec_off - (E_si_per + E_si_par); 
    q_phi_an = Ec_off - (E_si_per + E_si_par) - q*F_ox*tox;
    
    gamma_an = q_phi_an * (1 - q_phi_an/Eg_ox);
    gamma_an_prime = 1 - 2*(q_phi_an/Eg_ox);
    gamma_cat = q_phi_cat * (1 - q_phi_cat/Eg_ox);    
    gamma_cat_prime = 1 - 2*(q_phi_cat/Eg_ox);
    
    v_ox_cat = (1/gamma_cat_prime)*sqrt((2*gamma_cat)/mox);
    v_ox_an = (1/gamma_an_prime)*sqrt((2*gamma_an)/mox);


    T_wkb = exp(((Eg_ox*sqrt(2 * mox))/(4 * hbar * q * F_ox))*...
                (2*gamma_cat_prime*sqrt(gamma_cat) - 2*gamma_an_prime*sqrt(gamma_an)+...
                sqrt(Eg_ox)*(asin(gamma_cat_prime) - asin(gamma_an_prime))));
    
    
    TR = ((4 * v_si_per_cat * v_ox_cat)/...
            (v_si_per_cat^2 + v_ox_cat^2))*...
            ((4 * v_si_per_an * v_ox_an)/...
            (v_si_per_an^2 + v_ox_an^2));
        
    J_without_reflection = Q*f*T_wkb;
    J_with_reflection = Q*f*T_wkb* TR;
    
    Js_without_reflection = [Js_without_reflection J_without_reflection];
    Js_with_reflection = [Js_with_reflection J_with_reflection];
end
figure(1);
%semilogy(Foxs * 1e-2 , Js_without_reflection * 1e-4, 'Linewidth' , 2);
semilogy(Foxs * 1e-2, Js_with_reflection * 1e-4, 'Linewidth' , 2);

grid on;
xlabel('Oxide Field F_{ox} (MV-cm^{-1})');
ylabel('Gate Current J_g (A-cm^{-2})');
title('Gate Current vs Oxide Field');

