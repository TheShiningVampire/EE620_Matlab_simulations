caps = C;

cox = caps(1);
cmin = caps(end);

t_ox_extraction = eps_ox/cox;

cdep = 1/( 1/cmin - 1/cox );
wdep = eps_si / cdep ; 


ni = 1.5e10*1e6;
% Initial guess
Na_extraction = 1e17*1e6;

for i = 1:15
    phi_b_extraction = Vt*log(Na_extraction/ni);
    Na_extraction = 4*eps_si*phi_b_extraction/q/wdep^2;
end

phi_b_extraction = Vt*log(Na_extraction/ni);
Eg = q*1.1;
Vfb_extraction = -Eg/(2*q) - phi_b_extraction;

phi_sub_extraction=(Ea+Eg/2+phi_b_extraction*q);
phi_m_extraction = q*V_fb + phi_sub_extraction;
phi_m_extraction = phi_m_extraction/q;

