% Solving Brews model for the following NMOSFET
% L = 1 um, Tox = 10 nm, n+ poly-Si gate, Na (choose) so that V_T -> 0.8V
% range for the Vg and Vd sweeps is 0-5V

% Plot the Id - Vg characteristics for different values of Vd (transfer)
% Plot the Id - Vd characteristics for different values of Vg ( Vg > Vt )

% Compare it with brews and piece-wise model 

% Re-do above for a PMOSFET, with p+ poly Si gate
% use const eff. mobility values of = 200 cm^2 / Vs (e-) and 100 (holes)

% Constants 
q = 1.6e-19;
eps_0 = 8.85e-12;
kT = 26e-3*q;

% Semiconductor
Nsub = -4e17*1e6;       % -ve p-substrate -> nmos
k_si = 12;
ni = 1.5e10*1e6;
Eg = 1.1*q;
eps_si = k_si*eps_0;
chi_si = 4.05*q;
Na = abs(Nsub);
L = 1e-6;
W = 1e-6;
if Nsub < 0
    muf = 200e-4;
else 
    muf = 100e-4;
end

%Oxide
tox = 10e-9;
k_ox = 4;
eps_ox = k_ox*eps_0;
Cox = eps_ox/tox;

% Metal 
if Nsub < 0
    phi_m = chi_si/q;
else 
    phi_m = chi_si/q + Eg/q;
end

phi_b = -sign(Nsub) * kT/q * log(abs(Nsub)/ni);
phi_s = chi_si/q + Eg/(2*q) + phi_b;
Vfb = phi_m - phi_s;

% %% Brews Solution at Id = f_B(Vg, Vd) 
% Vg = 1;
% Vd = 2;
% 
% QT = @(psi_s) Cox*(Vg - Vfb - psi_s);
% QD = @(psi_s) (2*eps_si*q*Na*psi_s).^0.5;
% QI = @(psi_s) QT(psi_s) - QD(psi_s);
% dVdpsi_s = @(psi_s) 1 + 2*kT/q * (Cox*QT(psi_s) + eps_si*q*Na)./ ((QT(psi_s)).^2 - (QD(psi_s)).^2);
% Vgf = @(psi_s, V) Vfb + psi_s + 1/Cox*(2*eps_si*kT*Na)^0.5*(q*psi_s/kT + ni^2/Na^2*exp(q*(psi_s-V)/kT)).^0.5;
% 
% psi_s_min = -Vg-Vd-abs(Vfb);
% psi_s_max = Vg+Vd+abs(Vfb);
% dpsi_s = 1e-3;
% psi_svec = psi_s_min:dpsi_s:psi_s_max;
% 
% psi_ss = interp1(real(Vgf(psi_svec,0)), real(psi_svec),Vg);
% psi_sd = interp1(real(Vgf(psi_svec,Vd)), real(psi_svec), Vg);
% 
% intf = @(psi_s) QI(psi_s).*dVdpsi_s(psi_s);
% 
% Id = muf*W/L*integral(intf, psi_ss, psi_sd);

%% Plot the Id - Vg characteristics for different values of Vd (transfer)
% Vd = 0.1 V (Lin)
% Vd = 1 V (Sat)
Vd = 0.1;
Id_vec = [];
Vg_vec = 0:50e-3:5;
dpsi_s = 1e-3;

for Vg = 0:50e-3:5
    QT = @(psi_s) Cox*(Vg - Vfb - psi_s);
    QD = @(psi_s) (2*eps_si*q*Na*psi_s).^0.5;
    QI = @(psi_s) QT(psi_s) - QD(psi_s);
    dVdpsi_s = @(psi_s) 1 + 2*kT/q * (Cox*QT(psi_s) + eps_si*q*Na)./ ((QT(psi_s)).^2 - (QD(psi_s)).^2);
    Vgf = @(psi_s, V) Vfb + psi_s + 1/Cox*(2*eps_si*kT*Na)^0.5*(q*psi_s/kT + ni^2/Na^2*exp(q*(psi_s-V)/kT)).^0.5;

    psi_s_min = -Vg-Vd-abs(Vfb);
    psi_s_max = Vg+Vd+abs(Vfb);  
    psi_svec = psi_s_min:dpsi_s:psi_s_max;

    psi_ss = interp1(real(Vgf(psi_svec,0)), real(psi_svec),Vg);
    psi_sd = interp1(real(Vgf(psi_svec,Vd)), real(psi_svec), Vg);

    intf = @(psi_s) QI(psi_s).*dVdpsi_s(psi_s);
    Id = muf*W/L*integral(intf, psi_ss, psi_sd);
    
    Id_vec = [Id_vec Id];
end

figure(2);
yyaxis left;
semilogy(Vg_vec, Id_vec*1e6);
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);

yyaxis right; 
plot(Vg_vec, Id_vec*1e6);
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);

hold on;

Vd = 1;
Id_vec = [];
Vg_vec = 0:50e-3:5;
dpsi_s = 1e-3;

%functions

for Vg = 0:50e-3:5
    QT = @(psi_s) Cox*(Vg - Vfb - psi_s);
    QD = @(psi_s) (2*eps_si*q*Na*psi_s).^0.5;
    QI = @(psi_s) QT(psi_s) - QD(psi_s);
    dVdpsi_s = @(psi_s) 1 + 2*kT/q * (Cox*QT(psi_s) + eps_si*q*Na)./ ((QT(psi_s)).^2 - (QD(psi_s)).^2);
    Vgf = @(psi_s, V) Vfb + psi_s + 1/Cox*(2*eps_si*kT*Na)^0.5*(q*psi_s/kT + ni^2/Na^2*exp(q*(psi_s-V)/kT)).^0.5;

    psi_s_min = -Vg-Vd-abs(Vfb);
    psi_s_max = Vg+Vd+abs(Vfb);  
    psi_svec = psi_s_min:dpsi_s:psi_s_max;

    psi_ss = interp1(real(Vgf(psi_svec,0)), real(psi_svec),Vg);
    psi_sd = interp1(real(Vgf(psi_svec,Vd)), real(psi_svec), Vg);

    intf = @(psi_s) QI(psi_s).*dVdpsi_s(psi_s);
    Id = muf*W/L*integral(intf, psi_ss, psi_sd);
    
    Id_vec = [Id_vec Id];
end

figure(2);
yyaxis left;
semilogy(Vg_vec, Id_vec*1e6,'--');
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);

yyaxis right; 
plot(Vg_vec, Id_vec*1e6,'--');
axes = gca;
axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
xlabel('Vg'); ylabel('Id (\mu A/\mu m)');
lines = axes.Children;
set(lines, 'LineWidth', 2);

hold off; 

% %% Plot Id - Vd characteristics for different values of Vd (transfer)
% Vd_vec = 0:50e-3:5;
% dpsi_s = 1e-3;
% 
% for Vg = 0.5:0.5:2
%     Id_vec = [];
%     for Vd = 0:50e-3:5
%         QT = @(psi_s) Cox*(Vg - Vfb - psi_s);
%         QD = @(psi_s) (2*eps_si*q*Na*psi_s).^0.5;
%         QI = @(psi_s) QT(psi_s) - QD(psi_s);
%         dVdpsi_s = @(psi_s) 1 + 2*kT/q * (Cox*QT(psi_s) + eps_si*q*Na)./ ((QT(psi_s)).^2 - (QD(psi_s)).^2);
%         Vgf = @(psi_s, V) Vfb + psi_s + 1/Cox*(2*eps_si*kT*Na)^0.5*(q*psi_s/kT + ni^2/Na^2*exp(q*(psi_s-V)/kT)).^0.5;
% 
%         psi_s_min = -Vg-Vd-abs(Vfb);
%         psi_s_max = Vg+Vd+abs(Vfb);  
%         psi_svec = psi_s_min:dpsi_s:psi_s_max;
% 
%         psi_ss = interp1(real(Vgf(psi_svec,0)), real(psi_svec),Vg);
%         psi_sd = interp1(real(Vgf(psi_svec,Vd)), real(psi_svec), Vg);
% 
%         intf = @(psi_s) QI(psi_s).*dVdpsi_s(psi_s);
%         Id = muf*W/L*integral(intf, psi_ss, psi_sd);
% 
%         Id_vec = [Id_vec Id];
%     end
%     
%     figure(3);
%     plot(Vd_vec, Id_vec*1e6);
%     axes = gca;
%     axes.LineWidth = 1; axes.FontSize = 14; axes.FontWeight = 'bold'; axes.Box = 'on';
%     xlabel('Vd'); ylabel('Id (\mu A/\mu m)');
%     lines = axes.Children;
%     set(lines, 'LineWidth', 2);
%     hold on;
% end
% figure(3);
% lgd = legend('0.5','1','1.5','2');
% title(lgd,'V_G')
% hold off;



    




















