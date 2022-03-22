
clc
clear all

%initialState = [-65, 0.05, 0.6, 0.32]; % V, m, h, n
SimulationTime = 10;
deltaT=.001;
timeInterval = 0:deltaT:SimulationTime;
currentLevel = 20e-12;


Tau1 = 3*50;
Tau2 = 3*450;
Tau3 = 3*800;
b = 1.5*3800;
%C_m = 0.02;

Idark(1:numel(timeInterval)) = -40e-12;

A(1:2000) =  0;
%A(20001:6000) = currentLevel;
A(2001:numel(timeInterval)) = currentLevel;
for i=1:numel(timeInterval)
    iphoto(i) = Idark(i) + A(i)*((1 - exp(-i/Tau1)) - 1.9*(1/(1 + (exp(-(i - b)/Tau2)))) + (1 - exp(-i/Tau3)));
    %Iphoto(i) = Idark(i) + A(i)*(1 - exp(-i/Tau3))*(1-(1/(1 + (exp(-(i - b)/Tau2)))));
    %Vm(i+1) = Vm(i) + deltaT*(-Iphoto(i)/Cm);
    
    %ikv = gbar_kv*(mkv.^3)*hkv*(Vm(i) - Ekv);
    
end


% Constants
    Cm = 0.02e-9;
    alpha1 = 50;
    alpha2 = 0.0003;
    alpha3 = 0.03;
    eps = 0.5e-6;
    Ttot = 1000e-6;
    beta1 = 2.5;
    Tau1 = 0.2e-6;
    Tau2 = 5;
    PDEtot = 100;
    b = 0.25e-18;
    Jmax = 5040e-12;
    eT = 500.0e-6;
    GammaCa = 50;
    Ca0 = 0.1e-6;
    k1 = 0.2e-6;
    k2 = 0.8;
    Amax = 65.6e-6;
    Kc = 100e-6;
    sigma = 1e-6;
    Vbar = 0.4; 
    
    gbar_kv = 2.0e-9;
    Ekv = -80;
    
    gbar_h = 0.305e-9;
    Eh = -32;
    
    gbar_ca = 1.2e-9;
    CA0 = 1600e-6;
    
    gbar_cacl = 6.5e-9;
    Ecacl = -45;
    
    gbar_kca = 0.5e-9;
    Ek = -80;
    
    gl = 0.5e-9;
    El = -55;
    
    F = 9.648e4;
    V1 = 3.812e-13;
    V2 = 5.236e-13;
    delta = 5.9e-5;
    S1 = 3.142e-8;
    Lb1 = 0.4e-6;
    Lb2 = 0.2e-6;
    Hb1 = 100e-6;
    Hb2 = 90e-6;
    Bl = 500e-6;
    BH = 200e-6;
    Jex = 20e-12;
    Jex2 = 20e-12;
    Dca = 6e-8;
    Cae = 0.01e-6;
    
 currentLevels = 10000;
%{
Jhv(1:200000) =  0;
Jhv(200001:202000) = currentLevels;
Jhv(202001:204000) = 0;
Jhv(204001:206000) = currentLevels;
Jhv(208001:208000) = 0;
Jhv(208001:210000) = currentLevels;
Jhv(210001:212000) = 0;
Jhv(212001:214000) = currentLevels;
Jhv(214001:216000) = 0;
Jhv(216001:218000) = currentLevels;
Jhv(218001:220000) = 0;
Jhv(220001:222000) = currentLevels;
Jhv(222001:224000) = 0;
Jhv(224001:226000) = currentLevels;
Jhv(226001:228000) = 0;
Jhv(228001:230000) = currentLevels;
Jhv(230001:232000) = 0;
Jhv(232001:234000) = currentLevels;
Jhv(234001:236000) = 0;
Jhv(236001:238000) = currentLevels;
Jhv(238001:numel(timeInterval)) = 0;
%}

%{
Jhv(1:200) =  0;
Jhv(201:220) = currentLevels;
Jhv(221:240) = 0;
Jhv(241:260) = currentLevels;
Jhv(261:280) = 0;
Jhv(281:310) = currentLevels;
Jhv(311:numel(timeInterval)) = 0;
%}
%{
Jhv(2121:2140) = currentLevels;
Jhv(2141:2160) = 0;
Jhv(2161:2180) = currentLevels;
Jhv(2181:2200) = 0;
Jhv(2201:2220) = currentLevels;
Jhv(2221:2240) = 0;
Jhv(2241:2260) = currentLevels;
Jhv(2261:2280) = 0;
Jhv(2281:2300) = currentLevels;
Jhv(2301:2320) = 0;
Jhv(2321:2340) = currentLevels;
Jhv(2341:2360) = 0;
Jhv(2361:2380) = currentLevels;
Jhv(2381:numel(timeInterval)) = 0;
%}

% Initial Values
V = -36.18;
Rh = 0;
Rhi = 0;
Tr = 0;
PDE = 0;
Ca2 = 0.3e-6;
Cab = 34.88e-6;
cGMP = 2.0;

% Ikv constant
mkv = 0.430;
hkv = 0.999;

%Ih constant
c1 = 0.646;
c2 = 0.298;
o1 = 0.0517;
o2 = 0.00398;
o3 = 0.000115;

%Ica constant
CaS = 0.0966e-6;
mca = 0.436;

mkca = 0.642;
CaF = 0.0966e-6;
Cabls = 80.929e-6;
Cabhs = 29.068e-6;
Cablf = 80.929e-6;
Cabhf = 29.068e-6;


% Iphoto
j = (Jmax*(cGMP.^3))/((cGMP.^3) + 1000);
%iphoto = -j*(1 - exp((V - 8.5)/17.0));

% Ikv
almkv = (5*(100 - V))/(exp((100-V)/42)-1);
bmkv = 9*(exp(-(V-20)/40));
alhkv = 0.15*(exp(-V/22));
bhkv = 0.4125/(exp((10-V)/7) + 1);
ikv = gbar_kv*(mkv.^3)*hkv*(V - Ekv)*1e-5;

% Ih
alh = 3/(exp((V + 88)/14)+1);
bh = 18/(exp(-(V + 18)/19)+1);
ih = gbar_h*(o1+o2+o3)*(V - Eh)*1e-2;

%Ica
almca = 3*(80-V)/(exp((80-V)/25) - 1);
bmca = 10/(1+(exp(V+38/7)));
hca = exp((40-V)/18)/(1+exp((40-V)/18));
Eca = -12.9*log(CaS/CA0);
ica = gbar_ca*(mca.^4)*hca*(V-Eca)*1e-5;

%Icacl
mcacl = 1/(1 + (exp((0.37 - CaS)/0.09)));
iclca = gbar_cacl*mcacl*(V - Ecacl)*1e-3;

%Ikca
almkca = 15*(80-V)/(exp((80-V)/40) - 1);
bmkca = 20*exp(-V/35);
mkcas = CaS/(CaS + 0.3);
ikca = gbar_kca*(mkca.^2)*mkcas*(V - Ek)*1e4;

%Il
il = gl*(V - El)*1e-3;

% Calcium System
iex = Jex*exp(-(V+14)/70)*((CaS - Cae)/(CaS - Cae + 2.3));
iex2 = Jex2*((CaS - Cae)/(CaS - Cae + 0.5));
jev = 0;


%Synapse
Vth = -51;
Vslope = 10;
gsyn_max = 2.56e-9;
tau = 0.1;
s_inf = tanh((V - Vth)/Vslope);
S = 0.09;



for i=1:numel(timeInterval)-1 %Compute coefficients, currents, and derivates at each time step
   
    %---calculate the coefficients---%
    %Equations here are same as above, just calculating at each time step
    j(i) = (Jmax*power(cGMP(i), 3))/(power(cGMP(i), 3) + power(10, 3));
    %iphoto(i) = -j(i) *(1 - exp((V(i) - 8.5)/17.0));
    jev(i) = 1 - exp((V(i) - 8.5)/17.0);
    
    almkv(i) = (5*(100 - V(i)))/(exp((100-V(i))/42)-1);
    bmkv(i) = 9*(exp(-(V(i)-20)/40));
    alhkv(i) = 0.15*(exp(-V(i) /22));
    bhkv(i) = 0.4125/(exp((10-V(i))/7) + 1);
    ikv(i) = gbar_kv*power(mkv(i), 3)*hkv(i) *(V(i) - Ekv)*1e-5;
    
    alh(i) = 8/(exp((V(i) + 78)/14)+1);
    bh(i) = 20/(exp(-(V(i) + 8)/19) + 1);
    ih(i) = gbar_h*(o1(i) + o2(i) + o3(i))*(V(i) - Eh)*1e-2;
    
    almca(i) = 3*(80-V(i))/(exp((80-V(i))/25) - 1);
    bmca(i) = 10/(1+(exp((V(i)+38)/7)));
    hca(i) = exp((40-V(i))/18)/(1+exp((40-V(i))/18));
    Eca(i) = -12.5*log(CaS(i)/CA0);
    ica(i) = gbar_ca*power(mca(i) , 4)*hca(i) *(V(i) - Eca(i))*1e-5;
    
    mcacl(i) = 1/(1 + (exp((0.37 - CaS(i))/0.09)));
    iclca(i) = gbar_cacl*mcacl(i)*(V(i) - Ecacl)*1e-3;
    
    almkca(i) = 15*(80-V(i))/(exp((80-V(i))/40) - 1);
    bmkca(i) = 20*exp(-V(i)/35);
    mkcas(i) = CaS(i)/(CaS(i) + 0.3);
    ikca(i) = gbar_kca*power(mkca(i), 2)*mkcas(i)*(V(i) - Ek)*1e4;
    
    
    il(i) = gl*(V(i) - El)*1e-3;
    
    iex(i) = Jex*exp(-(V(i)+14)/70)*((CaS(i) - Cae)/(CaS(i) - Cae + 2.3));
    iex2(i) = Jex2*((CaS(i) - Cae)/(CaS(i) - Cae + 0.5));

    s_inf(i) = tanh((V(i) - Vth)/Vslope);

    %---calculate the currents---%
    I_ALL = -(iphoto(i) + ih(i) + ikv(i) + il(i) + ikca(i) + iclca(i) + iex(i) + iex2(i) + ica(i));
    %
    
   
    %---calculate the derivatives using Euler first order approximation---%
    V(i+1) = V(i) + (deltaT*I_ALL*(1e3)/Cm);
    
    
    % Photo Current
    Rh(i+1) = Rh(i) + deltaT*(100 - alpha1*Rh(i) + alpha2*Rhi(i));
    Rhi(i+1) = Rhi(i) + deltaT*(alpha1*Rh(i) - (alpha2+alpha3)*Rhi(i));
    Tr(i+1) = Tr(i) + deltaT*(eps*Rh(i)*(Ttot - Tr(i)) - beta1*Tr(i) + Tau2*PDE(i) - Tau1*Tr(i)*(PDEtot - PDE(i)));
    PDE(i+1) = PDE(i) + deltaT*(Tau1*Tr(i)*(PDEtot - PDE(i)) - Tau2*PDE(i));
    Ca2(i+1) = Ca2(i) + deltaT*(b*j(i) - GammaCa*(Ca2(i) - Ca0) - k1*(eT - Cab(i))*Ca2(i) + (k2*Cab(i)));
    Cab(i+1) = Cab(i) + deltaT*(k1*(eT - Cab(i))*Ca2(i) - (k2*Cab(i)));
    cGMP(i+1) = cGMP(i) + deltaT*((Amax/(1 + power((Ca2(i)/Kc), 4))) - (cGMP(i)*(Vbar + sigma*PDE(i))));
    

    % Hyperpolarize Current
    c1(i+1) = c1(i) + deltaT*(-4*alh(i)*c1(i) + bh(i)*c2(i));
    c2(i+1) = c2(i) + deltaT*(4*alh(i)*c1(i) - 3*(alh(i) + bh(i))*c2(i) + 2*bh(i)*o1(i));
    o1(i+1) = o1(i) + deltaT*(((3*alh(i)*c2(i)) - ((2*alh(i) + 2*bh(i))*o1(i)) + (3*bh(i)*o2(i))));
    o2(i+1) = o2(i) + deltaT*(2*alh(i)*o1(i) - (alh(i) + 3*bh(i))*o2(i) + 4*bh(i)*o3(i));
    o3(i+1) = o3(i) + deltaT*( alh(i)*o2(i) - 4*bh(i)*o3(i));
    
    % Delayed Rectifier Current
    mkv(i+1) = mkv(i) + deltaT*(almkv(i)*(1-mkv(i)) - bmkv(i)*mkv(i));
    hkv(i+1) = hkv(i) + deltaT*(alhkv(i)*(1-hkv(i)) - bhkv(i)*hkv(i));
    
    % Calcium Current
    mca(i+1) = mca(i) + deltaT*(almca(i)*(1 - mca(i)) - bmca(i)*mca(i));
    
    % Calcium Activated Potassium Current
    mkca(i+1) = mkca(i) + deltaT*(almkca(i)*(1 - mkca(i)) - bmkca(i)*mkca(i));
    
    % Intracellular Calcium System
    CaS(i+1) = CaS(i) + deltaT*(-(((ica(i) + iex(i) + iex2(i))/(2*F*V1))*1e-6) - (Dca*(S1/(delta*V1))*(CaS(i)-CaF(i))) - (Lb1*CaS(i)*(Bl - Cabls(i))) + Lb2*Cabls(i) - (Hb1*CaS(i)*(BH - Cabhs(i))) + Hb2*Cabhs(i));
    CaF(i+1) = CaF(i) + deltaT*((Dca*(S1/(delta*V2))*(CaS(i)-CaF(i))) - (Lb1*CaF(i)*(Bl - Cablf(i))) + Lb2*Cablf(i) - (Hb1*CaF(i)*(BH - Cabhf(i))) + Hb2*Cabhf(i));
    Cabls(i+1) = Cabls(i) + deltaT*(Lb1*CaS(i)*(Bl - Cabls(i)) - Lb2*Cabls(i));
    Cabhs(i+1) = Cabhs(i) + deltaT*(Hb1*CaS(i)*(BH - Cabhs(i)) - Hb2*Cabhs(i));
    Cablf(i+1) = Cablf(i) + deltaT*(Lb1*CaF(i)*(Bl - Cablf(i)) - Lb2*Cablf(i));
    Cabhf(i+1) = Cabhf(i) + deltaT*(Hb1*CaF(i)*(BH - Cabhf(i)) - Hb2*Cabhf(i));

    % Chemical Synapse
    S(i+1) = S(i) + deltaT*((s_inf(i) - S(i)))/(((1 - s_inf(i))*tau*S(i)));

end


%V = V-65; %Set resting potential to -70mv

%il(numel(il)+1) = il(numel(il));
plot(timeInterval, A, 'r')
%xlim([0 10])
legend('Step Response of Light Stimulus vs time')

%{
figure
plot(timeInterval, S, 'r', 'linewidth', 1)
%xlim([0 0.1])
legend('S(t) vs time')

figure
s_inf(numel(s_inf)+1) = s_inf(numel(s_inf));
plot(timeInterval, s_inf, 'r', 'linewidth', 1)
%xlim([0 0.1])
legend('S_inf vs time')


%}

%figure
%iphoto(numel(iphoto)+1) = iphoto(numel(iphoto));
plot(timeInterval, iphoto, 'r')
%xlim([0 10])
legend('Iphoto vs time')
hold on

figure
%timeInterval
plot(timeInterval, V, 'b', 'linewidth', 1)
%xlim([0 10])
legend('Vm vs time')
hold on



figure
ih(numel(ih)+1) = ih(numel(ih));
plot(timeInterval, ih, 'g')
%xlim([0 10])
legend('Ih vs time')
%plot(timeInterval, gbar_kv*(n.^3), 'r')
%legend('GKV vs time')
hold on


%figure
%plot(timeInterval, Iinj, 'r')
%legend('Jhv vs time')

figure
ikv(numel(ikv)+1) = ikv(numel(ikv));
plot(timeInterval, ikv, 'r')
%xlim([0 10])
legend('Ikv vs time')
hold on

%{
figure
%Jhv(numel(Jhv)+1) = Jhv(numel(Jhv));
plot(timeInterval, Jhv, 'r')
%xlim([0 10])
legend('Jhv vs time')
%}

figure
ica(numel(ica)+1) = ica(numel(ica));
plot(timeInterval, ica, 'r')
%xlim([0 10])
legend('Ica vs time')
hold on


%{
figure
plot(timeInterval, Rh, 'r')
%xlim([0 10])
legend('Rh vs time')

figure
plot(timeInterval, Rhi, 'r')
%xlim([0 10])
legend('Rhi vs time')

figure
plot(timeInterval, Tr, 'r')
%xlim([0 10])
legend('Tr vs time')

figure
plot(timeInterval, PDE, 'r')
%xlim([0 10])
legend('PDE vs time')

figure
plot(timeInterval, cGMP, 'r')
%xlim([0 10])
legend('cGMP vs time')


figure
plot(timeInterval, Ca2, 'r')
%xlim([0 10])
legend('Ca vs time')


figure
plot(timeInterval, Cab, 'r')
%xlim([0 10])
legend('Cab vs time')
%}

figure
ikca(numel(ikca)+1) = ikca(numel(ikca));
plot(timeInterval, ikca, 'r')
%xlim([0 10])
legend('Ikca vs time')
hold on


figure
iclca(numel(iclca)+1) = iclca(numel(iclca));
plot(timeInterval, iclca, 'r')
%xlim([0 10])
legend('Iclca vs time')
hold on


figure
il(numel(il)+1) = il(numel(il));
plot(timeInterval, il, 'r')
%xlim([0 10])
legend('Il vs time')
hold on

figure
almkv(numel(almkv)+1) = almkv(numel(almkv));
plot(timeInterval, almkv, 'r')
%xlim([0 10])
legend('mkv vs time')
%}    
   