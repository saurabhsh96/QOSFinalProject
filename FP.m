%Final presentation project Quasi Optical
clear;
%% Defining inputs

%Operation frequency
freq = 150e9;

%Antenna Spec
Da = 30e-2;
Fa = 40e-2;
Ra = Da/2;

%Incidence %TM
drad = pi/180;
th_i = 5*drad;
ph_i = 0*drad;

%Lens specification
er = 11.9;
majAxis = 15e-3; %Major Axis
a = majAxis/2; %semi major axis

%EM Spec
c = 3e8;
lam = c/freq;
k0 = 2*pi/lam;
kd = sqrt(er).*k0;
n = 2; %Exponent in Field's equation

%Meshgrid definitions
%Mesh for Lens field calculations
drad = pi/180;
dthl = drad;
dphl = drad;
[thl, phl] = meshgrid(eps:dthl:pi/2, eps:dphl:2*pi);

%Temperory observation point
r_temp = 1000*lam;

%Observation Mesh
dtho = drad/2;
dpho = pi/180;
[theta_obs, phi_obs] = meshgrid(eps:dtho:pi/2-dtho, eps:dpho:2*pi);

%% Calculation of lense dimensions from given dimensions

e = 1/sqrt(er);
%Taking critical angle as the maximum angle ? Question ?
thC = pi/2 - asin(e);
Th0 = (thC);
b = tan(Th0).*e.*a; %minor axis of the lens
cent = a.*e; %center of ellipse
h = a.*(1+e); %height of the lens
rmin = a.*(1-e^2)./(1-e.*cos(Th0));
Dl = 2*rmin.*sin(Th0); %Diameter of the apparent aperture

jf = [1; 0; 0];
%% Part 0: Radiation pattern of the lens

[EthInL, EphInL, PradFeed, rho, phi_s, JsX, JsY, JsZ, ...
    Axx, Ayy, EFLx, EFLy, EFLz, Prad, etaRad] = LensPatterns(freq, er, n, Dl, Th0, jf, ...
    dtho, dpho, theta_obs, phi_obs);

% field inside the lens
% [Eth, Ephi] = FieldMedia(n, kd, thl, phl, r_temp);
% Zero component
% zeroComp = zeros(size(thl));
% 
% [Emagl, Emaxl] = absVal({Eth, Ephi, zeroComp});
% Plotting
% figure(1);
% name = "Electric field";
% plotReq([-thl(1,size(thl, 2):-1:1)./drad thl(1,:)./drad],...
%     [mag2db(Emagl(1,size(thl, 2):-1:1)./Emaxl), mag2db(Emagl(1,:)./Emaxl)] ...
%     ,"\theta","Normalized Electric Field","Electric field inside lens, \phi = 0", ...
%     [-60, 0], name);
% legend show;

%% Reflector dimentions calculation

fNum = Fa/Da;
Th0ant = 2*atan(Da/(4*Fa)); %Max Angle
theta = linspace(eps, Th0ant, 150);
phi_req = linspace(eps, 2*pi, 400);
[th, ph] = meshgrid(theta, phi_req);
dth = th(1,2) - th(1,1);
dph = ph(2,1) - ph(1,1);

E0 = 1; %Assume

%Calculating Vgo
[Vgoth0, Vgoph0, Egoth0, Egoph0] = GOField(E0, th, ph, Fa, k0);
EgoMag = sqrt(abs(Egoth0).^2 + abs(Egoph0).^2);
EgoMax = max(max(EgoMag));

%Calculating Vth, now, we already have field by the feed, trauncate it to
%required Th0ant
%% Calculate Ethfeed for above Th0ant

[~, ~, PradFeed1, ~, ~, ~, ~, ~, ...
    ~, ~, EFLx1, EFLy1, EFLz1, Prad1, etaRad1] = LensPatterns(freq, ...
    er, n, Dl, Th0, jf, dth, dph, th, ph);

%% Calculating vth from field recieved

%r_obs to remove the constant term
r_obs = 1000*lam;

%Constant term, phase and amp
remove_term = exp(-1j*k0*r_obs)/r_obs;
constant_term = exp(-1j*k0*Fa)/Fa;

EthLReq = EFLx1.*cos(th).*cos(ph).*(constant_term/remove_term) ...
    + EFLy1.*cos(th).*sin(ph).*(constant_term/remove_term) ...
    - EFLz1.*sin(th).*(constant_term/remove_term);
EphLReq = - EFLx1.*sin(ph).*(constant_term/remove_term) ...
    + EFLy1.*cos(ph).*(constant_term/remove_term);

VthL = EthLReq./constant_term;
VphL = EphLReq./constant_term;

%% Power calculation for normal incidence

%[PrxN, etaAN] = PowerRec(Vath0, Vaph0, Vgoth0, Vgoph0, Prad2, zeta0, th, ph, Rr, E0);
delTh = (1 - cos(th))./(1+cos(th));
%In th, ph space, basically the reflector space
kRhox = k0.*sin(th).*cos(ph);
kRhoy = k0.*sin(th).*sin(ph);

%Defining DelK
delKx = k0.*sin(th_i).*cos(ph_i);
delKy = k0.*sin(th_i).*sin(ph_i);

%Define DelRho
delRhox = (Fa/k0).*delKx;
delRhoy = (Fa/k0).*delKy;

%Dot products
kDel = delRhox.*kRhox + delRhoy.*kRhoy;

%Exponent Term
expTerm1 = exp(-1j.*kDel);
expTerm2 = exp(-1j.*kDel.*delTh);

%GO-field and voltage at given incident angle
EgothDel = Egoth0.*expTerm1.*expTerm2;
EgophDel = Egoph0.*expTerm1.*expTerm2;
VgothDel = Vgoth0.*expTerm1.*expTerm2;
VgophDel = Vgoph0.*expTerm1.*expTerm2;

%Change dx accordingly
dx = -40*lam:0.0001:40*lam;
PrxD = zeros(size(dx));
etaAD = zeros(size(dx));

kx = k0.*sin(th).*cos(ph);
zeta0 = 377;

for ind = 1:size(dx, 2)
    VthDx =  VthL.*exp(1j*kx*dx(ind));
    VphDx =  VphL.*exp(1j*kx*dx(ind));
    [PrxD(ind), etaAD(ind)] = PowerRec(VthDx, VphDx, VgothDel, VgophDel ...
        , Prad, zeta0, th, ph, Ra, E0);
end

%%
figure;
%plot(dx, pow2db(PrxD))
name = "Power Recieved";
plotReq(dx.*10^3,...
    pow2db(PrxD./(max(PrxD))) ...
    ,"\Delta d (in mm)","Normalized Recieve Power (in dBm)","Power Recieved by Quasi Optical System", ...
    [-60, 0], name);
%legend show;
%% Power received vs. theta incident

%Change the angle of incidence Th, keeping the phi cons
%Defining meshgrid for incident angles
thi_arr = linspace(0*drad, 12*drad, 220);
phi_arr = linspace(eps, 2*pi, 360);
[th_in, ph_in] = meshgrid(thi_arr, phi_arr);
dthi = th_in(1,2) - th_in(1,1);
dphi = ph_in(2,1) - ph_in(1,1);

loc = dx(PrxD == max(PrxD));

VthLoc =  VthL.*exp(1j*kx*loc);
VphLoc =  VphL.*exp(1j*kx*loc);

%Recieved power calculations
PrxRDx = zeros(size(th_in));
etaARDx = zeros(size(th_in));

for indPh = 1:size(th_in, 1)
    for ind = 1:size(th_in, 2)
        %Defining DelK
        delKx_ = k0.*sin(th_in(indPh, ind)).*cos(ph_in(indPh, ind));
        delKy_ = k0.*sin(th_in(indPh, ind)).*sin(ph_in(indPh, ind));

        %Define DelRho
        delRhox_ = (Fa/k0).*delKx_;
        delRhoy_ = (Fa/k0).*delKy_;

        %Dot products
        kDel_ = delRhox_.*kRhox + delRhoy_.*kRhoy;

        %Exponent Term
        expTerm1_ = exp(-1j.*kDel_);
        expTerm2_ = exp(-1j.*kDel_.*delTh);

        %GO-field and voltage at given incident angle
        VgothDel_ = Vgoth0.*expTerm1_.*expTerm2_;
        VgophDel_ = Vgoph0.*expTerm1_.*expTerm2_;

        [PrxRDx(indPh, ind), etaARDx(indPh, ind)] = PowerRec(VthLoc, ...
            VphLoc, VgothDel_, VgophDel_, Prad, zeta0, th, ph, Ra, E0);
    end
end

figure;
%plot(dx, pow2db(PrxD))
name = "Power Recieved";
plotReq(th_in(1,:)./drad,...
    pow2db(PrxRDx(1,:)./max(max(PrxRDx))) ...
    ,"\theta_{incident} deg","Normalized Recieve Power (in dBm)" ...
    ,"Recieved power vs. Angle of incidence \phi = 0 deg", ...
    [-50, 0], name);

figure;
%plot(dx, pow2db(PrxD))
name = "Power Recieved";
plotReq(th_in(1,:)./drad,...
    pow2db(PrxRDx(91,:)./max(max(PrxRDx))) ...
    ,"\theta_{incident} deg","Normalized Recieve Power (in dBm)" ...
    ,"Recieved power vs. Angle of incidence \phi = 90 deg", ...
    [-50, 0], name);

%% Directivity and efficiencies of the given FO system

%Rx Approach
%Directivity
denom = sum(PrxRDx.*sin(th_in), 'all').*dthi.*dphi;
Dir = 4.*pi.*PrxRDx./(denom);

Dmax = max(max(Dir));

Gain = max(max(etaARDx)).*Dir;
Gmax = max(max(Gain));

%Spillover efficiency; PradTotal by Prad at Th0ant %check once again
%etaS = Prad1/Prad;
%For full half space, VthL
EthLReqR = EFLx.*cos(theta_obs).*cos(phi_obs).*(constant_term/remove_term) ...
    + EFLy.*cos(theta_obs).*sin(phi_obs).*(constant_term/remove_term) ...
    - EFLz.*sin(theta_obs).*(constant_term/remove_term);
EphLReqR = - EFLx.*sin(phi_obs).*(constant_term/remove_term) ...
    + EFLy.*cos(phi_obs).*(constant_term/remove_term);

VthLR = EthLReqR./constant_term;
VphLR = EphLReqR./constant_term;

VaMagR = (abs(VthLR)).^2 + (abs(VphLR)).^2;
%theta0
VaMag = (abs(VthL)).^2 + (abs(VphL)).^2;

etaS = ((sum(sum(VaMag.*sin(th)))).*dth.*dph)./ ...
    (sum(sum(VaMagR.*sin(theta_obs))).*dtho.*dpho);

taperS = ((abs(sum((VthL.*Vgoth0 + VphL.*Vgoph0).*sin(th), 'all').*dth.*dph)).^2)./ ...
    ((pi*(Da/2).^2).*(sum(sum(VaMag.*sin(th))).*dth.*dph));

apS = etaS.*taperS;
apS1 = max(max(etaARDx));

%Tx Approach
dRho = 0.0001;
rho_dash1R = eps:dRho:Ra;
theta0R = 2*atan(Ra*2/(4*Fa));
[rho_dashR, phi_dashR] = meshgrid(rho_dash1R, eps:drad:2*pi);
th_dashR = 2*atan(rho_dashR./(2*Fa));
r_dashR = Fa.*(1+(tan(th_dashR/2)).^2);

[tapSTx, ~, Area] = Spillover(freq, 1, jf, Dl/2, 1000*lam, ...
            th_dashR, phi_dashR, thl, phl, rho_dashR, Fa, Ra);

etaStx = Prad1/Prad;
etaATx = tapSTx.*etaStx;

%Maximum directiviy
DirTxM = (4*pi/lam^2).*Area;

%Directivity
DirTx = DirTxM.*tapSTx;

%Gain
GTx = DirTxM.*etaATx;

figure;
%plot(dx, pow2db(PrxD))
name = "Directivity";
plotReq(th_in(1,:)./drad,...
    pow2db(Dir(1,:)) ...
    ,"\theta_{incident} deg","Directivity (in dBi)" ...
    ,"Directivity, gain vs. Angle of incidence \phi = 0 deg", ...
    [0, 50], name); hold on;

name = "Gain";
plotReq(th_in(1,:)./drad,...
    pow2db(Gain(1,:)) ...
    ,"\theta_{incident} deg","Gain (in dBi)" ...
    ,"Directivity, gain vs. Angle of incidence \phi = 0 deg", ...
    [0, 50], name);
legend show;
grid on;
%% Changing the dimensions of the lens and reflector

