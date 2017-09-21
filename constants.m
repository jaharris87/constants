%% Fundamental, astrophysical, and numerical constants

%% numerical constants

% integers
zero      = 0.0;
one       = 1.0;
two       = 2.0;
three     = 3.0;
four      = 4.0;
five      = 5.0;
six       = 6.0;
seven     = 7.0;
eight     = 8.0;
nine      = 9.0;
ten       = 10.0;
eleven    = 11.0;
twelve    = 12.0;
fifteen   = 15.0;
sixteen   = 16.0;
sixty     = 60.0;

% fractions
half      = 0.5;
third     = one / three;
fourth    = 0.25;
fifth     = 0.2;
sixth     = one / six;
seventh   = one / seven;
eighth    = 0.125;
ninth     = one / nine;
tenth     = 0.1;
twelfth   = one / twelve;
fifteenth = one / fifteen;
sixteenth = 0.0625;
sixtieth  = one / sixty;
two3rd    = two * third;
four3rd   = four * third;
five3rd   = five * third;
three4th  = 0.75;
two5th    = 0.4;
three5th  = 0.6;
four5th   = 0.8;
five6th   = five * sixth;
three8th  = 0.375;
five8th   = 0.625;
seven8th  = 0.875;
five12th  = five * twelfth;
seven12th = seven * twelfth;
two15th   = two * fifteenth;
eight15th = eight * fifteenth;

% 71 digits of pi (256-bit precision)
pi        = 3.1415926535897932384626433832795028841971693993751058209749445923078164;
pi2       = pi * pi;
pi3       = pi * pi2;
pi4       = pi * pi3;
pi5       = pi * pi4;
twpi      = pi * two;
thpi      = pi * three;
frpi      = pi * four;
etpi      = pi * eight;
sxtnpi    = pi * sixteen;
frpith    = frpi * third;
thpifrth  = thpi * fourth;
pi_inv    = one / pi;
twpi_inv  = pi_inv * half;
frpi_inv  = pi_inv * fourth;
etpi_inv  = pi_inv * eighth;

% logarithms
log_e     = log10( exp(one) );
ln_2      = log( two );
ln_10     = log( ten );


%% fundamental physical constants (CODATA recommended 2014 values)
%   P. J. Mohr, D. B. Newell, and B. N. Taylor, Rev. Mod. Phys. 88, 035009 (2016)
%   http://physics.nist.gov/cuu/Constants/Table/allascii.txt
% * whenever possible, derived constants are calculated explicitly for floating-point consistency

% common unit conversion factors
eV2MeV = 1.0e-6;
MeV2eV = 1.0e+6;

eV2GeV = 1.0e-9;
GeV2eV = 1.0e+9;

MeV2GeV = 1.0e-3;
GeV2MeV = 1.0e+3;

erg2J = 1.0e-7;
J2erg = 1.0e+7;

cm2fm = 1.0e+13;
fm2cm = 1.0e-13;

cm32fm3 = cm2fm * cm2fm * cm2fm;
fm32cm3 = fm2cm * fm2cm * fm2cm;

% speed of light in vacuum
c_light = 2.99792458e+10; % [cm s^{-1}]
c_light_inv = one / c_light;
c_light2 = c_light * c_light;
c_light2_inv = one / c_light2;
c_light3 = c_light * c_light2;
c_light3_inv = one / c_light3;

% common mass-energy conversion factors
g2erg = c_light2; % = 8.987551787368176e+20;
erg2g = c_light2_inv; % = 1.112650056053618e-21;

% elementary charge
q_e_C = 1.6021766208e-19; % [C]
q_e = q_e_C * tenth * c_light; % = 4.80320467e-10; % [esu] or [erg^{1/2} cm^{1/2}]
q_e2 = q_e * q_e;

% common electronvolt conversion factors
eV2J = q_e_C; % = 1.6021766208e-19;
J2eV = one / eV2J; % = 6.241509126e+18;

eV2erg = eV2J * J2erg; % = 1.6021766208e-12;
erg2eV = erg2J * J2eV; % = 6.241509126e+11;

eV2g = eV2erg * erg2g; % = 1.782661907e-33;
g2eV = g2erg * erg2eV; % = 5.609588650e+32;

% Planck constant
h = 6.626070040e-27; % [erg s]
h_MeV = h * erg2eV * eV2MeV; % = 4.135667662e-21; % [MeV s]
h_inv = one / h;
h2 = h * h;
h2_inv = one / h2;
h3 = h * h2;
h3_inv = one / h3;

% common inverse length conversion factors
cm_inv2g = h * c_light_inv; % = 2.210219057e-37;
g2cm_inv = c_light * h_inv; % = 4.524438411e+38;

cm_inv2erg = h * c_light; % = 1.986445824e-16;
erg2cm_inv = h_inv * c_light_inv; % = 5.034116651e+15;

cm_inv2eV = cm_inv2erg * erg2eV; % = 1.2398419739e-4;
eV2cm_inv = eV2erg * erg2cm_inv; % = 8.065544005e+3;

% Planck constant, reduced; Dirac constant
hbar = twpi_inv * h; % = 1.054571800e-27; % [erg s] 
hbar_MeV = hbar * erg2eV * eV2MeV; % = 6.582119514e-22; % [MeV s]
hbar_inv = one / hbar;
hbar2 = hbar * hbar;
hbar2_inv = one / hbar2;
hbar3 = hbar * hbar2;
hbar3_inv = one / hbar3;

% common conversion factor
hbarc = hbar * c_light; % [erg cm] 
hbarc_MeV = hbar_MeV * c_light * cm2fm; % = 197.3269788; % [MeV fm] 
hbarc_inv = one / hbarc;

% Newtonian constant of gravitation
G_N = 6.67408e-8; % [cm^3 g^-1 s^-2]
G_N_inv = one / G_N;

% Boltzmann constant
k_B = 1.38064852e-16; % [erg K^{-1}]
k_B_MeV = k_B * erg2eV * eV2MeV; % = 8.6173303e-11; % [MeV K^{-1}]
k_B_inv = one / k_B;
k_B4 = k_B * k_B * k_B * k_B;
k_B4_MeV = k_B_MeV * k_B_MeV * k_B_MeV * k_B_MeV;

% common temperature-energy conversion factors
K2erg = k_B; % = 1.38064852e-16;
erg2K = one / K2erg; % = 7.2429731e+15;

K2eV = k_B * erg2eV; % = 8.6173303e-5;
eV2K = one / K2eV; % = 1.16045221e+4;

% Avogadro constant
N_A = 6.022140857e+23; % [mol^{-1}]

% atomic mass unit
m_u = one / N_A; % = 1.660539040e-24; % [g]
m_u_eV = m_u * g2eV; % = 931.4940954e+6; % [eV]
m_u_MeV = m_u_eV * eV2MeV; % = 931.4940954; % [MeV]

% common atomic-mass-unit conversion factors
u2g = m_u; % = 1.660539040e-24;
g2u = N_A; % = 6.022140857e+23;

u2eV = u2g * g2eV; % = 931.4940954e+6 
eV2u = eV2g * g2u; % = 1.0735441105e-9 

% electron mass
m_e = 9.10938356e-28; % [g]
m_e_erg = m_e * g2erg; % = 8.18710565e-7; % [erg]
m_e_eV = m_e * g2eV; % = 0.5109989461e+6; % [eV]
m_e_MeV = m_e_eV * eV2MeV; % = 0.5109989461; % [MeV]
m_e_u = m_e * g2u; % = 5.48579909070e-4; % [u]

% proton mass
m_p = 1.672621898e-24; % [g]
m_p_erg = m_p * g2erg; % = 1.503277593e-3; % [erg]
m_p_eV = m_p * g2eV; % = 938.2720813e+6; % [eV]
m_p_MeV = m_p_eV * eV2MeV; % = 938.2720813; % [MeV]
m_p_u = m_p * g2u; % = 1.007276466879; % [u]

% neutron mass
m_n = 1.674927471e-24; % [g]
m_n_erg = m_n * g2erg; % = 1.505349739e-3; % [erg]
m_n_eV = m_n * g2eV; % = 939.5654133e+6; % [eV]
m_n_MeV = m_n_eV * eV2MeV; % = 939.5654133; % [MeV]
m_n_u = m_n * g2u; % = 1.00866491588; % [u]

% fine-structure constant
fine = q_e2 * hbarc_inv; % = 7.2973525664e-3;
fine_inv = one / fine; % = 137.035999139;

% classical electron radius
r_e = q_e2 / m_e_erg; % = 2.8179403227e-13; % [cm]

% electron Compton wavelength, reduced
lambdabar_e = r_e * fine_inv; % = 3.8615926764e-11 [cm];

% electron Compton wavelength
lambda_e = lambdabar_e * twpi; % = 2.4263102367e-10 [cm];

% Bohr radius
a_0 = r_e * fine_inv * fine_inv; % = 5.2917721067e-9; % [cm]
r_Bohr = a_0;

% Rydberg constant
R_inf = frpi_inv * fine / a_0; % = 1.0973731568508e+5; % [cm^{-1}]

% Rydberg energy unit
Ry_cm_inv = R_inf; % = 1.0973731568508e+5; % [cm^{-1}]
Ry_erg = R_inf * h * c_light; % = 2.179872325e-11; % [erg]
Ry_eV = Ry_erg * erg2eV; % = 13.605693009; % [eV]
Ry = Ry_eV;

% Stefan-Boltzmann constant
sigma_SB = two15th * pi5 * k_B4 * h3_inv * c_light2_inv; % = 5.670367e-5; % [erg s^{-1} cm^{-2} K^{-4}]

% radiation constant
a_rad = four * sigma_SB * c_light_inv; % = 7.565723e-15; % [erg cm^{-3} K^{-4}]
a_rad_eV = a_rad * erg2eV * eV2K * eV2K * eV2K * eV2K; % = 8.563456e+13; % [eV cm^{-3} eV^{-4}]
a_rad_MeV = a_rad_eV * MeV2eV * MeV2eV * MeV2eV * fm32cm3; % = 8.563456e-8; % [MeV fm^{-3} MeV^{-4}]


%% electro-weak constants
%   C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
%     http://pdg.lbl.gov/2016/reviews/rpp2016-rev-phys-constants.pdf

% weak coupling constant
G_F = 1.1663787e-5; % [GeV^{-2}]
G_F_MeV = G_F * MeV2GeV * MeV2GeV; % [MeV^{-2}]

% W+- boson mass
m_W = 80.385; % [GeV]
m_W_MeV = m_W * GeV2MeV; % [MeV]

% Z0 boson mass
m_Z = 91.1876; % [GeV]
m_Z_MeV = m_Z * GeV2MeV; % [MeV]

% weak mixing angle
%sin2W = 0.22336;
sin2W = 0.23129;
%sin2W = 0.23148;
%sin2W = 0.23865;
%sin2W = 0.23152;

% Gamow-Teller beta decay constant
% anomalous proton magnetic moment
% anomalous neutron magnetic moment
% vector current neutrino-proton coupling constant
% vector current neutrino-neutron coupling constant
% axial vector current neutrino-proton coupling constant
% axial vector current neutrino-neutron coupling constant
% neutron and proton form factors
% neutrino number and energy density and flux coefficients


%% astrophysical constants
%   C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
%     http://pdg.lbl.gov/2016/reviews/rpp2016-rev-astrophysical-constants.pdf

% astronomical unit
au = 1.4959787070000e+13; % [cm]
au2 = au * au;
au3 = au * au2;

% parsec
pc = 3.08567758149e+18; % [cm]

% sidereal year
yr = 3.15581498e+7; % [s]
yr_inv = one / yr;
yr2_inv = yr_inv * yr_inv;

% mean sidereal day
day = 8.616409053e+4; % [s]

% Solar mass parameter; heliocentric gravitational constant
GM_solar = 1.3271244e+26; % [cm^3 s^{-2}]

% Solar mass
M_solar = GM_solar * G_N_inv; % = 1.98848e+33; % [g]

% nominal Solar equatorial radius
R_solar = 6.957e+10; % [cm]

% Earth mass parameter; geocentric gravitational constant
GM_earth = 3.986004e+20; % [cm^3 s^{-2}]

% Earth mass
M_earth = GM_earth * G_N_inv; % = 5.9724e+27; % [g]

% nominal Earth equatorial radius
R_earth = 6.3781e+8; % [cm]