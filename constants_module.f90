! Fundamental, astrophysical, and numerical constants

module constants_module

  implicit none
  public

  ! kind parameters for floating-point (real) number precision
  !   if 128-bit (qp) not supported, use 64-bit (dp)
  integer, parameter :: sp = selected_real_kind(6,37)
  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: qp = max( selected_real_kind(33,4931), dp )
  ! integer, parameter :: sp = real32
  ! integer, parameter :: dp = real64
  ! integer, parameter :: qp = max( real128, real64 )

  ! kind parameters for integer precision
  !   if 64-bit (i8) not supported, use 32-bit (i4)
  integer, parameter :: i2 = selected_int_kind(4)
  integer, parameter :: i4 = selected_int_kind(9)
  integer, parameter :: i8 = selected_int_kind(18)
  ! integer, parameter :: i2 = int16
  ! integer, parameter :: i4 = int32
  ! integer, parameter :: i8 = int64


  ! numerical constants

  ! integers
  real (dp), parameter :: zero      = 0.0_dp
  real (dp), parameter :: one       = 1.0_dp
  real (dp), parameter :: two       = 2.0_dp
  real (dp), parameter :: three     = 3.0_dp
  real (dp), parameter :: four      = 4.0_dp
  real (dp), parameter :: five      = 5.0_dp
  real (dp), parameter :: six       = 6.0_dp
  real (dp), parameter :: seven     = 7.0_dp
  real (dp), parameter :: eight     = 8.0_dp
  real (dp), parameter :: nine      = 9.0_dp
  real (dp), parameter :: ten       = 10.0_dp
  real (dp), parameter :: eleven    = 11.0_dp
  real (dp), parameter :: twelve    = 12.0_dp
  real (dp), parameter :: fifteen   = 15.0_dp
  real (dp), parameter :: sixteen   = 16.0_dp
  real (dp), parameter :: sixty     = 60.0_dp

  ! fractions
  real (dp), parameter :: half      = 0.5_dp
  real (dp), parameter :: third     = one / three
  real (dp), parameter :: fourth    = 0.25_dp
  real (dp), parameter :: fifth     = 0.2_dp
  real (dp), parameter :: sixth     = one / six
  real (dp), parameter :: seventh   = one / seven
  real (dp), parameter :: eighth    = 0.125_dp
  real (dp), parameter :: ninth     = one / nine
  real (dp), parameter :: tenth     = 0.1_dp
  real (dp), parameter :: twelfth   = one / twelve
  real (dp), parameter :: fifteenth = one / fifteen
  real (dp), parameter :: sixteenth = 0.0625_dp
  real (dp), parameter :: sixtieth  = one / sixty
  real (dp), parameter :: two3rd    = two * third
  real (dp), parameter :: four3rd   = four * third
  real (dp), parameter :: five3rd   = five * third
  real (dp), parameter :: three4th  = 0.75_dp
  real (dp), parameter :: two5th    = 0.4_dp
  real (dp), parameter :: three5th  = 0.6_dp
  real (dp), parameter :: four5th   = 0.8_dp
  real (dp), parameter :: five6th   = five * sixth
  real (dp), parameter :: three8th  = 0.375_dp
  real (dp), parameter :: five8th   = 0.625_dp
  real (dp), parameter :: seven8th  = 0.875_dp
  real (dp), parameter :: five12th  = five * twelfth
  real (dp), parameter :: seven12th = seven * twelfth
  real (dp), parameter :: two15th   = two * fifteenth
  real (dp), parameter :: eight15th = eight * fifteenth

  ! 71 digits of pi (256-bit precision)
  real (dp), parameter :: pi        = 3.1415926535897932384626433832795028841971693993751058209749445923078164_dp
  real (dp), parameter :: pi2       = pi * pi
  real (dp), parameter :: pi3       = pi * pi2
  real (dp), parameter :: pi4       = pi * pi3
  real (dp), parameter :: pi5       = pi * pi4
  real (dp), parameter :: twpi      = pi * two
  real (dp), parameter :: thpi      = pi * three
  real (dp), parameter :: frpi      = pi * four
  real (dp), parameter :: etpi      = pi * eight
  real (dp), parameter :: sxtnpi    = pi * sixteen
  real (dp), parameter :: frpith    = frpi * third
  real (dp), parameter :: thpifrth  = thpi * fourth
  real (dp), parameter :: pi_inv    = one / pi
  real (dp), parameter :: twpi_inv  = pi_inv * half
  real (dp), parameter :: frpi_inv  = pi_inv * fourth
  real (dp), parameter :: etpi_inv  = pi_inv * eighth

  ! logarithms
  real (dp), parameter :: log_e     = log10( exp(one) )
  real (dp), parameter :: ln_2      = log( two )
  real (dp), parameter :: ln_10     = log( ten )


  ! fundamental physical constants (CODATA recommended 2014 values)
  !   P. J. Mohr, D. B. Newell, and B. N. Taylor, Rev. Mod. Phys. 88, 035009 (2016)
  !   http://physics.nist.gov/cuu/Constants/Table/allascii.txt
  ! * whenever possible, derived constants are calculated explicitly for floating-point consistency

  ! common unit conversion factors
  real (dp), parameter :: eV2MeV = 1.0e-6_dp
  real (dp), parameter :: MeV2eV = 1.0e+6_dp

  real (dp), parameter :: eV2GeV = 1.0e-9_dp
  real (dp), parameter :: GeV2eV = 1.0e+9_dp

  real (dp), parameter :: MeV2GeV = 1.0e-3_dp
  real (dp), parameter :: GeV2MeV = 1.0e+3_dp

  real (dp), parameter :: erg2J = 1.0e-7_dp
  real (dp), parameter :: J2erg = 1.0e+7_dp

  real (dp), parameter :: cm2fm = 1.0e+13_dp
  real (dp), parameter :: fm2cm = 1.0e-13_dp

  real (dp), parameter :: cm32fm3 = cm2fm * cm2fm * cm2fm
  real (dp), parameter :: fm32cm3 = fm2cm * fm2cm * fm2cm

  ! speed of light in vacuum
  real (dp), parameter :: c_light = 2.99792458e+10_dp ! [cm s^{-1}]
  real (dp), parameter :: c_light_inv = one / c_light
  real (dp), parameter :: c_light2 = c_light * c_light
  real (dp), parameter :: c_light2_inv = one / c_light2
  real (dp), parameter :: c_light3 = c_light * c_light2
  real (dp), parameter :: c_light3_inv = one / c_light3

  ! common mass-energy conversion factors
  real (dp), parameter :: g2erg = c_light2 ! = 8.987551787368176e+20
  real (dp), parameter :: erg2g = c_light2_inv ! = 1.112650056053618e-21

  ! elementary charge
  real (dp), parameter :: q_e_C = 1.6021766208e-19_dp ! [C]
  real (dp), parameter :: q_e = q_e_C * tenth * c_light ! = 4.80320467e-10 ! [esu] or [erg^{1/2} cm^{1/2}]
  real (dp), parameter :: q_e2 = q_e * q_e

  ! common electronvolt conversion factors
  real (dp), parameter :: eV2J = q_e_C ! = 1.6021766208e-19
  real (dp), parameter :: J2eV = one / eV2J ! = 6.241509126e+18

  real (dp), parameter :: eV2erg = eV2J * J2erg ! = 1.6021766208e-12
  real (dp), parameter :: erg2eV = erg2J * J2eV ! = 6.241509126e+11

  real (dp), parameter :: eV2g = eV2erg * erg2g ! = 1.782661907e-33
  real (dp), parameter :: g2eV = g2erg * erg2eV ! = 5.609588650e+32

  ! Planck constant
  real (dp), parameter :: h = 6.626070040e-27_dp ! [erg s]
  real (dp), parameter :: h_MeV = h * erg2eV * eV2MeV ! = 4.135667662e-21 ! [MeV s]
  real (dp), parameter :: h_inv = one / h
  real (dp), parameter :: h2 = h * h
  real (dp), parameter :: h2_inv = one / h2
  real (dp), parameter :: h3 = h * h2
  real (dp), parameter :: h3_inv = one / h3

  ! common inverse length conversion factors
  real (dp), parameter :: cm_inv2g = h * c_light_inv ! = 2.210219057e-37
  real (dp), parameter :: g2cm_inv = c_light * h_inv ! = 4.524438411e+38

  real (dp), parameter :: cm_inv2erg = h * c_light ! = 1.986445824e-16
  real (dp), parameter :: erg2cm_inv = h_inv * c_light_inv ! = 5.034116651e+15

  real (dp), parameter :: cm_inv2eV = cm_inv2erg * erg2eV ! = 1.2398419739e-4
  real (dp), parameter :: eV2cm_inv = eV2erg * erg2cm_inv ! = 8.065544005e+3

  ! Planck constant, reduced; Dirac constant
  real (dp), parameter :: hbar = twpi_inv * h ! = 1.054571800e-27 ! [erg s] 
  real (dp), parameter :: hbar_MeV = hbar * erg2eV * eV2MeV ! = 6.582119514e-22 ! [MeV s]
  real (dp), parameter :: hbar_inv = one / hbar
  real (dp), parameter :: hbar2 = hbar * hbar
  real (dp), parameter :: hbar2_inv = one / hbar2
  real (dp), parameter :: hbar3 = hbar * hbar2
  real (dp), parameter :: hbar3_inv = one / hbar3

  ! common conversion factor
  real (dp), parameter :: hbarc = hbar * c_light ! [erg cm] 
  real (dp), parameter :: hbarc_MeV = hbar_MeV * c_light * cm2fm ! = 197.3269788 ! [MeV fm] 
  real (dp), parameter :: hbarc_inv = one / hbarc

  ! Newtonian constant of gravitation
  real (dp), parameter :: G_N = 6.67408e-8_dp ! [cm^3 g^-1 s^-2]
  real (dp), parameter :: G_N_inv = one / G_N

  ! Boltzmann constant
  real (dp), parameter :: k_B = 1.38064852e-16_dp ! [erg K^{-1}]
  real (dp), parameter :: k_B_MeV = k_B * erg2eV * eV2MeV ! = 8.6173303e-11 ! [MeV K^{-1}]
  real (dp), parameter :: k_B_inv = one / k_B
  real (dp), parameter :: k_B4 = k_B * k_B * k_B * k_B
  real (dp), parameter :: k_B4_MeV = k_B_MeV * k_B_MeV * k_B_MeV * k_B_MeV

  ! common temperature-energy conversion factors
  real (dp), parameter :: K2erg = k_B ! = 1.38064852e-16
  real (dp), parameter :: erg2K = one / K2erg ! = 7.2429731e+15

  real (dp), parameter :: K2eV = k_B * erg2eV ! = 8.6173303e-5
  real (dp), parameter :: eV2K = one / K2eV ! = 1.16045221e+4

  ! Avogadro constant
  real (dp), parameter :: N_A = 6.022140857e+23_dp ! [mol^{-1}]

  ! atomic mass unit
  real (dp), parameter :: m_u = one / N_A ! = 1.660539040e-24 ! [g]
  real (dp), parameter :: m_u_eV = m_u * g2eV ! = 931.4940954e+6 ! [eV]
  real (dp), parameter :: m_u_MeV = m_u_eV * eV2MeV ! = 931.4940954 ! [MeV]

  ! common atomic-mass-unit conversion factors
  real (dp), parameter :: u2g = m_u ! = 1.660539040e-24
  real (dp), parameter :: g2u = N_A ! = 6.022140857e+23

  real (dp), parameter :: u2eV = u2g * g2eV ! = 931.4940954e+6 
  real (dp), parameter :: eV2u = eV2g * g2u ! = 1.0735441105e-9 

  ! electron mass
  real (dp), parameter :: m_e = 9.10938356e-28_dp ! [g]
  real (dp), parameter :: m_e_erg = m_e * g2erg ! = 8.18710565e-7 ! [erg]
  real (dp), parameter :: m_e_eV = m_e * g2eV ! = 0.5109989461e+6 ! [eV]
  real (dp), parameter :: m_e_MeV = m_e_eV * eV2MeV ! = 0.5109989461 ! [MeV]
  real (dp), parameter :: m_e_u = m_e * g2u ! = 5.48579909070e-4 ! [u]

  ! proton mass
  real (dp), parameter :: m_p = 1.672621898e-24_dp ! [g]
  real (dp), parameter :: m_p_erg = m_p * g2erg ! = 1.503277593e-3 ! [erg]
  real (dp), parameter :: m_p_eV = m_p * g2eV ! = 938.2720813e+6 ! [eV]
  real (dp), parameter :: m_p_MeV = m_p_eV * eV2MeV ! = 938.2720813 ! [MeV]
  real (dp), parameter :: m_p_u = m_p * g2u ! = 1.007276466879 ! [u]

  ! neutron mass
  real (dp), parameter :: m_n = 1.674927471e-24_dp ! [g]
  real (dp), parameter :: m_n_erg = m_n * g2erg ! = 1.505349739e-3 ! [erg]
  real (dp), parameter :: m_n_eV = m_n * g2eV ! = 939.5654133e+6 ! [eV]
  real (dp), parameter :: m_n_MeV = m_n_eV * eV2MeV ! = 939.5654133 ! [MeV]
  real (dp), parameter :: m_n_u = m_n * g2u ! = 1.00866491588 ! [u]

  ! fine-structure constant
  real (dp), parameter :: fine = q_e2 * hbarc_inv ! = 7.2973525664e-3
  real (dp), parameter :: fine_inv = one / fine ! = 137.035999139

  ! classical electron radius
  real (dp), parameter :: r_e = q_e2 / m_e_erg ! = 2.8179403227e-13 ! [cm]

  ! electron Compton wavelength, reduced
  real (dp), parameter :: lambdabar_e = r_e * fine_inv ! = 3.8615926764e-11 [cm]

  ! electron Compton wavelength
  real (dp), parameter :: lambda_e = lambdabar_e * twpi ! = 2.4263102367e-10 [cm]

  ! Bohr radius
  real (dp), parameter :: a_0 = r_e * fine_inv * fine_inv ! = 5.2917721067e-9 ! [cm]
  real (dp), parameter :: r_Bohr = a_0

  ! Rydberg constant
  real (dp), parameter :: R_inf = frpi_inv * fine / a_0 ! = 1.0973731568508e+5 ! [cm^{-1}]

  ! Rydberg energy unit
  real (dp), parameter :: Ry_cm_inv = R_inf ! = 1.0973731568508e+5 ! [cm^{-1}]
  real (dp), parameter :: Ry_erg = R_inf * h * c_light ! = 2.179872325e-11 ! [erg]
  real (dp), parameter :: Ry_eV = Ry_erg * erg2eV ! = 13.605693009 ! [eV]
  real (dp), parameter :: Ry = Ry_eV

  ! Stefan-Boltzmann constant
  real (dp), parameter :: sigma_SB = two15th * pi5 * k_B4 * h3_inv * c_light2_inv ! = 5.670367e-5 ! [erg s^{-1} cm^{-2} K^{-4}]

  ! radiation constant
  real (dp), parameter :: a_rad = four * sigma_SB * c_light_inv ! = 7.565723e-15 ! [erg cm^{-3} K^{-4}]
  real (dp), parameter :: a_rad_eV = a_rad * erg2eV * eV2K * eV2K * eV2K * eV2K ! = 8.563456e+13 ! [eV cm^{-3} eV^{-4}]
  real (dp), parameter :: a_rad_MeV = a_rad_eV * MeV2eV * MeV2eV * MeV2eV * fm32cm3 ! = 8.563456e-8 ! [MeV fm^{-3} MeV^{-4}]


  ! electro-weak constants
  !   C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
  !     http://pdg.lbl.gov/2016/reviews/rpp2016-rev-phys-constants.pdf

  ! weak coupling constant
  real (dp), parameter :: G_F = 1.1663787e-5_dp ! [GeV^{-2}]
  real (dp), parameter :: G_F_MeV = G_F * MeV2GeV * MeV2GeV ! [MeV^{-2}]

  ! W+- boson mass
  real (dp), parameter :: m_W = 80.385_dp ! [GeV]
  real (dp), parameter :: m_W_MeV = m_W * GeV2MeV ! [MeV]

  ! Z0 boson mass
  real (dp), parameter :: m_Z = 91.1876_dp ! [GeV]
  real (dp), parameter :: m_Z_MeV = m_Z * GeV2MeV ! [MeV]

  ! weak mixing angle
  !real (dp), parameter :: sin2W = 0.22336
  real (dp), parameter :: sin2W = 0.23129
  !real (dp), parameter :: sin2W = 0.23148
  !real (dp), parameter :: sin2W = 0.23865
  !real (dp), parameter :: sin2W = 0.23152

  ! Gamow-Teller beta decay constant
  ! anomalous proton magnetic moment
  ! anomalous neutron magnetic moment
  ! vector current neutrino-proton coupling constant
  ! vector current neutrino-neutron coupling constant
  ! axial vector current neutrino-proton coupling constant
  ! axial vector current neutrino-neutron coupling constant
  ! neutron and proton form factors
  ! neutrino number and energy density and flux coefficients


  ! astrophysical constants
  !   C. Patrignani et al. (Particle Data Group), Chin. Phys. C, 40, 100001 (2016)
  !     http://pdg.lbl.gov/2016/reviews/rpp2016-rev-astrophysical-constants.pdf

  ! astronomical unit
  real (dp), parameter :: au = 1.4959787070000e+13_dp ! [cm]
  real (dp), parameter :: au2 = au * au
  real (dp), parameter :: au3 = au * au2

  ! parsec
  real (dp), parameter :: pc = 3.08567758149e+18_dp ! [cm]

  ! sidereal year
  real (dp), parameter :: yr = 3.15581498e+7_dp ! [s]
  real (dp), parameter :: yr_inv = one / yr
  real (dp), parameter :: yr2_inv = yr_inv * yr_inv

  ! mean sidereal day
  real (dp), parameter :: day = 8.616409053e+4_dp ! [s]

  ! Solar mass parameter; heliocentric gravitational constant
  real (dp), parameter :: GM_solar = 1.3271244e+26_dp ! [cm^3 s^{-2}]

  ! Solar mass
  real (dp), parameter :: M_solar = GM_solar * G_N_inv ! = 1.98848e+33 ! [g]

  ! nominal Solar equatorial radius
  real (dp), parameter :: R_solar = 6.957e+10_dp ! [cm]

  ! Earth mass parameter; geocentric gravitational constant
  real (dp), parameter :: GM_earth = 3.986004e+20_dp ! [cm^3 s^{-2}]

  ! Earth mass
  real (dp), parameter :: M_earth = GM_earth * G_N_inv ! = 5.9724e+27 ! [g]

  ! nominal Earth equatorial radius
  real (dp), parameter :: R_earth = 6.3781e+8_dp ! [cm]

 
end module constants_module