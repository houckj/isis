% These values are in cgs system
% 1998 CODATA recommended values

% The following are exact
variable Const_c = 2.99792458e10;	% [cm/s] speed of light
%

variable Const_h = 6.62606876e-27;	% [erg s] Planck const
variable Const_G = 6.673e-8;   		% [cm^3/g/s^2] Newtonian const
variable Const_m_e = 9.10938188e-28;	% [g] electron mass
variable Const_m_p = 1.67262158e-24;	% [g] proton mass
variable Const_k = 1.3806503e-16;	% [erg/K] Boltzmann const

% Conversion factors
variable Const_eV = 1.602176462e-12;	% [erg]
variable Const_u = 1.66053873e-24;	% [g] atomic mass unit

% hc in [keV Angstrom] 
variable Const_keV_A = Const_h * Const_c * 1.e8 / (1.e3 * Const_eV);
variable Const_hc = Const_keV_A;

