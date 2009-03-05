%; Time-stamp: <2001-10-23 15:18:39 dph> 
%; MIT Directory: ~dph/libisis
%; File: apec_contin_fit.sl
%; Author: D. Huenemoerder
%; Original version: 2001.10.19
%;====================================================================
% version: 1.0
%
% Name: aped_contin_fit( xlo, xhi, par )
%
% Purpose: define an apec model which computes the apec continuum
%	   spectrum for 1 temperature component.  Does not
%	   support continuua by ion or for variable abundances
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Name: fast_aped_contin_fit( xlo, xhi, par )
%
% Purpose: Use a global continuum model variable for
%	   a re-normalizable continuum shape model.  If the norm is
%	   frozen, then it is an absolute continuum model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Name: fast_aped_contin_init(nhist)
%
% Purpose: Read a continuum model (from a fit) and store into a global
%	   variable for use by fast_aped_contin_fit().
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
static variable Contin_X_lo, Contin_X_hi, Contin_Value;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
define fast_aped_contin_init (nhist)
{  %+ define
   % ASSUMES that a continuum model has been evaluated for histogram
   % nhist.  This should be called with an MEG histogram, in order to
   % provide enough bandpass for both HEG and MEG.
   % Resolution is not very important, since it is a continuum.

% Setup is not done here.  First, the equivalent of the following
% should be done:
%%   fit_fun("aped_contin(1) + aped_contin(2) + ... + aped_contin(n)");
%%   set_par(...);
%%   fit_counts;
%%   % etc.
%% Prepare to grab the final continuum array:
%%   ignore([lo:hi]);
%%   notice(hist);
%%   group_data(nhist,0);     % undo any grouping or binning
%%   rebin_data(nhist,0);
%%   eval_flux;		    % re-evaluate flux, unbinned.

% Here we  get the continuum array and set global variables:
   variable y;
   y=get_model_flux(nhist);
   Contin_Value = y.value;
   Contin_X_lo = y.bin_lo;
   Contin_X_hi = y.bin_hi;
%
}   %- define

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now we simply rebin the Contin_Value to the specified grid,
%% allowing a re-normalization:
%
define fast_aped_contin_fit (xlo, xhi, par)
{  %+ define
   variable y =
      rebin (xlo, xhi, Contin_X_lo, Contin_X_hi, Contin_Value) * par[0];
   return y;
}   %- define


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define an aped continuum function.
% Proper use of this requires ignoring line features.
%
define aped_contin_fit ( xlo, xhi, par )
{ %+ define
  % xlo, xhi = grid of wavelength bins;
  % par[0] = norm;              10^{-14} * VEM/ (4piD^2)
  %				10^{-14}~ int{dV n_e n_h} / (4piD^2)
  % par[1] = electron temperature [K];

  variable y;
  variable KNORM=1.e14;

  y = (get_contin(xlo, xhi, par[1])).pseudo;
  y = (y + (get_contin(xlo, xhi, par[1])).true) * par[0]*KNORM;
    
  return y;

} %- define

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

add_slang_function("aped_contin", ["norm","T"]);
add_slang_function ("fast_aped_contin", ["norm"]);




