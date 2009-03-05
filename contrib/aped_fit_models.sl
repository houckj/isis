%; Time-stamp: <2002-05-15 10:53:37 dph> 
%; MIT Directory: ~dph/libisis
%; File: aped_fit_models.sl
%; Author: D. Huenemoerder
%; Original version: 2001.10.23
%;====================================================================
%
% purpose:  Define isis user models for aped thermal plasmas.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example scenarios:
%

% Single-temperature plasma, reduced iron:
%
% fit_fun("apedx(1)");
% set_par("apedx(1).norm", 0.001, 0, 0.0, 1.0);
% set_par("apedx(1).temp", 10.^7.2, 0, 1.e5, 8.e8);
% list_par;
% eval_counts;


% Multi-thermal plasma:
%
%  multi_Temps = 10.^[6.5:8.0:0.1];        %  temperature grid
%  multi_wgts  = Float_Type[ length(multi_Temps) ] + 0.001 ;   % weights grid
%  elem_list = [Fe, O, Ne];
%  elem_abund = [0.1, 0.5, 3.0];
%  aped_set_abund( elem_list, elem_abund );
%
% aped_multi_setup( multi_Temps, multi_wgts );
% list_par;      % a zillion....
% eval_counts;   % wait a while...
% plot_model_counts(n);

%%%%%%%%%
% Implicit line and continuum models (modal behavior!)
%  ... continuing w/ above multi-thermal...
% Aped_Filter.type = MODEL_LINES;
% Aped_Filter.lines = where(el_ion(Fe,17));
% (x1, x2) = linear_grid (1, 1+8191*0.005, 8192);
% y_Fe17 = eval_fun( x1, x2 );
% Aped_Filter.type = MODEL_CONTIN;
% y_cont = eval_fun( x1, x2 );
% Aped_Filter.type = NULL;             % RESET



%%%%%%%%%%
% Explicit continuum models:
%
% Single temperature continuum, low oxygen abundance:
%
% fit_fun("aped_contin(1)");
% set_par("aped_contin(1).norm", 0.001, 0, 0.0, 1.0);
% set_par("aped_contin(1).T",  10.^7.2, 0, 1.e5, 8.e8);
%  aped_set_abund( O, 0.0 );

% Multi temperature continuum, normal abundances:
%
% fit_fun("aped_multi_contin(1)");
%  aped_set_abund();
%  aped_multi_contin_setup( 10.^[6.8, 7.2, 7.6], [1.0, 10.0, 5.0] );
%  eval_counts;
%
% Cache that and define the fast version w/ normalization parameter only:
%
% ignore(all_data) ; 
% xnotice( n_MEG, 1.7, 25 );  % notice an MEG histogram
% group_data( n_MEG , 0 ) ;   % highest resolution grid.
% eval_flux;                  % compute the model
% fast_aped_contin_init( n_MEG );     % cache it.
% fit_fun("fast_aped_contin(1) + gauss(1)");  % define line + continuum model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% version 4.0 2002.05.15  Added revised continuum models.
%                         Unify parameter names:  norm, temp
%
% version 3.2 2002.05     fixed aped_fit - collided w/ aped the
%                          structure; renamed apedx_fit.
%
% version 3.1 2002.04.03  added aped_fit - basic model, one component.
%
% version: 3.0 2002.01.06 add global variable and switch to support 
%                         qualifier of model_spectrum, and use for
%                         evaluating a model for specific lines or
%                         contin. 
%
% version: 2.0 2001.11.29 include rv version of model
%                fixed bug in aped_multi_fit, w/ array counting/indexing.
% version: 1.1 2001.10.29 added is_a_norm flags to
%			  add_slang_function() call
% version: 1.0
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   public variable Aped_Filter = struct {type, lines};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
provide("apedx_fit");
provide("aped_multi_fit");
provide("aped_multi_init");
provide("aped_multi_setup");
provide("aped_set_abund");

provide("aped_rv_fit");
provide("aped_multi_rv_fit");
provide("aped_multi_rv_init");
provide("aped_multi_rv_setup");
provide("aped_multi_2rv_setup");

provide("aped_contin_fit");
provide("aped_multi_contin_fit");
provide("fast_aped_contin_fit");
provide("fast_aped_contin_init");
provide("aped_multi_contin_setup");
provide("aped_multi_contin_init");


% pre-define functions:
%
define aped_multi_fit();       % multple T isis model
define aped_multi_init();      % defines the model for given # components
define aped_multi_setup();     % establishes fit function, sets the
			       %  default model parameters
define aped_set_abund();       % sets/resets local static abundance arrays.
%
% ditto, but w/ radial velocity as a parameter (single for all components.) 
%
define aped_multi_rv_fit();       % multple T isis model
define aped_multi_rv_init();      % defines the model for given # components
define aped_multi_rv_setup();     % establishes fit function, sets the
	 		          %  default model parameters

define aped_contin_fit();
define aped_multi_contin_fit();
define fast_aped_contin_fit();
define fast_aped_contin_init();
define aped_multi_contin_setup();
define aped_multi_contin_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% define local variables for use in storing alternative abundance
% mixes for the plasma models.  These are not used as parameters of a
% fit, but for ad hoc tweaks of abundances. See aped_set_abund();
%
%
static variable Aped_Elem_List,  % an array of Z values; e.g [Fe, Ne] or [26, 10]
		Aped_Abund_List; % an array of factores relative to
				 % cosmic, e.g., [0.1, 2.0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
define apedx_fit ( xlo, xhi, par )
{
%
% PURPOSE:  isis user model, APED thermal plasma; simple interface to
% the plasma state structure.
%	    
% xlo, xhi = grid of wavelength bins;
% 
% params are 
%
%    norm = 1
%    temperature = 1e+07
%    vturb = 0
%    vrad = 0
%

  variable ii,                        	     % loop counter
	   y;                         	     % return value

  variable mdl = default_plasma_state();  % plasma model definition

% NOTE: density, absorption not currently supported....

    mdl.norm = par[0];     
    mdl.temperature = par[1];
    mdl.vturb = par[2];              % given in km/s?
    mdl.redshift = par[3]/2.99792458e5;   % param is in velocity, km/s.

    if ( __is_initialized( &Aped_Elem_List ) )   % apply different abundances, if set
    {
        mdl.elem = Aped_Elem_List;
        mdl.elem_abund = Aped_Abund_List;
    }

  define_model(mdl);                      % define the model (use list_model; to see).

   switch (Aped_Filter.type)
   {
       case MODEL_LINES:
       if (Aped_Filter.lines == NULL)
         y = model_spectrum (xlo, xhi, MODEL_LINES);
       else
         y = model_spectrum (xlo, xhi, MODEL_LINES, Aped_Filter.lines);
   }
   {
       case MODEL_CONTIN:
       y = model_spectrum (xlo, xhi, MODEL_CONTIN);
   }
   {
       case NULL:
       y = model_spectrum (xlo, xhi);
   }
   
  return y;

}  
%-

variable params = ["norm", "temp", "vturb", "vrad"];
add_slang_function( "apedx", params );




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
define aped_multi_fit ( xlo, xhi, par )
{
%
% PURPOSE:  isis user model, APED thermal plasma, multi-temperature,
%	    cosmic abundances, low density limit.

% xlo, xhi = grid of wavelength bins;
% 

  variable ii,                        	     % loop counter
	   y,                         	     % return value
           ncomps = length(par)/2,    	     % # model components
	   pstate = default_plasma_state();  % plasma model definition

  variable mdl = Struct_Type[ncomps];
  variable norms=par[ [0:ncomps-1] ];
  variable temps=par[ [ncomps:2*ncomps-1] ];	   	   

  for (ii = 0; ii < ncomps; ii++ )
  {
    mdl[ii] =             @pstate;
    mdl[ii].temperature = temps[ii];
    mdl[ii].norm =        norms[ii];

    if ( __is_initialized( &Aped_Elem_List ) )   % apply different abundances, if set
    {
        mdl[ii].elem = Aped_Elem_List;
        mdl[ii].elem_abund = Aped_Abund_List;
    }

  }

  define_model(mdl);                      % define the model (use list_model; to see).

   switch (Aped_Filter.type)
   {
       case MODEL_LINES:
       if (Aped_Filter.lines == NULL)
         y = model_spectrum (xlo, xhi, MODEL_LINES);
       else
         y = model_spectrum (xlo, xhi, MODEL_LINES, Aped_Filter.lines);
   }
   {
       case MODEL_CONTIN:
       y = model_spectrum (xlo, xhi, MODEL_CONTIN);
   }
   {
       case NULL:
       y = model_spectrum (xlo, xhi);
   }
   
  return y;

}  
%-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
define aped_multi_init ( ncomps )
{
% Dynamically define a model w/ ncomps temperature compents.

  variable i, func_name, T_params, Norm_params, params;

  params = "norm_1";

  for (i=2; i <= ncomps; i++)
    params = [params, "norm_"+string(i)];

  for (i=1; i <= ncomps; i++)
    params = [params, "temp_"+string(i)];

  add_slang_function( "aped_multi", params, [1:ncomps] );

}
%-



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+

define aped_multi_setup(Temps, Norms)
{
  variable i;

  aped_multi_init( length(Temps) );

  fit_fun("aped_multi(1)");

  for ( i = 0; i < length(Temps); i++)
  {
    set_par( i+1, Norms[i], 0, 0, 100*Norms[i] );
    set_par( i+length(Temps)+1, Temps[i], 0, 1.e4, 7.95e8 );
  }

%  list_par;
}  

%-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+

define aped_set_abund()       % sets local static abundance arrays.
{
  variable  elems;

% example usage:
% aped_set_abund([Fe, Ne], [0.1, 2.0]);   % change abunds of iron and  neon
% aped_set_abund([Fe, Ne]);   % reset abunds of iron and  neon to  cosmic 
% aped_set_abund();           % reset current list to cosmic.

  switch (_NARGS)
  {
    case 0:           % reset list to cosmic (1.0)
    if (__is_initialized (&Aped_Abund_List) )
      Aped_Abund_List = Aped_Abund_List*0 + 1;     % set to cosmic;
  }
  {
    case 1:   % have element list, so set them to cosmic.
    elems=();
    variable i;
    for ( i = 0;  i < length(elems); i++)
      Aped_Abund_List[where(Aped_Elem_List == elems[i])] = 1.;    % reset list to 1
  }
  {
    case 2:
    Aped_Abund_List = ();
    Aped_Elem_List = ();
  }
  {
    vmessage("Bad arguments.  USAGE: aped_set_abund([elems [,abunds]]);");
  }
}
%-



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+  Version w/ radial velocity delta.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
define aped_multi_rv_fit ( xlo, xhi, par )
{
%
% PURPOSE:  isis user model, APED thermal plasma, multi-temperature,
%	    cosmic abundances, low density limit, w/ velocity shift 

% xlo, xhi = grid of wavelength bins;
% 

  variable ii,                        	     % loop counter
	   y,                         	     % return value
           ncomps = (length(par)-1) / 2,     % # model temperature components
	   pstate = default_plasma_state();  % plasma model definition

  variable rv = par[0];               % radial velocity, km/s
  variable mdl = Struct_Type[ncomps];
  variable norms = par[ [1 : ncomps] ];  % normalizations, 1.e-14*EM/4piD^2
  variable temps=par[ [1+ncomps : 2*ncomps]];  % temperatures, K

  for (ii = 0; ii < ncomps; ii++ )
  {
    mdl[ii] =             @pstate;
    mdl[ii].temperature = temps[ii];
    mdl[ii].norm =        norms[ii];
    mdl[ii].redshift =    rv / 2.99792458e5;

    if ( __is_initialized( &Aped_Elem_List ) )   % apply different abundances, if set
    {
        mdl[ii].elem = Aped_Elem_List;
        mdl[ii].elem_abund = Aped_Abund_List;
    }

  }

  define_model(mdl);                      % define the model (use list_model; to see).

   switch (Aped_Filter.type)
   {
       case MODEL_LINES:
       if (Aped_Filter.lines == NULL)
         y = model_spectrum (xlo, xhi, MODEL_LINES);
       else
         y = model_spectrum (xlo, xhi, MODEL_LINES, Aped_Filter.lines);
   }
   {
       case MODEL_CONTIN:
       y = model_spectrum (xlo, xhi, MODEL_CONTIN);
   }
   {
       case NULL:
       y = model_spectrum (xlo, xhi);
   }

  return y;

}  
%-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
define aped_multi_rv_init ( ncomps )
{
% Dynamically define a model w/ ncomps temperature components.

  variable i, func_name, T_params, Norm_params, params;

  params = "vrad";

  for (i=1; i <= ncomps; i++)
    params = [params, "norm_"+string(i)];

  for (i=1; i <= ncomps; i++)
    params = [params, "temp_"+string(i)];

  add_slang_function( "aped_multi_rv", params, [2:ncomps+1] );

}
%-



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+

define aped_multi_rv_setup(Vrad, Temps, Norms)
{
  variable i;

  aped_multi_rv_init( length(Temps) );

  fit_fun("aped_multi_rv(1)");

  set_par( 1, Vrad, 0, Vrad-600, Vrad+600);

  for ( i = 0; i < length(Temps); i++)
  {
    set_par( i+2, Norms[i], 0, 0, 100*Norms[i] );
    set_par( i+length(Temps)+2, Temps[i], 0, 1.e4, 7.95e8 );
  }

%  list_par;
}  

%-


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+

define aped_multi_2rv_setup(Vrads, Temps, Norms1, Norms2)
{
  variable i, ncomps;

  ncomps = length(Temps);

  aped_multi_rv_init( ncomps );

  fit_fun("aped_multi_rv(1)+aped_multi_rv(2)");

  set_par( 1, Vrads[0], 0, Vrads[0]-600, Vrads[0]+600);

  for ( i = 0; i < ncomps; i++)
  {
    set_par( i+2, Norms1[i], 0, 0, 100*Norms1[i] );
    set_par( i+ncomps+2, Temps[i], 0, 1.e4, 7.95e8 );
  }


  set_par( ncomps*2+2, Vrads[1], 0, Vrads[1]-600, Vrads[1]+600);

  for ( i = 0; i < ncomps; i++)
  {
    set_par( i+2+2*ncomps+1, Norms2[i], 0, 0, 100*Norms2[i] );
    set_par( i+ncomps+2+2*ncomps+1, Temps[i], 0, 1.e4, 7.95e8 );
  }

%  list_par;
}  

%-

% added following to this file 2002.05.15...
%
%; File: apec_contin_fit.sl
%; Author: D. Huenemoerder
%; Original version: 2001.10.19
%;====================================================================
%; version 3.0  2002.05.15; revise to use model_spectrum's line or
%;                          contin parameters, and merge into
%;                          aped_fit_models.sl. 
%;
%; version 2.0  2002.03.12; fix regrid -> rebin
% version: 1.0
%%%%%%%%%%%%%%%%%%%%%%%%
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
%% Now we simply regrid the Contin_Value to the specified grid,
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
  variable mdl = default_plasma_state(); 

  mdl.norm = par[0];     
  mdl.temperature = par[1];

  if ( __is_initialized( &Aped_Elem_List ) )   % apply different abundances, if set
  {
      mdl.elem = Aped_Elem_List;
      mdl.elem_abund = Aped_Abund_List;
  }

  define_model(mdl);
  y = model_spectrum (xlo, xhi, MODEL_CONTIN);
    
  return y;

} %- define

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define an aped multi-T continuum function.
% Proper use of this requires ignoring line features.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
define aped_multi_contin_init ( ncomps )
{
% Dynamically define a continuum model w/ ncomps temperature compents.

  variable i, func_name, T_params, Norm_params, params;

  params = "norm_1";

  for (i=2; i <= ncomps; i++)
    params = [params, sprintf("norm_%d",i)];

  for (i=1; i <= ncomps; i++)
    params = [params, sprintf("temp_%d", i)];

  add_slang_function( "aped_multi_contin", params, [1:ncomps] );

}
%-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+

define aped_multi_contin_setup(Temps, Norms)
{
  variable i;

  aped_multi_contin_init( length(Temps) );

  fit_fun("aped_multi_contin(1)");

  for ( i = 0; i < length(Temps); i++)
  {
    set_par( i+1, Norms[i], 0, 0, 100*Norms[i] );
    set_par( i+length(Temps)+1, Temps[i], 0, 1.e4, 7.95e8 );
  }

%  list_par;
}  

%-

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+
define aped_multi_contin_fit ( xlo, xhi, par )
{
%
% PURPOSE:  isis user model, APED thermal plasma continuum, multi-temperature,
%	    cosmic abundances, low density limit.

% xlo, xhi = grid of wavelength bins;
% 
  variable ii,                        	     % loop counter
	   y,                         	     % return value
           ncomps = length(par)/2,    	     % # model components
	   pstate = default_plasma_state();  % plasma model definition

  variable mdl = Struct_Type[ncomps];
  variable norms=par[ [0:ncomps-1] ];
  variable temps=par[ [ncomps:2*ncomps-1] ];	   	   

  for (ii = 0; ii < ncomps; ii++ )
  {
    mdl[ii] =             @pstate;
    mdl[ii].temperature = temps[ii];
    mdl[ii].norm =        norms[ii];

    if ( __is_initialized( &Aped_Elem_List ) )   % apply different abundances, if set
    {
        mdl[ii].elem = Aped_Elem_List;
        mdl[ii].elem_abund = Aped_Abund_List;
    }

  }

  define_model(mdl);                      % define the model (use list_model; to see).

  y = model_spectrum (xlo, xhi, MODEL_CONTIN);

  return y;
}  
%-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


add_slang_function("aped_contin", ["norm","temp"]);
add_slang_function ("fast_aped_contin", ["norm"]);
