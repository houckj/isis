% Time-stamp: <2001-10-16 14:05:15 dph> 
% MIT Directory: ~dph/libisis
% CfA Directory: ~dph/libisis
% File: fspec_summary.sl
% Author: D. Huenemoerder
% Original version: 2001.10.12
%====================================================================

%_debug_info=1;             % turn these on to debug.
%_traceback=1;

% COMMENTARY

% PURPOSE
%  Make a summary plot of a spectrum from a pha II file.
%  Plots photon flux vs wavelength to postscript file.
%  If a model was specified, mark line positions.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROVIDES:
% 
% hfspec_summary() : make flux summary spec plot;  +-1, rebin MEG and/or HEG
% lfspec_summary() : make flux summary spec plot;  +-1, rebin LEG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE:

%  (mspec,mmspec, hspec,mhspec)=hfspec_summary( pha_file, order, outfile, harf, marf[, binfact[, mdl_temps[, mdl_wgts]]]);");
%    pha_file = file name"
%    order     = diffraction order (tg_m); integer; e.g., -1
%    outfile   = filename for output ps file.
%    binfact   = optional rebinning factor, pixels
%            (default=0 => no rebinning; HEG is binned by same amount.)
%    harf  = HEG arf filename corresponding to the order.
%    marf  = MEG arf filename corresponding to the order.
%    mdl_temps = optional log10 temperatures of model components (for line labeling)
%    mdl_wgts =  optional weights of model components.
%  Returns meg_flux, meg_model, heg_flux, heg_model

% (lspec,mlspec)=lfspec_summary( pha_file, order, outfile, arf[, binfact[, mdl_temps[, mdl_wgts]]]);
%    pha_file =  file name
%    order     = diffraction order (tg_m); integer; only -1 or +1 for LETG/HRC-S
%    outfile   = filename for output ps file.
%    arf       = filename for arf corresponding to the order.
%    binfact   = optional rebinning factor, pixels
%                (default=0 => no rebinning; HEG is binned by same amount.)
%    mdl_temps = optional log10 temperatures of model components (for line labeling).
%    mdl_wgts =  optional weights of model components.
%  Returns leg_flux, leg_model, writes a postscript file.



% EXAMPLES:

%
% LETG spectrum summary; since it's not an emission line source, no model given.
%
% isis> (sl,ml)=lfspec_summary("Mrk_421/Pha/hrcf01715N001_pha2.fits",
% isis> -1, "Mrk_421_letgs.ps", Arf/hrcf01715_002N001LEG_-1_garf.fits,12);
% isis> ! gv -noswap -noantialias Tst2_leg.ps  # inspect ps file
%
% Plot the result in isis:
% isis> xrange(1.5,20); yrange(0); hplot(sl.bin_lo, sl.bin_hi, sl.value);

% +1 for HETG/ACIS-S observation, specify multi-temperature model.
%
% isis> fpha1 = "./AR_Lac/Pha/acism00006_004N000_pha2.fits"; % file name
% isis> temps= [6.5,6.8,7.0,7.5,7.8,8.0];    % 6 Temperature model
% isis> wgts= [1.0,2.0,10.0,10.0,5.0,1.0];   % weights 10.e-14 n_e n_H VEM/4piD^2
% isis> harf= "acism00006_004N000HEG_1_garf.fits"
% isis> marf= "acism00006_004N000MEG_1_garf.fits"
%
% isis> plasma(aped);
% isis> (cm,mm,ch,mh)=hfspec_summary(fpha, 1, "ARLac_hetgs_summary.ps",harf,marf, 2, temps,wgts);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%


% Subsidiary functions:

%
%  pick_el()             Given array of line indices return subset
%  p_groups()            Line label plotting utility
%   

%
%  NOTE: lfspec_summary() might work for both LETG/HRC-S and
%        LETG/ACIS-S.  Hard-wired max wavelength limit might cause a burp.

%
% RESTRICTIONS:  Will not work for pha Type I files, because in this
%		 case, ISIS assumes that there is no wavelength array.
%
%                If you specify a model, the plasma database must
%                be loaded (e.g., isis> plasma(aped); )
%                The reason it isn't done by this program is because
%                it takes many seconds to load, and would be a
%                hindrance to interactive use.
%
%                ISIS 0.9.50 or later for default_plasma_state(),
%                define_model() functions.


% SIDE AFFECTS: Leaves the plot window unselected.
%               Scribbles many "failed retrieving info for line"
%               messages if a model is specified. These are harmless.
%               They occur when some data, which isn't needed here, is
%               missing from the plasma database (such as
%               spectroscopic line labels).
%
%               Writes a postscript file in the current directory.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%+ utility functions :           

%%%%%%%%%%%%%%%%%%% pick ion from an index array %%%%%%%%%%%%%%%%%%%%%%
%+  
define pick_el()
  {
  % given array of line indices, a, return the subset of the array
  % which are lines of element, el, and ion, ion.
  % example:
  %     a = brightest(30, where(wl(10,12))); % brightest in a region
  %     a_Fe = pick_el(a,Fe);                % pick the iron lines in the list
  % If there aren't any, returns the array, [-1]

  variable a, ion, el;

  switch (_NARGS)

  {
  case 3:
  ion = ();
  el = ();
  a = ();
  }

  {
  case 2:
  el = ();
  a = ();
  ion = NULL;
  }

  {
  message("USAGE: g = pick_el(array, elem[, ion]);");
  }

  variable ii, s, n = length(a);
  variable result = -1;         % dummy value

  for (ii=0; ii<n; ii++)           % each line in array
    {
    s = line_info(a[ii]);

    switch (ion)

     {
     case NULL:
     if (s.Z == el)
       result = [result, a[ii]];
     }

     {
     % default case:
     if ( (s.Z == el) and (s.ion == ion) )
       result = [result, a[ii]];
     }

    } % end for each line.

  % were any lines found?
  if (length(result) > 1)
    {
    result = result[ [1:length(result)-1] ];  % truncate the -1
    }
  return result;
  }
%-  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% +   line label plotting utility; plot element groups in region;
%                                  colorize, offset
define p_groups (xlo, xhi)
  {
  variable ii, nbright=50;

  % ISIS pre-defined colors:
  % red=2, green=3, blue=4, purple=6, yellow=7, orange=8, grey=15
  % Define some additional colors (done by experiment) :

  variable light_blue=5, bluegreen=10, periwinkle=11,
           violet=12, magenta=13, slategray=14;
  variable p_elems = [Fe, Ar, O, Si, Mg, Ne, S, Ca, Ni, Al];
  variable p_color =
       [green, light_blue, magenta, orange, purple, red, bluegreen,
           periwinkle, violet, grey];
  variable p_top = [0.8, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.425, 0.40];
  variable p_bot = p_top - 0.1;

  variable grp = brightest(nbright, where( wl(xlo,xhi) ) );
  variable s = line_label_default_style();
  variable egrp;

  for (ii=0; ii<length(p_elems); ii++)
    {
     egrp = pick_el(grp, p_elems[ii]);
     if (egrp[0] != -1)
      {
      s.top_frac = p_top[ii];
      s.bottom_frac = p_bot[ii];
      plot_group([egrp], p_color[ii], s);
      }
    } % endfor p_elems.

   % label any other elements, ions in the region:
   plot_group(brightest(5, where( wl(xlo,xhi)  and not (el_ion(p_elems)))), slategray);
  }
%-


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+     MAIN functions        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
define hfspec_summary( )
{ %+define
  variable pha_file, order, outfile, harf, marf, binfact, mdl_temps, mdl_wgts;
  variable ii ;                         % loop counter
  variable hists ;                      % histogram indices array
  variable d;                           % temp data struct
  variable HEG=1, MEG=2;                % TG_PART mapping
  variable HAVE_MODEL=0;

  if (_NARGS < 5)
    {
    message("");
    message("Purpose:  make flux summary spec plot; rebinned MEG and HEG");
    message("USAGE:  (mspec,mmspec, hspec, mhspec, mhspec)=hfspec_summary( pha_file, order, outfile, harf, marf[, binfact[, mdl_temps[, mdl_wgts]]]);");
    message("    pha_file = file name");
    message("    order     = diffraction order (tg_m); integer; e.g., -1");
    message("    outfile   = filename for output ps file.");
    message("    harf  = HEG arf filename corresponding to the order.");
    message("    marf  = MEG arf filename corresponding to the order.");
    message("    binfact   = optional rebinning factor, pixels;");
    message("                (default=0 => no rebinning; HEG is binned by same amount.)");  
    message("    mdl_temps = optional log10 temperatures of model components (for line labeling; defaul=6.8)");
    message("    mdl_wgts =  optional weights of model components (default=1).");
    message("  Returns meg_flux, meg_model, heg_flux, heg_model");
    return 1;
    }
 
  binfact = 1;    % default - don't rebin.

%+ pop arguments

   switch(_NARGS)

    {
    case 8:
    mdl_wgts = ();
    mdl_temps = ();
    binfact = ();
    marf = ();
    harf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=1;
    }
     
    {
    case 7:
    mdl_temps = ();
    mdl_wgts = mdl_temps-mdl_temps + 1.0;
    binfact = ();
    marf = ();
    harf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=1;
    }
     
    {
    case 6:
    binfact = ();
    marf = ();
    harf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=0;
    }

    {
    case 5:
    marf = ();
    harf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=0;
    }
%-


%+ load the data

  hists = load_data(pha_file);
  variable n_harf = load_arf(harf);
  variable n_marf = load_arf(marf);
  variable hists_info = get_data_info(hists);     % get info on histograms.
  
%-

%+ find order; rebin

  variable h_idx;
  variable xm1, xm2, cm, HAVE_MEG=0;   % histogram bin coordinates and counts for MEG...
  variable xh1, xh2, ch, HAVE_HEG=0;   % ... and HEG.

  plot_bin_density;

  % find the HEG order from the histogram lists:
  %
  h_idx = where( (hists_info.part == HEG) and (hists_info.order == order) );
  h_idx=h_idx[0] ;   % better be scalar!
  if (length(h_idx) > 0)
  {
    HAVE_HEG = 1;
    variable heg_idx = hists[h_idx];
    assign_arf(n_harf, heg_idx);
    group_data(heg_idx, binfact);     % rebin
    flux_corr(heg_idx);
    variable hfspec = get_data_flux(heg_idx);
    xh1 = hfspec.bin_lo;
    xh2 = hfspec.bin_hi;
    ch = hfspec.value / (xh2-xh1);
  }    


  % find the MEG 1st orders from the histogram lists:
  %
  h_idx = where( (hists_info.part == MEG) and (hists_info.order == order) );
  h_idx=h_idx[0] ;   % better be scalar!
  if (length(h_idx) > 0)
  {
    HAVE_MEG = 1;
    variable meg_idx = hists[h_idx];
    assign_arf(n_marf, meg_idx);
    group_data(meg_idx, binfact);    % rebin
    flux_corr(meg_idx);
    variable mfspec = get_data_flux(meg_idx);
    xm1 = mfspec.bin_lo;
    xm2 = mfspec.bin_hi;
    cm = mfspec.value / (xm2-xm1);
  }

%-

  if (_isis_version < 950)
  {
    message("ISIS version does not support type of model. Continuing without.");
    HAVE_MODEL=0;      % unset it.
  }
  if (HAVE_MODEL)
  { %+ have model
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %+    Construct model
  % the following is for isis version 0.9.50 and higher only.
  %  for prev, need to load model from a table.

    mdl_temps = 10.0^mdl_temps;
    variable mdl = Struct_Type[length(mdl_temps)];
    variable pstate = default_plasma_state();
  
  %  plasma(aped);   % do externally? (it takes time to load)
  
    for (ii=0; ii<length(mdl_temps); ii++)
    {
      mdl[ii] = @pstate;
      mdl[ii].temperature = mdl_temps[ii];
      mdl[ii].norm = mdl_wgts[ii];
    }
	
    define_model(mdl);
    list_model;     % optional; for info

  %- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %+ evaluate model on grids:
    if (HAVE_HEG)
    {
	vmessage("Evaluating model on HEG grid........");     %%%%%%%% DEBUG
	variable fh = model_spectrum(xh1, xh2);   % evaluate the model for HEG...
    }
  
    if (HAVE_MEG)
    {
	vmessage("Evaluating model on MEG grid........");     %%%%%%%% DEBUG
	variable fm = model_spectrum(xm1, xm2);   %  and MEG
			% (do MEG last for "brightest" to work
			%     over full wavelength range)
    }
  } %- have model
%-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%+ set limits, make plot;

  variable pid;
  variable xmin, xmax, ymin,ymax, xlo, xhi, dx;

  pid = open_plot(outfile+"/cps");      % color landscape postscript

  resize(10*2.54, 7./10.);              % 10 by 7 inches
  xmin=1.5;  xmax= 25.5; dx = 4.0;      % full range, and range per plot
  if (HAVE_HEG and  not HAVE_MEG) xmax=17.;

  % get info about the bin sizes to put in the label:
  variable hbinsize="", mbinsize="", sbinsize="";
  if (HAVE_HEG)
  {
    hbinsize = string(xh1[1]-xh1[0]) + "H";
  }
  if (HAVE_MEG)
  {
    mbinsize = string(xm1[1]-xm1[0]) + "M";
  }

  sbinsize = "(" + mbinsize + "/" +  hbinsize + "A)" ; 

  for (ii=0; ii<(xmax-xmin)/dx; ii++)   % make each page
  {
    %+ Plot the data:

    label("Wavelength [Angstroms]",
        "Flux [photons/cm^2/s/A] " + sbinsize,
	pha_file);

    xlo = xmin + ii*dx;
    xhi = xlo + dx;
    xrange(xlo - 0.25, xhi);  ;  % minus a bit from min so plots overlap
    ymin=0;

    if (HAVE_MEG)
      ymax = max( cm[where( (xm1 > xlo) and (xm2 < xhi))]) * 1.5;
    else if (HAVE_HEG)
      ymax = max( ch[where( (xh1 > xlo) and (xh2 < xhi))]) * 1.5;
    
    yrange(ymin,ymax);

    if (HAVE_MEG and HAVE_HEG)
    {
      plot_data_flux(heg_idx,2);
      oplot_data_flux(meg_idx,3);
    }
    else if (HAVE_MEG)
    {
      plot_data_flux(meg_idx);
    }
    else if (HAVE_HEG)
    {
      plot_data_flux(heg_idx);
    }
    
    if ( (HAVE_HEG or HAVE_MEG) and HAVE_MODEL) p_groups(xlo, xhi);
%-

    if (HAVE_MODEL)
    {
    %+ make plot of the model:
    % 
      label("Wavelength [Angstroms]",
	  "Flux [photons/cm^2/s/A] " + sbinsize,
	  "Model for " + pha_file);
  
      if (HAVE_MEG)
	ymax = max( fm[where( (xm1 > xlo) and (xm2 < xhi))]) * 1.5 /(xm2[0]-xm1[0]);
      else if (HAVE_HEG)
	ymax = max( fh[where( (xh1 > xlo) and (xh2 < xhi))]) * 1.5/(xh2[0]-xh1[0]);
  
      yrange(ymin,ymax);
  
      if (HAVE_MEG and HAVE_HEG)
      {
	hplot(xh1,xh2,fh,4);
	ohplot(xm1,xm2,fm,5);
      }
      else if (HAVE_MEG)
      {
	hplot(xm1,xm2,fm,2);
      }
      else if (HAVE_HEG)
      {
	hplot(xh1,xh2,2);
      }
      
      if (HAVE_HEG or HAVE_MEG)  p_groups(xlo, xhi);
    }
  %-
  }     %- Plot the data:

  close_plot(pid);
  

%+ cleanup this histogram arrays.
%   (this is important for interactive sessions, to leave the histograms
%          as they were when this started)

  delete_data(hists);
  delete_arf(n_marf);
  delete_arf(n_harf);
%
%-

%+ return structures w/ flux and model spectra:

%             flux          model
%           MEG    HEG    MEG     HEG
  variable mspec, hspec, mmspec, mhspec;

  mspec = struct                % define structure
    {
    bin_lo, bin_hi, value
    };

  mmspec=@mspec;                 % initialize
  hspec=@mspec;
  mhspec=@mspec;

  if (HAVE_MEG)
  {
    mspec.bin_lo = xm1;           % meg flux
    mspec.bin_hi = xm2;
    mspec.value = cm;

    if (HAVE_MODEL)
    {
      mmspec.bin_lo = xm1;
      mmspec.bin_hi = xm2;
      mmspec.value = fm;
    }
  }

  if (HAVE_HEG)
  {
    hspec.bin_lo = xh1;
    hspec.bin_hi = xh2;
    hspec.value = ch;

    if (HAVE_MODEL)
    {
      mhspec.bin_lo = xh1;
      mhspec.bin_hi = xh2;
      mhspec.value = fh;
    }
  }

  return mspec,mmspec, hspec, mhspec;
    
%-

}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
define lfspec_summary( )
{ %+define
  variable pha_file, order, outfile, larf, binfact, mdl_temps, mdl_wgts;
  variable ii ;                         % loop counter
  variable hists ;                      % histogram indices array
  variable d;                           % temp data struct
  variable LEG=3;                % TG_PART mapping
  variable HAVE_MODEL=0;

  if (_NARGS < 4)
    {
    message("");
    message("Purpose:  make flux summary spec plot; rebinned LEG");
    message("USAGE:  (mspec,mmspec, hspec, mhspec, mhspec)=hfspec_summary( pha_file, order, outfile, harf, marf[, binfact[, mdl_temps[, mdl_wgts]]]);");
    message("    pha_file = file name");
    message("    order     = diffraction order (tg_m); integer;  -1 or +1 for LEG");
    message("    outfile   = filename for output ps file.");
    message("    larf  = LEG arf filename corresponding to the order.");
    message("    binfact   = optional rebinning factor, pixels;");
    message("                (default=0 => no rebinning; HEG is binned by same amount.)");  
    message("    mdl_temps = optional log10 temperatures of model components (for line labeling; defaul=6.8)");
    message("    mdl_wgts =  optional weights of model components (default=1).");
    message("  Returns leg_flux, leg_model");
    return 1;
    }
 
  binfact = 1;    % default - don't rebin.

%+ pop arguments

   switch(_NARGS)

    {
    case 7:
    mdl_wgts = ();
    mdl_temps = ();
    binfact = ();
    larf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=1;
    }
     
    {
    case 6:
    mdl_temps = ();
    mdl_wgts = mdl_temps-mdl_temps + 1.0;
    binfact = ();
    larf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=1;
    }
     
    {
    case 5:
    binfact = ();
    larf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=0;
    }

    {
    case 4:
    larf = ();
    outfile = ();
    order = ();
    pha_file = ();
    HAVE_MODEL=0;
    }
%-


%+ load the data

  hists = load_data(pha_file);
  variable n_larf = load_arf(larf);
  variable hists_info = get_data_info(hists);     % get info on histograms.
  
%-

%+ find order; rebin

  variable h_idx;
  variable x1, x2, c, HAVE_LEG=0;   % histogram bin coordinates and counts for LEG...

  plot_bin_density;

  % find the LEG order from the histogram lists:
  %
  h_idx = where( (hists_info.part == LEG) and (hists_info.order == order) );
  h_idx=h_idx[0] ;   % better be scalar!
  if (length(h_idx) > 0)
  {
    HAVE_LEG = 1;
    variable leg_idx = hists[h_idx];
    assign_arf(n_larf, leg_idx);
    group_data(leg_idx, binfact);     % rebin
    flux_corr(leg_idx);
    variable lfspec = get_data_flux(leg_idx);
    x1 = lfspec.bin_lo;
    x2 = lfspec.bin_hi;
    c = lfspec.value / (x2-x1);
  }    

%-

  if (_isis_version < 950)
  {
    message("ISIS version does not support type of model. Continuing without.");
    HAVE_MODEL=0;      % unset it.
  }
  if (HAVE_MODEL)
  { %+ have model
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %+    Construct model
  % the following is for isis version 0.9.50 and higher only.
  %  for prev, need to load model from a table.

    mdl_temps = 10.0^mdl_temps;
    variable mdl = Struct_Type[length(mdl_temps)];
    variable pstate = default_plasma_state();
  
  %  plasma(aped);   % do externally? (it takes time to load)
  
    for (ii=0; ii<length(mdl_temps); ii++)
    {
      mdl[ii] = @pstate;
      mdl[ii].temperature = mdl_temps[ii];
      mdl[ii].norm = mdl_wgts[ii];
    }
	
    define_model(mdl);
    list_model;     % optional; for info

  %- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %+ evaluate model on grids:
    if (HAVE_LEG)
    {
	vmessage("Evaluating model on LEG grid........");     %%%%%%%% DEBUG
	variable f = model_spectrum(x1, x2);   % evaluate the model for LEG...
    }

  } %- have model
%-
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%+ set limits, make plot;

  variable pid;
  variable xmin, xmax, ymin,ymax, xlo, xhi, dx;

  pid = open_plot(outfile+"/cps");      % color landscape postscript

  resize(10*2.54, 7./10.);              % 10 by 7 inches
  xmin=1.5;  xmax= 180; dx = 20.0;      % full range, and range per plot

  % get info about the bin sizes to put in the label:
  variable sbinsize="";
  if (HAVE_LEG)
    sbinsize = string(x1[1]-x1[0]);

  for (ii=0; ii<(xmax-xmin)/dx; ii++)   % make each page
  {    %+ Plot the data:

    label("Wavelength [Angstroms]",
        "Flux [photons/cm^2/s/A] " + sbinsize,
	pha_file);

    xlo = xmin + ii*dx;
    xhi = xlo + dx;
    xrange(xlo - 0.25, xhi);  ;  % minus a bit from min so plots overlap
    ymin=0;

    if (HAVE_LEG)
    {
      ymax = max( c[where( (x1 > xlo) and (x2 < xhi))]) * 1.5;
      if (ymax > 0) yrange(ymin,ymax);
      plot_data_flux(leg_idx,2);
      if (HAVE_MODEL) p_groups(xlo, xhi);
    }
%-

    if (HAVE_MODEL and HAVE_LEG)
    {
    %+ make plot of the model:
    % 
      label("Wavelength [Angstroms]",
	  "Flux [photons/cm^2/s/A] " + sbinsize,
	  "Model for " + pha_file);
  
      ymax = max( f[where( (x1 > xlo) and (x2 < xhi))]) * 1.5 /(x2[0]-x1[0]);
      yrange(ymin,ymax);
      hplot(x1,x2,f,4);
      p_groups(xlo, xhi);
    }
  %-
  }     %- Plot the data:

  close_plot(pid);
  

%+ cleanup this histogram arrays.
%   (this is important for interactive sessions, to leave the histograms
%          as they were when this started)

  delete_data(hists);
  delete_arf(n_larf);
%
%-

%+ return structures w/ flux and model spectra:

%             flux          model
%           MEG    HEG    MEG     HEG
  variable lspec, mlspec;

  lspec = struct                % define structure
    {
    bin_lo, bin_hi, value
    };

  mlspec=@lspec;                 % initialize

  if (HAVE_LEG)
  {
    lspec.bin_lo = x1;           % leg flux
    lspec.bin_hi = x2;
    lspec.value = c;

    if (HAVE_MODEL)
    {
      mlspec.bin_lo = x1;
      mlspec.bin_hi = x2;
      mlspec.value = f;
    }
  }


  return lspec,mlspec;
}    
%-

