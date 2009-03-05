% Time-stamp: <2001-10-16 14:08:07 dph> 
% MIT Directory: ~dph/libisis
% CfA Directory: ~dph/libisis
% File: cspec_summary.sl
% Author: D. Huenemoerder
% Original version: 2001.10.12
%====================================================================

%_debug_info=1;             % turn these on to debug.
%_traceback=1;

% COMMENTARY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROVIDES:
% 
% hcspec_summary() : make counts summary spec plot; sum +-1, rebin MEG and/or HEG
% lcspec_summary() : make counts summary spec plot; sum +-1, rebin LEG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USAGE:

%  (mspec,hspec, mmspec, mhspec)=hcspec_summary( pha_files, outfile[, binfact[, mdl_temps[, mdl_wgts]]]);");
%    pha_files = string array of file names"
%    outfile   = filename for output ps file.
%    binfact   = optional rebinning factor, pixels
%            (default=0 => no rebinning; HEG is binned by same amount.)
%    mdl_temps = optional log10 temperatures of model components (for line labeling)
%    mdl_wgts =  optional weights of model components.
%  Returns meg_counts, meg_model, heg_counts, heg_model

% (lspec,mlspec)=lcspec_summary( pha_files, outfile[, binfact[, mdl_temps[, mdl_wgts]]]);
%    pha_files = string array of file names
%    outfile   = filename for output ps file.
%    binfact   = optional rebinning factor, pixels
%                (default=0 => no rebinning; HEG is binned by same amount.)
%    mdl_temps = optional log10 temperatures of model components (for line labeling).
%    mdl_wgts =  optional weights of model components.
%  Returns leg_counts, leg_model, writes a file.



% EXAMPLES:

%
% LETG spectrum summary; since it's not an emission line source, no model given.
%
% isis> (sl,ml)=lcspec_summary("Mrk_421/Pha/hrcf01715N001_pha2.fits", "Mrk_421_letgs.ps",12);
% isis> ! gv -noswap -noantialias Tst2_leg.ps  # inspect ps file
%
% Plot the result in isis:
% isis> xrange(1.5,20); yrange(0); hplot(sl.bin_lo, sl.bin_hi, sl.value);

% Sum +-1 for 6 HETG/ACIS-S observations, specify multi-temperature
% model.
%
% isis> fpha1 = "./AR_Lac/Pha/acis";       % root of the file names
% isis> fpha3 = "_pha2.fits";              % suffix of the names
% isis> fpha2 = ["f00007_005N001",         % middles of the names
% isis> 	 "f00008_005N001",
% isis> 	 "f00009_003N001",
% isis> 	 "f00010_003N001",
% isis> 	 "f00011_003N001",
% isis> 	 "m00006_004N000"];
% isis> fphas = fpha1 + fpha2 + fpha3;       % concat the strings.
% isis> temps= [6.5,6.8,7.0,7.5,7.8,8.0];    % 6 Temperature model
% isis> wgts= [1.0,2.0,10.0,10.0,5.0,1.0];   % weights 10.e-14 n_e n_H VEM/4piD^2
%
% isis> plasma(aped);
% isis> (cm,mm,ch,mh)=hcspec_summary(fphas, "ARLac_hetgs_summary.ps",2, temps,wgts);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%


% Subsidiary functions:

%
%  spec_sum_add_hists()  Add counts histograms
%  pick_el()             Given array of line indices return subset
%  p_groups()            Line label plotting utility
%   

%
%  NOTE: lcspec_summary() might work for both LETG/HRC-S and
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+ 
define spec_sum_add_hists(spec_nums)
  {
  variable d, c, ii, xlo, xhi;

  d = get_data_counts(spec_nums[0]);
  c = d.value;
  xlo = d.bin_lo;
  xhi = d.bin_hi;
   
  for (ii=1; ii<length(spec_nums); ii++)
    {
    c = c + ( get_data_counts(spec_nums[ii]) ).value;
    }

  return xlo, xhi, c;
  }
%-   

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
define hcspec_summary( )
{ %+define
  variable pha_files, outfile, binfact, mdl_temps, mdl_wgts;
  variable ii ;                         % loop counter
  variable hists ;                      % histogram indices array
  variable d;                           % temp data struct
  variable HEG=1, MEG=2;                % TG_PART mapping
  variable HAVE_MODEL=0;

  if (_NARGS < 2)
    {
    message("");
    message("Purpose:  make counts summary spec plot; sum MEG and rebinned HEG");
    message("USAGE:  (mspec,hspec, mmspec, mhspec)=hcspec_summary( pha_files, outfile[, binfact[, mdl_temps[, mdl_wgts]]]);");
    message("    pha_files = string array of file names");
    message("    outfile   = filename for output ps file.");
    message("    binfact   = optional rebinning factor, pixels;");
    message("                (default=0 => no rebinning; HEG is binned by same amount.)");  
    message("    mdl_temps = optional log10 temperatures of model components (for line labeling; defaul=6.8)");
    message("    mdl_wgts =  optional weights of model components (default=1).");
    message("  Returns meg_counts, heg_counts, meg_model, heg_model");
    return 1;
    }
 
  binfact = 1;    % default - don't rebin.

%  mdl_temps = [6.8];     % default model temperature
%  mdl_wgts  = [1.0];
  
%+ pop arguments

   switch(_NARGS)

    {
    case 5:
    mdl_wgts = ();
    mdl_temps = ();
    binfact = ();
    outfile = ();
    pha_files = ();
    HAVE_MODEL=1;
    }
     
    {
    case 4:
    mdl_temps = ();
    mdl_wgts = mdl_temps-mdl_temps + 1.0;
    binfact = ();
    outfile = ();
    pha_files = ();
    HAVE_MODEL=1;
    }
     
    {
    case 3:
    binfact = ();
    outfile = ();
    pha_files = ();
    HAVE_MODEL=0;
    }

    {
    case 2:
    outfile = ();
    pha_files = ();
    HAVE_MODEL=0;
    }
%-


%+ load the data

  if (length(pha_files) == 1)        % if only one pha file...
   {
   pha_files = [pha_files];          % ...force to array.
   }

  hists = -1;                                % initializa variable to dummy value
  for (ii=0; ii<length(pha_files); ii++)     % load all files
    {
    hists = [hists,load_data(pha_files[ii])];     % concatenate list
    }
  hists = hists[[1:length(hists)-1:1]];	          % truncate dummy element
  variable hists_info = get_data_info(hists);     % get info on histograms.
%-


%+ get and sum counts; rebin

  variable h_idx;
  variable xm1, xm2, cm, HAVE_MEG=0;   % histogram bin coordinates and counts for MEG...
  variable xh1, xh2, ch, HAVE_HEG=0;   % ... and HEG.

  % find the HEG 1st orders from the histogram lists:
  %
  h_idx = where( (hists_info.part == HEG) and (abs(hists_info.order) == 1) );
  if (length(h_idx) > 0)
  {
    HAVE_HEG = 1;
    for (ii=0; ii<length(h_idx); ii++)
      {
      group_data( (hists[h_idx])[ii], binfact);     % rebin
      }
    (xh1,xh2,ch) = spec_sum_add_hists(hists[h_idx]);   % sum them.
  }    


  plot_bin_integral;

  % find the MEG 1st orders from the histogram lists:
  %
  h_idx = where( (hists_info.part == MEG) and (abs(hists_info.order) == 1) );
  if (length(h_idx) > 0)
  {
    HAVE_MEG = 1;
    for (ii=0; ii<length(h_idx);ii++)
      {
      group_data( (hists[h_idx])[ii], binfact);    % rebin
      }
    (xm1,xm2,cm) = spec_sum_add_hists(hists[h_idx]);   % sum
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
        "Flux [Counts/bin] " + sbinsize,
	pha_files[0]+" + ...");

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
      hplot( xm1, xm2, cm, 2);
      ohplot(xh1, xh2, ch, 3);
    }
    else if (HAVE_MEG)
    {
      hplot( xm1, xm2, cm, 2);
    }
    else if (HAVE_HEG)
    {
      hplot( xh1, xh2, ch, 2);
    }
    
    if ( (HAVE_HEG or HAVE_MEG) and HAVE_MODEL) p_groups(xlo, xhi);
%-

    if (HAVE_MODEL)
    {
    %+ make plot of the model:
    % 
      label("Wavelength [Angstroms]",
	  "Flux [photons/cm^2/s/bin] " + sbinsize,
	  "Model for " + pha_files[0]+" + ...");
  
      if (HAVE_MEG)
	ymax = max( fm[where( (xm1 > xlo) and (xm2 < xhi))]) * 1.5;
      else if (HAVE_HEG)
	ymax = max( fh[where( (xh1 > xlo) and (xh2 < xhi))]) * 1.5;
  
      yrange(ymin,ymax);
  
      if (HAVE_MEG and HAVE_HEG)
      {
	hplot( xm1, xm2, fm, 4);
	ohplot(xh1, xh2, fh, 5);
      }
      else if (HAVE_MEG)
      {
	hplot( xm1, xm2, fm, 2);
      }
      else if (HAVE_HEG)
      {
	hplot( xh1, xh2, fh, 2);
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
%
%-

%+ return structures w/ counts and model spectra:

%             counts          model
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
    mspec.bin_lo = xm1;           % meg counts
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%+         LETG CASE
%
define lcspec_summary( )
  {  %+ define
  variable pha_files, outfile, binfact, mdl_temps, mdl_wgts;
  variable ii ; % loop counter
  variable hists ; % histogram indices array
  variable d; % temp data struct
  variable LEG=3; % TG_PART mapping
  variable HAVE_MODEL=0;

  if (_NARGS < 2)
    {
    message("");
    message("Purpose:  make counts summary spec plot; sum MEG and rebinned HEG");
    message("USAGE:  (lspec,mlspec)=lcspec_summary( pha_files, outfile[, binfact[, mdl_temps[, mdl_wgts]]]);");
    message("    pha_files = string array of file names");
    message("    outfile   = filename for output ps file.");
    message("    binfact   = optional rebinning factor, pixels;");
    message("                (default=0 => no rebinning; HEG is binned by same amount.)");  
    message("    mdl_temps = optional log10 temperatures of model components (for line labeling; defaul=6.8)");
    message("    mdl_wgts =  optional weights of model components (default=1).");
    message("  Returns leg_counts, leg_model");
    return 1;
    }
 
  binfact = 1;    % default - don't rebin.

  mdl_temps = [6.8];     % default model temperature
  mdl_wgts  = [1.0];
  
%+ pop arguments

   switch(_NARGS)

    {
    case 5:
    mdl_wgts = ();
    mdl_temps = ();
    binfact = ();
    outfile = ();
    pha_files = ();
    HAVE_MODEL=1;
    }
     
    {
    case 4:
    mdl_temps = ();
    mdl_wgts = mdl_temps-mdl_temps + 1.0;
    binfact = ();
    outfile = ();
    pha_files = ();
    HAVE_MODEL=1;
    }
     
    {
    case 3:
    binfact = ();
    outfile = ();
    pha_files = ();
    HAVE_MODEL=0;
    }

    {
    case 2:
    outfile = ();
    pha_files = ();
    HAVE_MODEL=0;
    }
%-


%+ load the data

  if (length(pha_files) == 1)        % if only one pha file...
   {
   pha_files = [pha_files];          % ...force to array.
   }

  hists = -1;                                % initializa variable to dummy value
  for (ii=0; ii<length(pha_files); ii++)     % load all files
    {
    hists = [hists,load_data(pha_files[ii])];     % concatenate list
    }
  hists = hists[[1:length(hists)-1:1]];	          % truncate dummy element
  variable hists_info = get_data_info(hists);     % get info on histograms.
%-

%+ get and sum counts; rebin

  variable h_idx;
  variable x1, x2, c, HAVE_LEG=0;   % histogram bin coordinates and counts for MEG...

  % find the LEG "1st" orders from the histogram lists:
  %
  h_idx = where( (hists_info.part == LEG) and (abs(hists_info.order) == 1) );
  if (length(h_idx) > 0)
  {
    for (ii=0; ii<length(h_idx); ii++)
      {
      group_data( (hists[h_idx])[ii], binfact);     % rebin
      }
    (x1,x2,c) = spec_sum_add_hists(hists[h_idx]);   % sum them.
    HAVE_LEG = 1;
  }    
%-

  if (HAVE_MODEL)
  {
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
      variable f = model_spectrum(x1, x2);   % evaluate the model for HEG...
    }
  }

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
  {
    sbinsize = string(x1[1]-x1[0]);
  }
  sbinsize = "(" + sbinsize + "A)" ; 

  for (ii=0; ii<(xmax-xmin)/dx; ii++)   % make each page
    {
    %+ Plot the data:

    label("Wavelength [Angstroms]",
        "Flux [Counts/bin] " + sbinsize,
	pha_files[0]+" + ...");

    xlo = xmin + ii*dx;
    xhi = xlo + dx;
    xrange(xlo - 0.25, xhi);  ;  % minus a bit from min so plots overlap
    ymin=0; yrange(ymin);

    if (HAVE_LEG)
    {
      ymax = max( c[where( (x1 > xlo) and (x2 <= xhi))]) * 1.5;
      if (ymax > 0) yrange(ymin,ymax);
      hplot( x1, x2, c, 2);
      if (HAVE_MODEL) p_groups(xlo, xhi);
    }
%-

    if (HAVE_MODEL)
    {
    %+ make plot of the model:
    %   
      label("Wavelength [Angstroms]",
	  "Flux [photons/cm^2/s/bin] " + sbinsize,
	  "Model for " + pha_files[0]+" + ...");
  
      if (HAVE_LEG)
      {
	ymax = max( f[where( (x1 > xlo) and (x2 < xhi))]) * 1.5;
	yrange(ymin,ymax);
	hplot( x1, x2, f, 4);
	p_groups(xlo, xhi);
      }
    }
  %-
    } %- Plot the data

  close_plot(pid);


%+ cleanup this histogram arrays.
%   (this is important for interactive sessions, to leave the histograms
%          as they were when this started)

  delete_data(hists);
%
%-

%+ return structures w/ counts and model spectra:

  variable lspec, mlspec;

  lspec = struct
    {
    bin_lo, bin_hi, value
    };

  mlspec = @lspec;

  if (HAVE_LEG)           % leg counts
  {
    lspec.bin_lo = x1;
    lspec.bin_hi = x2;
    lspec.value = c;

    if (HAVE_MODEL)                % leg model
    {
      mlspec.bin_lo = x1;
      mlspec.bin_hi = x2;
      mlspec.value = f;
    }
  }

  return lspec,mlspec;

%-

}  %-define


