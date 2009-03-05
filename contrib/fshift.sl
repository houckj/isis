% Time-stamp: <2001-10-17 11:35:39 dph> 
% MIT Directory: ~dph/libisis
% File: fshift.sl
% Author: D. Huenemoerder
% Original version: 2001.07.09
%====================================================================
% v.1

% shift an array by fractional bin amount, using phase shift in
% Fourier space.
%
% Use fft1d
% Assumes uniform grid.
% Does not worry about integral form of y.  Should it?

% NOTES:  
%
%  Empirical tests seem to indicate that the Hanning filter isn't
%  needed.
%  Tests also show that a high-frequency cutoff filter is probably
%  needed. A narrow square feature rings for fractional bin shifts.  A
%  narrow gaussian is OK.



% EXAMPLE use:
% > ()=evalfile("Sl/gsmooth.sl");
% > ()=evalfile("Sl/hanning.sl");
% >  dpix=22.3;    
% >  % get your data into arrays...  xlo, xhi, counts
% >  scounts = fshift(xlo,xhi,counts, dx);
% >  plot(xlo,xhi,counts,4); ohplot(xlo, xhi, scounts,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example test case.
% x1 must be defined.  e.g., [0:8191:0.005]+1.  % MEG xlo.
%
%   evalfile(ilib+"/fshift.sl");
%   d=(sin(x1*PI*2/0.5)+5.) *exp(-x1/10);
%   d[[4000:5000:1]]=2.0;
%   d[[4500:4510:1]]=10.;
%   d[[4100:4120:1]]= exp(-0.5*([4100:4120:1]-4110.)^2/2.^2)*10. + 2.;
%   sd=fshift(x1,x2,d,1.3);
%   xrange(0,45);hplot(x1,x2,d,4);ohplot(x1,x2,sd,3);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

require("hanning");   % window function

define fshift()
{
  variable xlo, xhi, y, dx, npix;
  variable rffty, iffty, cffty, sy, isy, nrm;
  variable rfftf, ifftf, cfftf, sf, isf;

  if (_NARGS != 4)
    {
    message("");
    message("USAGE: shifted_y = fshift(x1, x2, y, dx);");
    message("");
    message("Shift bins of y by amount dx.");
    message("METHOD: phase shift FFT, after applying Hanning filter.");
    message("RESTRICTIONS: grid must be uniform.");
    return 1;
    }

%+ pop arguments:
    dx=();
    y=();
    xhi=();
    xlo=();
%-  

%%+ apply Hanning window (cosine)
   variable f;
   f = hanning(length(y));
%   y = y*f;
%%-

%+ determine number of pixels to shift:
  npix = dx / (xhi[0]-xlo[0]);
%-  
   

%%+ forward fft
%
  (rffty, iffty) = fft1d(y, y*0., -1);               % fft of data
  cffty = rffty + iffty * 1i;     % convert to Slang complex types.
%%-

%%+ apply phase shift 
  variable phase_shift, freqs;

  % fft defn of frequency array:   1/N == 1 cycle per N bins.
  %                               (N/2)/N == Nyquist;
  %
  % [DC, 1/N, 2/N..., +-(N/2)/N, ... -3/N, -2/N, -1/N]

  freqs = [0:length(y)/2:1];   % DC to Nyquist
  freqs = 2.0*PI/length(y) * [freqs, -[length(y)/2-1:1:-1]];   
  phase_shift = freqs*(-npix);
  phase_shift = cos(phase_shift) + sin(phase_shift)*1i;
  cffty = cffty*phase_shift;
%%-  

%+ reverse fft
  (sy, isy) = fft1d(Real(cffty), Imag(cffty), 1);   % inverse fft;
%%-


% shift filter
%%+ forward fft filter
%
  (rfftf, ifftf) = fft1d(f, f*0., -1);               % fft of data
  cfftf = rfftf + ifftf * 1i;     % convert to Slang complex types.
%%-

%%+ apply phase shift to filter
  cfftf = cfftf*phase_shift;
%%-  

%+ reverse fft filter
  (sf, isf) = fft1d(Real(cfftf), Imag(cfftf), 1);   % inverse fft;
%%-



%
  variable lz;   % indices of where filter is != 0.
  lz = where( sf != 0.0);
  nrm = sum(y) / sum(sy[lz]); % ad hoc normalization; should be known a priori
  sy = sy * nrm;
%  sy[lz] = sy[lz] / sf[lz] ;   % remove filter.	

  return sy;

}

    

  

