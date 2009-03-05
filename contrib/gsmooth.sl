% Time-stamp: <2002-05-20 23:23:37 dph> 
% MIT Directory: ~dph/libisis
% File: gsmooth.sl
% Author: D. Huenemoerder
% Original version: 2000.10.29
%====================================================================
% 2000.10.29 dph
% v.1
% 2001.07.07
% v.2      ------ apply hanning filter; re-align on input x array.
%        NOTE: v.2 changed interface!
% v.2.1  fixed problem w/ shifting - had freq array wrong...
% v.2.2  fixed bug w/ odd-length arrays.


% gaussian smoothing
%
% convolve y by gaussian kernel, width sigma, in x's units
% Use fft1d
% Assumes uniform grid.
% Does not worry about integral form of y.  Should it?

% EXAMPLE use:
% > ()=evalfile("Sl/hanning.sl");
% > ()=evalfile("Sl/gsmooth.sl");
% >  sigma=0.004;    % approx HEG resolution
% >  % get your data into arrays...  xlo, xhi, counts
% >  scounts=gsmooth(xlo, xhi, counts, sigma);
% >  plot(xlo,xhi,counts,4); ohplot(xlo, xhi, scounts,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example test case.
% x1 must be defined.  e.g., [0:8191:0.005]+1.  % MEG xlo.
%
%   evalfile(ilib+"/hanning.sl");
%   evalfile(ilib+"/gsmooth.sl");
%   d=(sin(x1*PI*2/0.5)+5.) *exp(-x1/10);
%   d[[4000:5000:1]]=2.0;
%   d[4500]=10.;
%   sd=gsmooth(x1,x2,d,0.008);
%   xrange(0,45);hplot(x1,x2,d,4);ohplot(sx1,sx2,sd,3);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

require("hanning");

define gsmooth()
{
  variable xlo, xhi, y, sigma;

  variable k, karray, dx, nsig;
  variable nkbins;
  variable rffty, iffty, rfftk, ifftk, cffty, cfftk, c, sy, isy, nrm;
  variable xlo_out, xhi_out;
  variable ymean;

  if (_NARGS != 4)
    {
    message("");
    message("USAGE: smoothed_y = gsmooth(x_lo, x_hi, y, sigma)");
    message("");
    message("Apply gaussian smoothing filter of width sigma to histogram array y.");
    message("  x_lo, x_hi are histogram bin edge coordinate arrays.");
    message("  y is histogram values array (x_lo,x_hi,y must have same lengths");
    message("METHOD: Convolution is performed via FFT, after applying Hanning filter.");
    message("RESTRICTIONS: grid must be uniform.");
    return 1;
    }

%+ pop arguments:
    sigma=();
    y=();
    xhi=();
    xlo=();
%-  

  % make kernel
  dx =  xhi[0]-xlo[0];
  nsig = sigma/dx;          % # bins in sigma.
  nkbins = int(4.*nsig)+1;  % # bins to use in the gaussian kernel; +-2sigma.
  k=exp( -([0:nkbins-1:1]-nkbins/2.)^2/ nsig^2/2.) / (sqrt(2.*PI)*nsig); % eval kernel

%+
  % if arrays have odd # points, make even; discard at end:
  variable is_odd = 0;
  variable n = length(y);
  if (n mod 2) is_odd = 1;
  if (is_odd) 
  {
    y = [y,0];
    xlo = [xlo, xlo[n-1]+dx ];
    xhi = [xhi, xhi[n-1]+dx ];
  }
%-

% make array of kernel, same length as y:
%
  karray =  y*0.0;
  karray[ [0:nkbins-1:1] ] = k;

%%+ apply Hanning window (cosine)
   variable f;
   f = hanning( length(y) );
   y = y*f;
%%-

%%+ forward fft
%
  (rffty, iffty) = fft1d(y, y*0., -1);               % fft of data
  (rfftk, ifftk) = fft1d(karray, karray*0., -1);     % fft of kernel

  cffty = rffty + iffty * 1i;     % convert to Slang complex types.
  cfftk = rfftk + ifftk * 1i;

  c = cffty * cfftk ;	                    % convolve: fft product.
%%-

%%+ apply phase shift of half-kernel width in order to return
%
% smoothed array on same x array as input:
% shift by nkbins/2
%
  variable phase_shift, freqs;

  % fft defn of frequency array:   1/N == 1 cycle per N bins.
  %                               (N/2)/N == Nyquist;
  %
  % [DC, 1/N, 2/N..., +-(N/2)/N, ... -3/N, -2/N, -1/N]

  freqs = [0:length(y)/2:1];   % DC to Nyquist
  freqs = 2.0*PI/length(y) * [freqs, -[length(y)/2-1:1:-1]];   

  phase_shift = freqs*(nkbins/2.);
  phase_shift = cos(phase_shift) + sin(phase_shift)*1i;
  c = c*phase_shift;
%%-  

%+ reverse fft
  (sy, isy) = fft1d(Real(c), Imag(c), 1);   % inverse fft;
%%-


%
  variable lz;   % indices of where filter is != 0.
  lz = where( f != 0.0);
  nrm = sum(y) / sum(sy[lz]); % ad hoc normalization; should be known a priori
  sy = sy * nrm;
  sy[lz] = sy[lz] / f[lz] ;   % remove filter.	

  if (is_odd)               % truncate 1 bin padding...
    sy = sy[[0:n-1]];

  return sy;

}
provide("gsmooth");
