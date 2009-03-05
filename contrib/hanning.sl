% Time-stamp: <2001-10-17 11:36:49 dph> 
% MIT Directory: ~dph/libisis
% File: hanning.sl
% Author: D. Huenemoerder
% Original version: 2000.10.29
%====================================================================
%
% the Hanning window is:  1/2 * (1 - cos(2*PI*i/(N-1))), i=0..N-1
%
define hanning(length)
{
 variable result;

 result = 0.5 * (1.0 - cos( 2*PI/(length-1) * [0:length-1:1]) );

 return result;

} 

% modified hanning - ramp up to 1.0 quicker.
% n_end is number <= length/2.
% length is length of array
%
define mhanning(n_end, length)
{
  variable f1, result;

  f1=hanning(n_end*2);
  result = [0:length-1:1]*0.0 + 1.0;
  result[ [0:n_end-1:1] ] = f1[ [0:n_end-1:1] ];
  result[ [length-n_end:length-1:1] ] = f1[ [n_end:n_end*2-1:1] ];

  return result;
}
provide ("hanning");
