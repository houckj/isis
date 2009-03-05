variable e, s, n, lo, hi, t;
require ("xspec");

e = 10.0^[-1:1:0.005];                 % make a log energy grid
s = _mekal (5.0, e);                   % 5 keV, solar abundance
                                        % MEKAL spectrum

 n = length (s);                        % The e array has n+1 elements
lo = e[[0:n-1]];                       % low edge of histogram bin
hi = e[[1:n]];                         % high edge of histogram bin

variable id = open_plot ("xspec.ps/vcps");
resize(15);

xlog; ylog;
label ("log E (keV)", "log Flux", "");
hplot (lo, hi, s);                     % plot log-log

t = _wabs (0.01, e);                   % include Nh = 1.e20 cm^-2
ohplot (lo, hi, t * s);                % overplot absorbed spectrum

close_plot (id);
