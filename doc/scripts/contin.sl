
plasma(aped);

variable lo, hi, p, pid, ymax;

pid = open_plot ("contin.ps/vcps");
resize(15);

(lo, hi) = linear_grid(1,20, 3000);
p = get_contin (lo, hi, 3.e6);

print(p);

ymax = max ([p.true, p.pseudo]);

yrange (ymax*1.e-4, ymax);

ylog;

label ("Wavelength [Angstrom]", "log Flux [photon cm\\u3\\d s\\u-1\\d]", "T = 3.0e6 K");
hplot (lo, hi, p.true);

ohplot (lo, hi, p.pseudo);

variable q;
q = get_contin (lo, hi, 3.e6, 1.0e-3, 8);        % Contribution from Oxygen
ohplot (lo, hi, q.true);
ohplot (lo, hi, q.pseudo);

close_plot (pid);
