require("xspec");

variable names, lo, hi, id, fp;
fp = fopen (Xspec_Model_Names_File, "r");
names = fgetslines (fp);
() = fclose (fp);
names = array_map (String_Type, &strtrim, names);

(lo, hi) = linear_grid (1,48, 1024);

id = plot_open ("test.ps/cps", 5,5);
fp = fopen ("test.out", "w");
if (fp == NULL)
  exit ("Failed opening test.out for writing");

% need non-zero redshift for cflow models
variable fixup = Assoc_Type[];
fixup ["cflow"] = "set_par (6,0.1)";
fixup ["mkcflow"] = "set_par (5,0.1)";
fixup ["vmcflow"] = "set_par (18,0.1)";

% convolution models are different
variable conv = Assoc_Type[];
conv ["gsmooth"] = "";
conv ["lsmooth"] = "";
conv ["reflect"] = "";
conv ["rgsxsrc"] = "";
conv ["kdblur"] = "";
conv ["kdblur2"] = "";
conv ["rdblur"] = "";
%conv ["kerrconv"] = "";
conv ["partcov"] = "";
conv ["simpl"] = "";
conv ["cflux"] = "";

% FIXME - what's up with these?
variable skip = Assoc_Type[];
% skip ["equil"] = "";
% skip ["kerrd"] = "";
% skip ["plcabs"] = "";
% skip ["smaug"] = "";
skip ["rgsxsrc"] = "";
skip ["ezdiskbb"] = "";
skip ["gnei"] = "";
skip ["kerrdisk"] = "";
skip ["kerrconv"] = "";  % convolution model
skip ["nsmax"] = "";

%skip ["sedov"] = "";   % causes invalid read
%skip ["srcut"] = "";   % causes invalid read?
%skip ["sresc"] = "";   % causes invalid read?

% Note that c6mekl and c6vmekl compute a DEM which is
% not constrained to be positive.

() = fprintf (fp, "%20s  %12s\n", "Name", "Flux");

variable n, flux;

foreach (names)
{
   n = ();

   if (assoc_key_exists (skip,n))
     {
        vmessage ("skipping %s",n);
        continue;
     }
   else message (n);

   if (assoc_key_exists (conv, n))
     fit_fun (n + "(1,mekal(1))");
   else
     fit_fun (n + "(1)");

   if (assoc_key_exists (fixup, n))
     eval(fixup[n]);

   variable e = NULL;
   try (e)
     {
        flux = eval_fun (lo, hi);
     }
   catch AnyError:
     {
        vmessage ("failed evaluating %s", n);
        print(e);
     }

   () = fprintf (fp, "%20s  %12.4e\n", n, sum(flux));

   label ("Wavelength [Angstrom]", "Flux", n);
   hplot (lo, hi, flux);
}

() = fclose (fp);
plot_close (id);

vmessage ("DONE");
