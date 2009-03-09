
define slsh_main ()
{
   variable test_scripts =
     ["fit", "sys_err", "par_fun", "multi", "pileup", "readcol",
      "rebin", "xgroup", "renorm", "arrayops", "ds_combine", "hist",
      "user_grid_eval", "post_model_hook", "backscale",
      "opfun", "param_defaults", "constraint", "flux_corr",
      "confmap", "rebin_dataset", "region_stats", "stat", "yshift",
      "eval_fun2", "cache", "aped_models"];

#ifdef WITH_HEADAS
   test_scripts = [test_scripts, "xspec_import"];
#endif

   variable test_pgm, prefix;
   test_pgm = __argv[1];
   prefix = __argv[2];

   vmessage ("Testing %s:", test_pgm);

   foreach (test_scripts)
     {
        variable s =();
        variable run_cmd = sprintf ("%s -n --batch %s.sl", test_pgm, s);
        if (0 != system (run_cmd))
          exit(1);
     }

   if (system (sprintf ("./here_doc.sh %s", test_pgm)))
     exit(1);
   if (system (sprintf ("./module_import.sl %s", prefix)))
     exit(1);
}

