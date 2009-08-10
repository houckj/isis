%#! /usr/bin/env slsh

if (__argc != 2)
{
   vmessage ("Usage:  %s <isis-install-path>", __argv[0]);
   exit(0);
}

private define is_binary ()
{
   variable file = path_concat (path_dirname (__FILE__), "../.binary");
   return stat_file (file) != NULL;
}

#ifnexists require
$1 = _slang_install_prefix + "/share/slsh/require.sl";
if (stat_file($1) == NULL)
{
   if (is_binary())
     $1 = path_dirname (__FILE__) + "/../opt/share/slsh/require.sl";
   else $1 = "require";
}
() = evalfile ($1);
#endif

$0 = __argv[1];

$1 = get_slang_load_path();
$2 = path_concat ($0, "share");
set_slang_load_path (strjoin ([$2, $1], ":"));

$1 = get_import_module_path();
$2 = path_concat ($0, "lib/modules");
$3 = path_concat ($0, "src/elfobjs");
set_import_module_path (strjoin ([$3, $2, $1], ":"));

if (NULL == stat_file (path_concat ($2, "isis-module.so"))
    && NULL == stat_file (path_concat ($3, "isis-module.so")))
  exit(0);

() = fprintf (stderr, "testing isis module import.... ");
require ("isis");
() = fprintf (stderr, "ok\n");
