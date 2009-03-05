#! /usr/bin/env slsh

if (__argc != 2)
{
   vmessage ("Usage:  %s <isis-install-path>", __argv[0]);
   exit(0);
}

$0 = __argv[1];

$1 = get_slang_load_path();
$2 = path_concat ($0, "share");
set_slang_load_path (strjoin ([$2, $1], ":"));

$1 = get_import_module_path();
$2 = path_concat ($0, "lib/modules");
set_import_module_path (strjoin ([$2, $1], ":"));

if (NULL == stat_file (path_concat ($2, "isis-module.so")))
  exit(0);

() = fprintf (stderr, "testing isis module import.... ");
require ("isis");
() = fprintf (stderr, "ok\n");
