% This file exists only so that require("isis") 
% will import the isis module.
import("isis");

$1 = path_concat (getenv("HOME"), ".isisrc");
if (NULL != stat_file ($1))
  () = evalfile ($1, current_namespace());

provide("isis");
