%
%  Use this file to provide site-local customizations.
%  If present, this file is loaded by ISIS at startup.
%
%  Site-local documentation may be placed in local_help.txt
%  using the format indicated by the provided template.
%

% provide hooks to local spectroscopy database(s):

autoload ("aped", path_concat(_isis_srcdir,"etc/aped.sl"));
autoload ("apec_d", path_concat(_isis_srcdir,"etc/apec_d.sl"));

% include paths to local modules and packages

$1 = "/nfs/cxc/a1/setup/isis-setup.sl";
if (NULL != stat_file ($1)) 
  () = evalfile ($1, current_namespace());
