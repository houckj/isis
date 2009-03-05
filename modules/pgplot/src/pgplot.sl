import ("pgplot");

$1 = getenv ("PGPLOT_DIR");
if (NULL == stat_file (path_concat ($1, "pgxwin_server")))
{
   $2 = path_concat ($1, "../bin/pgxwin_server");
   if (NULL != stat_file($2))
     {
        $2 = getenv ("PATH");
        putenv (sprintf ("PATH=%s:%s/../bin", $2, $1));
     }
}

provide ("pgplot");

