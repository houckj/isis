% Note that there's nothing special about this function.
% The whole point is to load a structure with some database
% filenames.  Using a function to generate the structure just
% allows some extra versatility.
%
% Its ok if some of the files are missing -- ISIS should
% be able to work with whatever subset of these files you
% may have.

% Usage:
%
% isis>  atoms(aped);
% isis>  plasma(aped);

define get_version_string (dir);  % prototype
define get_database_directory ();

define aped ()
{
   variable db = struct
     {
	  dir,                  % The ATOMDB directory
	  atomic_data_filemap,  % An ascii file, usually "filemap"
	  abundance,            % Elemental abundances (FITS file)
	  ion_balance,          % Ionization balance (FITS file)
	  line_emissivity,      % Line emissivity tables (FITS file
	  continuum_emissivity  % Continuum emissivity tables (FITS file)
     };

   db.dir = get_database_directory ();
   if (db.dir == NULL)
     return NULL;

   variable version = get_version_string (db.dir);
   if (version == NULL)
     {
	vmessage ("*** Failed loading ATOMDB from %s", db.dir);
	return NULL;
     }

   if (Isis_Verbose >= 0)
     vmessage ("Reading database version %s from %s", version, db.dir);

   % versions >= 1.2.0 don't have this problem.
   
   if (version == "1.1.0")
     version = version[[[0:2],4]];

   % All these files are assumed to be in db.dir.
   % If the specified files are not all within db.dir or
   % one of its subdirectories, define db.dir=""
   % (e.g. the empty string) and provide the full
   % path to each individual file:
   db.atomic_data_filemap = "filemap";

   db.abundance = "APED/misc/Abundances.fits";

   variable ion_balance = "APED/ionbal/v" + version + "_ionbal.fits";
   if (NULL == stat_file (path_concat (db.dir, ion_balance)))
     {
        ion_balance = "APED/ionbal/MM98_ionbal.fits";
     }
   db.ion_balance = ion_balance;

   db.line_emissivity = "apec_v" + version + "_line.fits";
   db.continuum_emissivity = "apec_v" + version + "_coco.fits";

   return db;
}

define get_database_directory ()
{
   variable not_found = "\n*** ATOMDB not found.\n*** Is the ATOMDB environment variable set correctly?\n";

   % Note that db.dir is taken as the definition for the $ATOMDB
   % variable which appears in the atomic_data_filemap file.
   % If the $ATOMDB string doesn't appear in your 'filemap' file,
   % this may not work.  Sorry for any trouble that may cause,
   % but its not under my control.
   variable dir = getenv ("ATOMDB");
   if (orelse
       {dir == NULL}
       {stat_file (dir) == NULL})
     {
	message (not_found);
	return NULL;
     }

   return dir;
}

define get_version_string (dir)
{
   variable fp = fopen (dir + "/VERSION", "r");
   if (fp == NULL)
     {
	vmessage ("*** File not found:  %s/VERSION", dir);
	return NULL;
     }

   variable version;
   if (-1 == fgets (&version, fp))
     {
	vmessage ("*** Failed reading %s/Version", dir);
	return NULL;
     }

   return strtrim(version);
}
