
define apec_d ()
{
   variable db = struct 
     {
	dir, atomic_data_filemap,
	  abundance, ion_balance, 
	  line_emissivity, continuum_emissivity
     };
   
   db.dir = "/nfs/cxc/a1/share/atomdb/aped_2000oct04";
   
   db.abundance = "Abundances.fits";
   db.ion_balance = "MM98_ionbal.fits";
   
   db.line_emissivity = "APEC_TN_v101_line.fits";
   db.continuum_emissivity = "APEC_TN_v101_coco.fits";
   
   return db;
}

