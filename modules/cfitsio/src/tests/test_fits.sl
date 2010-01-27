private variable MODULE_NAME = "cfitsio";
prepend_to_slang_load_path (".");
set_import_module_path (".:" + get_import_module_path ());

require ("fits");

private variable Failed = 0;

private define warn ()
{
   variable args = __pop_args (_NARGS);
   () = fprintf (stderr, "**** Warning: %s\n",
		 sprintf (__push_args (args)));
   Failed++;
}

private define check_key_read_write (fptr, key, val)
{
   variable val1;

   fits_update_key (fptr, key, val, NULL);
   val1 = fits_read_key (fptr, key);
   if ((val != val1) or (typeof (val) != typeof (val1)))
     warn ("failed to write and then read key %S (%S != %S)", key, 
	   typeof (val), typeof (val1));
}

define is_identical (a, b)
{
   variable dims_a, dims_b;
   (dims_a,,) = array_info (a);
   (dims_b,,) = array_info (b);
   if (length (dims_a) != length (dims_b))
     return 0;
   if (length (where(dims_a != dims_b)))
     return 0;
   if (_typeof (a) != _typeof(b))
     return 0;
   if (length (where (a != b)))
     return 0;
   return 1;
}

define test_img (filename)
{
   variable fptr = fits_open_file (filename, "c");
   variable dims = [2,10];
   variable npixels = dims[0]*dims[1];

   fits_create_image_hdu (fptr, NULL, UInt32_Type, dims);
   variable card = 
     "key_prec= 'This keyword was written by fxprec' / comment goes here";
   
   fits_write_records (fptr, [card]);
   
   variable val = "1234567890123456789012345678901234567890"
     + "12345678901234567890123456789012345";
   check_key_read_write (fptr, "card1", val);
   
   
   check_key_read_write (fptr, "keyint", 1);
   check_key_read_write (fptr, "keydbl", 1.2);
   check_key_read_write (fptr, "tstring", "a string");
   
   fits_update_logical (fptr, "tlogical", 1, NULL);
   if (1 != fits_read_key (fptr, "tlogical"))
     warn ("failed to read and write logical");
   
   
   fits_write_comment (fptr, "  This keyword was written by fxpcom.");
   fits_write_history (fptr, "    This keyword written by fxphis (w/ 2 leading spaces).");
   fits_write_date (fptr);

   % Write data 
   % define the null value (must do this before writing any data) %
   fits_update_key (fptr, "BLANK", -99, "value to use for undefined pixels");

   variable array = typecast ([1:npixels], UChar_Type);
   fits_write_img (fptr, array);
   fits_close_file (fptr);
   
   array = typecast (array, UInt32_Type);
   reshape (array, dims);

   fptr = fits_open_file (filename, "w");
   variable img = fits_read_img (fptr);
   if (0 == is_identical (img, array))
     {
	warn ("Write then read image failed: %S vs %S", array, img);
     }
   fits_close_file (fptr);
}


private define test_bt (filename)
{
   variable nrows = 71;
   variable uint16s = UInt16_Type[nrows]; uint16s[*] = [1:nrows];
   variable uint32s = UInt32_Type[nrows]; uint32s[*] = [1:nrows];
   variable xs = [1:10:#nrows*3*2];
   reshape (xs, [nrows, 3, 2]);

   variable data = struct {u16 = uint16s, u32 = uint32s, x=xs};

   fits_write_binary_table (filename, "FOO", data);
   
   variable delete = 0;

   variable table = fits_read_table (filename + "[FOO]");
   if (0 == is_identical (table.u16, data.u16))
     {
	warn ("testbt: failed to read/write an unsigned 16 bit column");
	delete = 0;
     }
   if (0 == is_identical (table.u32, data.u32))
     {
	warn ("testbt: failed to read/write an unsigned 32 bit column");
	delete = 0;
     }
   if (0 == is_identical (table.x, data.x))
     {
	warn ("testbt: failed to read/write an array column");
	delete = 0;
     }

   variable fp = fits_open_file (filename + "[FOO]", "r");
   variable r = 0, dr = 7;
   while (r < nrows)
     {
	if (r + dr > nrows)
	  dr = nrows - r;
	variable x = fits_read_col (fp, "X"; row=r+1, num=dr);
	if (0 == is_identical (xs[[r:r+dr-1],*,*], x))
	  {
	     warn ("testbt: failed to read using the row qualifier");
	     delete = 0;
	     break;
	  }
	r += nrows;
     }
   fits_close_file (fp);

   if (delete) 
     () = remove (filename);
}

test_bt ("testbt.fit");
test_img ("testimg.fit");

if (Failed == 0)
  message ("Passed");
else
  message ("Failed");

