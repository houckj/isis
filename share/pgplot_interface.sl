%
%  pgplot interface for isis
%

% 'make check' must work before 'make install'
private variable Installed =
  (NULL == stat_file (path_concat (_isis_srcdir, "configure")));

ifnot (Installed)
{
   $1 = path_concat (_isis_srcdir, "modules/pgplot/src");
   prepend_to_isis_load_path ($1);
   prepend_to_isis_module_path ($1);
}

require ("pgplot");
require ("structfuns");

private define open_plot_device (device)
{
   variable id = _pgopen(device);
   if (id <= 0)
     return -1;
   _pgask(0);
   return id;
}

private define close_plot_device ()
{
   _pgclos ();
   return 0;
}

private define subdivide (nx, ny)
{
   _pgsubp (nx, ny);
   return 0;
}

private define select_window (id)
{
   _pgslct (id);
   return 0;
}

private define select_viewport (xmin, xmax, ymin, ymax)
{
   _pgsvp (xmin, xmax, ymin, ymax);
   return 0;
}

private define set_plot_limits (xmin, xmax, ymin, ymax)
{
   _pgswin (xmin, xmax, ymin, ymax);
   return 0;
}

private define query_plot_limits ()
{
   _pgqwin();
}

private define erase ()
{
   if ("XTERM" != _pgqinf("TYPE"))
     {
        _pgeras();
     }
   else
     {
        variable esq_seq = array_to_bstring
        ([0x1b, '[', '?', '3', '8', 'h', % to Tek mode = ESC [?38h
          0x1b, 0x0c,                    % clear screen = ESC FF
          0x1b, 0x03,                    % to VT100 mode = ESC ETX
          '\0']);

        if (-1 == fputs (esq_seq, stdout))
          return -1;
     }

   return 0;
}

private define update ()
{
   _pgupdt();
   return 0;
}

private define next_page ()
{
   _pgpage();
   return 0;
}

private define set_color (c)
{
   _pgsci(c);
   return 0;
}

private define get_color ()
{
   _pgqci();
}

private define set_line_style (s)
{
   _pgsls(s);
   return 0;
}

private define set_line_width (s)
{
   _pgslw(s);
   return 0;
}

private define set_clipping (s)
{
   _pgsclp(s);
   return 0;
}

private define plot_xy (x, y)
{
   _pgline(x, y);
   return 0;
}

private define plot_points (x, y, symbol)
{
   _pgpt(x, y, symbol);
   return 0;
}

private define plot_symbol_points (x, y, symbols)
{
   _pgpnts(x, y, symbols);
   return 0;
}

private define plot_histogram (lo, hi, y)
{
   % pad with two zero-width bins to workaround a pgplot bug
   % that otherwise plots the last bin with the wrong width
   variable xx = [lo, hi[-1], hi[-1]];
   variable yy = [ y,  y[-1],  y[-1]];

   % X is the low edge of the bin
   _pgbin (xx, yy, 0);
   return 0;
}

private define plot_y_errorbar (x, top, bot, terminal_length)
{
   _pgerry (x, top, bot, terminal_length);
   return 0;
}

private define set_viewer_size (width_cm, aspect)
{
   % pgplot wants the width in inches
   _pgpap (width_cm/2.54, aspect);
   return 0;
}

private define set_char_size (size)
{
   _pgsch (size);
   return 0;
}

private define draw_box (xopt, xtick, nxsub, yopt, ytick, nysub)
{
   _pgbox (xopt, xtick, nxsub, yopt, ytick, nysub);
   return 0;
}

private define label_axes (xlabel, ylabel, tlabel)
{
   _pglab (xlabel, ylabel, tlabel);
   return 0;
}

private define put_text_xy (x, y, angle, justify, text)
{
   _pgptxt (x, y, angle, justify, text);
   return 0;
}

private define put_text_offset (where, offset, ox, oy, text)
{
   _pgmtxt (where, offset, ox, oy, text);
   return 0;
}

private define default_axis ()
{
   return "BCNST";
}

private define configure_axis (opt, is_log, has_numbers)
{
   if (opt == NULL)
     opt = "BCST";

   if (is_log == 0)
     opt = str_delete_chars (opt, "Ll");
   else if (0 == is_substr (opt, "L") && 0 == is_substr (opt, "l"))
     opt += "L";

   if (has_numbers == 0)
     opt = str_delete_chars (opt, "Nn");
   else if (0 == is_substr (opt, "N") && 0 == is_substr (opt, "n"))
     opt += "N";

   return opt;
}

private define read_cursor (c)
{
   variable x, y, ch;
   if (0 == _pgband (c.mode, c.position, c.xanchor, c.yanchor, &x, &y, &ch))
     return;  % signal an error by returning nothing
   return (x, y, ch);
}

private define init_internal_plot_interface ()
{
   variable pli = _isis->_get_plot_library_interface ();

   pli.open = &open_plot_device;
   pli.close = &close_plot_device;
   pli.subdivide = &subdivide;
   pli.select_window = &select_window;
   pli.select_viewport = &select_viewport;
   pli.set_plot_limits = &set_plot_limits;
   pli.query_plot_limits = &query_plot_limits;
   pli.erase = &erase;
   pli.update = &update;
   pli.next_page = &next_page;
   pli.get_color = &get_color;
   pli.set_color = &set_color;
   pli.set_line_style = &set_line_style;
   pli.set_line_width = &set_line_width;
   pli.set_clipping = &set_clipping;
   pli.plot_xy = &plot_xy;
   pli.plot_points = &plot_points;
   pli.plot_symbol_points = &plot_symbol_points;
   pli.plot_histogram = &plot_histogram;
   pli.plot_y_errorbar = &plot_y_errorbar;
   pli.set_viewer_size = &set_viewer_size;
   pli.set_char_size = &set_char_size;
   pli.draw_box = &draw_box;
   pli.label_axes = &label_axes;
   pli.put_text_xy = &put_text_xy;
   pli.put_text_offset = &put_text_offset;
   pli.default_axis = &default_axis;
   pli.read_cursor = &read_cursor;
   pli.configure_axis = &configure_axis;

   _isis->_set_plot_library_interface (pli);
}

private define is_hardcopy ()
{
   return ("YES" == strup (_pgqinf("HARDCOPY")));
}

private define set_standard_viewport ()
{
   _pgvstd();
   return 0;
}

private define match_viewport_to_window (xmin, xmax, ymin, ymax)
{
   _pgwnad (xmin, xmax, ymin, ymax);
   return 0;
}

private define device_id ()
{
   _pgqid();
}

private define plot_image (a, amin, amax, tr)
{
   variable dims = array_shape(a);
   _pgimag (a, 0, dims[0]-1, 0, dims[1]-1, amin, amax, tr);
   return 0;
}

private define plot_contour (a, tr, levs, line_attr)
{
   % If number of contours is
   % negative => use current line attributes
   % positive => solid lines for (+) contours, dashed for (-)
   variable dims = array_shape(a);
   _pgcont (a, 0, dims[0]-1, 0, dims[1]-1, levs, tr, line_attr * length(levs));
   return 0;
}

#ifnexists Color_Table_Type
typedef struct
{
   l,r,g,b
}
Color_Table_Type;
#endif

%  NOTE:  I'm not sure, but these image color tables
%         may require 8-bit displays...

private define init_ct (l,r,g,b) %{{{
{
   variable ct = @Color_Table_Type;
   ct.l = l;
   ct.r = r;
   ct.g = g;
   ct.b = b;
   return ct;
}

private variable Ct = Color_Table_Type[5];
Ct[0] = init_ct  %G
  (
   [0.0, 1.0],
   [0.0, 1.0],
   [0.0, 1.0],
   [0.0, 1.0]
   );

Ct[1] = init_ct %R
  (
   [-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7],
   [ 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0],
   [ 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0],
   [ 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0]
   );

Ct[2] = init_ct %H
  (
   [0.0, 0.2, 0.4, 0.6, 1.0],
   [0.0, 0.5, 1.0, 1.0, 1.0],
   [0.0, 0.0, 0.5, 1.0, 1.0],
   [0.0, 0.0, 0.0, 0.3, 1.0]
   );

Ct[3] = init_ct %W
  (
   [0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0],
   [0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0],
   [0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0],
   [0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0]
   );

Ct[4] = init_ct %A
  (
   [0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
    0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0],
   [0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
   [0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
    0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0],
   [0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   );

%}}}

private define select_color_table (id, contrast, brightness)
{
   id = abs(id) mod length(Ct);
   _pgctab (Ct[id].l, Ct[id].r, Ct[id].g, Ct[id].b, contrast, brightness);
   return 0;
}

define isis_plot_library_interface ()
{
   variable pli = _isis->_get_plot_library_interface ();

   variable extra_fields =
     ["is_hardcopy", "set_standard_viewport", "match_viewport_to_window",
      "device_id", "plot_image", "plot_contour", "select_color_table"
     ];
   pli = struct_combine (pli, extra_fields);

   pli.is_hardcopy = &is_hardcopy;
   pli.set_standard_viewport = &set_standard_viewport;
   pli.match_viewport_to_window = &match_viewport_to_window;
   pli.device_id = &device_id;
   pli.plot_image = &plot_image;
   pli.plot_contour = &plot_contour;
   pli.select_color_table = &select_color_table;

   return pli;
}

init_internal_plot_interface();

% If we're running in X-windows, specify X as
% the default plot device, unless the user has
% already specified something.
if (NULL == getenv ("PGPLOT_DEV")
    and NULL != getenv ("DISPLAY"))
  plot_device ("/xw");
