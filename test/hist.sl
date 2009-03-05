% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing hist.... ");

define test_hist1d (n, m) %{{{
{
   variable pts = urand (n); pts = pts[array_sort(pts)];
   variable edges = 0.1 * [1:m];  % [0.1, 0.2, 0.3, 0.4, 0.5]
   variable hi_with_overflow = [edges[[1:m-1]], _isis->DBL_MAX];
   variable max_edges = max(edges);
   variable rev_indices;

   variable h = histogram (pts, edges, hi_with_overflow, &rev_indices);
   variable len = length (where (pts >= max_edges));
   if (len != h[-1])
     failed ("histogram(%S,%S); last bin: expect %d, found %d",
	     n, m, len, h[-1]);

   len = length (where (pts < edges[0]));
   if (len + sum (h) != length (pts))
     failed ("histogram: total number expect is wrong");

   _for (0, m-1, 1)
     {
	variable i = ();
	variable j = rev_indices[i];

	!if (length (j))
	  continue;

	if (length (where ((j < 0) or (j >= n))))
	  failed ("histogram: reverse index out of range");

	variable p = pts[j];

	if (length (where (p < edges[i])))
	  failed ("histogram: reverse index problem 2");

	if (i != m-1)
	  {
	     if (length (where (p >= edges[i+1])))
	       failed ("histogram: reverse index problem 2");
	  }
     }
}

%}}}

define test_hist2d (num, nr, nc) %{{{
{
   variable r = urand(num);
   variable c = urand(num);

   variable gr, gc;
   gr = [0.0:nr]/nr;
   gc = [0.0:nc]/nc;

   variable rev, img;
   img = histogram2d (r, c, gr, gc, &rev);

   % all data points got binned
   if (sum(img) != num)
     verror ("Failed:  wrong histogram sum");

   % the reverse indices include every point
   variable i, ir, ic;
   foreach ([0:num-1])
     {
	i = ();
	ir = int(r[i] * nr);
	ic = int(c[i] * nc);
	if (1 != length(where(rev[ir,ic] == i)))
	  verror ("Failed: wrong bin entry");
     }
   
   % the reverse index arrays match the image
   variable size = nr*nc;
   variable _r = _reshape(rev, size);
   variable _img = _reshape(img, size);
   foreach ([0:size-1])
     {
	i = ();
	if (length(_r[i]) != _img[i])
	  verror ("Failed: reverse indices don't match image");
     }
   
   % there are no extra reverse indices (redundant)
   variable m, s = 0;
   foreach (rev)
     {
	m = ();
	s += length(m);
     }   
   if (s != num)
     verror ("Failed: wrong number of reverse indices");
}

%}}}

test_hist1d (20, 5);
test_hist1d (20, 4);
test_hist1d (20, 3);
test_hist1d (20, 2);
test_hist1d (20, 1);

test_hist2d (20, 5, 1);
test_hist2d (20, 5, 2);
test_hist2d (20, 5, 3);
test_hist2d (20, 5, 4);
test_hist2d (20, 5, 5);

test_hist2d (20, 1, 5);
test_hist2d (20, 2, 5);
test_hist2d (20, 3, 5);
test_hist2d (20, 4, 5);
test_hist2d (20, 5, 5);

msg ("ok\n");
