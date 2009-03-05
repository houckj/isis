% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing arrayops.... ");

define test_cumsum (a, dims)
{
   variable x = _reshape(a, dims);
   variable num = length(dims);
   
   foreach ([0:num-1])
     {
	variable d = ();
	!if (max(cumsum(x,d)) == dims[d])
	  verror ("cumsum failed d = %d", d);
     }
}

variable a = ones(240);

test_cumsum (a, [240]);
test_cumsum (a, [3,80]);
test_cumsum (a, [40,6]);
test_cumsum (a, [8,5,6]);
test_cumsum (a, [2,20,6]);
test_cumsum (a, [4,5,6,2]);
test_cumsum (a, [1,3,2,4,10]);

msg ("ok\n");

