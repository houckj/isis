% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing assign_model.... ");

variable l, h, c, n = 100, value = 0.05;
(l, h) = linear_grid (1,10,n);
c = ones(n)*value;
() = define_counts (l, h, c, sqrt(c));
() = define_counts (l, h, c, sqrt(c));

define unequal_values (x, y)
{
   return howmany (x[0].value != y[0].value
                   or x[1].value != y[1].value);
}

define cnst_fit (l,h,p)
{
   variable v = @l;  v[*] = p[0];
   return v;
}
add_slang_function ("cnst", "factor");

define old_method ()
{
   return cnst (Isis_Active_Dataset) + cnst(3);
}

fit_fun ("old_method");
set_par ("cnst(1).factor", 0.8*value);
set_par ("cnst(2).factor", 1.2*value);
set_par ("cnst(3).factor", 1.0*value);

() = eval_counts;
variable m0 = array_map (Struct_Type, &get_model_counts, [1,2]);

% Test models defined by a string expression
fit_fun ("null");
assign_model (1, "cnst(1) + cnst(3)");
assign_model (2, "cnst(2) + cnst(3)");

() = eval_counts;
variable m1 = array_map (Struct_Type, &get_model_counts, [1,2]);

ifnot (0 == unequal_values (m0, m1))
   throw ApplicationError, "m1 test failed";

% Test models defined by a Ref_Type, WITH arg-passing
define alt_method (id, arg_not_used)
{
   return cnst(id) + cnst(3);
}
assign_model (1, &alt_method, 1, "notused");
assign_model (2, &alt_method, 2, "notused");

() = eval_counts;
variable m2 = array_map (Struct_Type, &get_model_counts, [1,2]);
ifnot (0 == unequal_values (m0, m2))
   throw ApplicationError, "m2 test failed";

% Test models defined by a Ref_Type, WITHOUT arg-passing
define alt_method2 ()
{
   return cnst(Isis_Active_Dataset) + cnst(3);
}
assign_model (1, &alt_method2);
assign_model (2, &alt_method2);

() = eval_counts;
variable m3 = array_map (Struct_Type, &get_model_counts, [1,2]);
ifnot (0 == unequal_values (m0, m3))
   throw ApplicationError, "m3 test failed";

msg ("ok\n");
