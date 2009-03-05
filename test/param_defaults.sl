% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing param defaults.... ");

define plaw_fit (lo, hi, par)
{
   variable norm, p, result;

   norm = par[0];
   p    = par[1] + 1;
   result = norm * (hi^p - lo^p) / p;       % integral over bin-width

   return result;
}
typedef struct {value, freeze, min, max} Param_Info;

define set_default (x, value, freeze, min, max)
{
   x.value = value;
   x.freeze = freeze;
   x.min = min;
   x.max = max;
}

variable Defaults = Param_Info[2];
set_default (Defaults[0],  1.0, 0,  0.0,  3.0);
set_default (Defaults[1],    0, 0, -5.0 , 1.0);

define plaw_default (i)
{
   variable d = Defaults[i];
   return (d.value, d.freeze, d.min, d.max);
}

define plaw_default_args (i, arg1, arg2)
{
   if (arg1 != "a" or arg2 != 2)
     failed ("wrong default args");
   return plaw_default (i);
}

define same_config (a,b)
{
   return (a.value == b.value and a.freeze == b.freeze
           and a.min == b.min and a.max == b.max);
}

define different_param_config ()
{
   variable n, p = get_params();
   n = length(p);
   if (n != length(Defaults))
     return 0;

   variable flags = array_map (Integer_Type, &same_config, Defaults, p);

   return any(flags == 0);
}

add_slang_function ("plaw", ["norm","power"]);
fit_fun ("plaw(1)");
!if (different_param_config())
  failed ("wrong defaults [1]");

set_param_default_hook ("plaw", &plaw_default);
fit_fun ("plaw(2)");
if (different_param_config())
  failed ("wrong defaults [2]");

set_param_default_hook ("plaw", NULL);
fit_fun ("plaw(3)");
!if (different_param_config())
  failed ("wrong defaults [3]");

set_param_default_hook ("plaw", "plaw_default_args", "a", 2);
fit_fun ("plaw(4)");
if (different_param_config())
  failed ("wrong defaults [4]");

msg ("ok\n");
