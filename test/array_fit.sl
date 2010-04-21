() = evalfile ("inc.sl");
msg ("testing array_fit.... ");

define fun (x, pars)
{
  return exp (- 0.5 * sqr((x-pars[0])/pars[1]));
}

variable
  pars = [5.0, 2.0],
  pars_min = [1.0, 0.01],
  pars_max = [10.0, 10.0];

variable x = [1.0:10.0:#100];
variable y = fun (x, pars);
variable wt = ones(length(y));

variable pars_guess = pars * (1.0 + 0.2 * urand(2));

variable stat, best;
(best, stat) = array_fit (x, y, wt, pars_guess, pars_min, pars_max, &fun);

if (length(best) != howmany(feqs(best, pars, 1.e-3)))
{
   writecol (stderr, best, pars);
   throw ApplicationError, "array_fit: fit failed";
}

msg ("ok\n");
