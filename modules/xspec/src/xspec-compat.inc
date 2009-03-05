/* -*- mode: C; mode: fold -*- */

/* This stuff exists purely to resolve unused symbols pulled
 * in by the XSPEC libraries.  It would be nice if the XSPEC libraries
 * were more cleanly separated so that this would be unnecessary.
 *
 * Sigh.
 */

int executingScript;
char* XS_PROMPT;
const int XCM_LEN = 0;

static void xspec_call_hook (void)
{
#ifndef __APPLE__
   fprintf (stderr,
            "*** isis: xspec-module: this should never happen\n");
#endif
}

void Tcl_DeleteInterp (void *foo)
{
   (void) foo;
   xspec_call_hook();
}

void *Tcl_CreateInterp (void)
{
   xspec_call_hook();
   return NULL;
}

void *Tcl_ConcatObj (int objc, void *objv[])
{
   (void) objc; (void) objv;
   xspec_call_hook();
   return NULL;
}

char *Tcl_GetStringResult (void *foo)
{
   (void) foo;
   xspec_call_hook();
   return NULL;
}

char *Tcl_GetStringFromObj (void *p, int *n)
{
   (void) p; (void) n;
   xspec_call_hook();
   return NULL;
}

char *Tcl_SetVar (void *interp, char *varName, char *newValue, int flags)
{
   (void) interp; (void) varName; (void) newValue; (void) flags;
   xspec_call_hook();
   return NULL;
}

char *Tcl_GetString (char *objPtr)
{
   (void) objPtr;
   xspec_call_hook ();
   return NULL;
}

void Tcl_ResetResult (char *interp)
{
   (void) interp;
   xspec_call_hook ();
}

int Tcl_RegExpExec (void *interp, void *regexp, char *str, char *start)
{
   (void) interp; (void) regexp; (void) str; (void) start;
   xspec_call_hook();
   return -1;
}

int Tcl_RegExpMatch (void *interp, char *str, char *pattern)
{
   (void) interp; (void) str; (void) pattern;
   xspec_call_hook();
   return -1;
}

void Tcl_RegExpRange (void *regexp, int indx, char **startPtr, char **endPtr)
{
   (void) regexp; (void) indx; (void) startPtr; (void) endPtr;
   xspec_call_hook();
}

int Tcl_Eval (void *interp, char *string)
{
   (void) interp; (void) string;
   xspec_call_hook();
   return -1;
}

void *Tcl_RegExpCompile (void *interp, char *string)
{
   (void) interp; (void) string;
   xspec_call_hook ();
   return NULL;
}

#define MODSUB(name) \
   static void name (int ikey, float *ear, int ne, float *param, int ifl, float *photar, float *photer) \
   {  \
      (void) ikey; (void) ear; (void) ne; (void) param; (void) ifl; (void) photar; (void)photer; \
      xspec_call_hook (); \
   }
MODSUB(xacnmod)
MODSUB(xaddmod)
MODSUB(xconmod)
MODSUB(xmulmod)
FCALLSCSUB7(xacnmod,ACNMOD,acnmod,INT,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
FCALLSCSUB7(xaddmod,ADDMOD,addmod,INT,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
FCALLSCSUB7(xconmod,CONMOD,conmod,INT,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
FCALLSCSUB7(xmulmod,MULMOD,mulmod,INT,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

static void xmixmod (int ikey, float *param, float *ear, int neear, int ne, float *flxar, float *flxer)
{
   (void) ikey; (void) param; (void) ear; (void) neear; (void) ne; (void) flxar; (void) flxer;
   xspec_call_hook ();
}
FCALLSCSUB7(xmixmod,MIXMOD,mixmod,INT,FLOATV,FLOATV,INT,INT,FLOATV,FLOATV)

static int isis_stub_old_result (char *s)
{
   (void) s;
   xspec_call_hook();
   return 0;
}
FCALLSCFUN1(INT,isis_stub_old_result,OLD_RESULT,old_result,STRING)

static int isis_stub_result (char *s)
{
   (void) s;
   xspec_call_hook();
   return 0;
}
FCALLSCFUN1(INT,isis_stub_result,RESULT,result,STRING)

static int xs_reset_result (void)
{
   xspec_call_hook();
   return -1;
}
FCALLSCFUN0(INT,xs_reset_result,RERSLT,rerslt)

static void xcramff(float fnum, int nspace, int ndig, char *cbuf, int lbuf)
{
   (void) fnum; (void) nspace; (void) ndig; (void) cbuf; (void) lbuf;
   xspec_call_hook ();
}
FCALLSCSUB5(xcramff,CRAMFF,cramff,FLOAT,INT,INT,STRING,INT)

static void xplt (float *y, int *iery, int mxrow, int npts, int nvec, char *cmd, int ncmd, int ier)
{
   (void) y; (void) iery; (void) mxrow; (void) npts; (void) nvec; (void) cmd; (void) ncmd; (void) ier;
   xspec_call_hook ();
}
FCALLSCSUB8(xplt,PLT,plt,FLOATV,INTV,INT,INT,INT,STRING,INT,INT)

int xs_append_old_result (char *s)
{
   (void) s;
   xspec_call_hook();
   return 0;
}

#ifndef WITH_XSPEC_STATIC_LINKED
char *readline (char *s)
{
   (void) s;
   xspec_call_hook();
   return NULL;
}
#endif

#ifndef WITH_XSPEC_STATIC_LINKED
void add_history (char *s)
{
   (void) s;
   xspec_call_hook();
}
#endif

char *getNextLine (char *s)
{
   (void) s;
   xspec_call_hook();
   return NULL;
}

int xs_tcl_read (char *a, char *b, int *c)
{
   (void) a; (void) b; (void) c;
   xspec_call_hook();
   return 0;
}

int xs_tcl_write (char *a, int b)
{
   (void) a;  (void) b;
   /* xspec_call_hook (); */
   return 0;
}

int Tcl_PkgProvide(void *interp, const char *name, const char *version)
{
   (void) interp; (void) name; (void) version;
   xspec_call_hook();
   return 0;
}
