#define MODULE_MAJOR_VERSION	0
#define MODULE_MINOR_VERSION	3
#define MODULE_PATCH_LEVEL	5
#define MKSTR1(x) #x
#define MKSTR(x) MKSTR1(x)
static char *Module_Version_String = MKSTR(MODULE_MAJOR_VERSION) "." \
   MKSTR(MODULE_MINOR_VERSION) "." MKSTR(MODULE_PATCH_LEVEL);

#define MODULE_VERSION_NUMBER	\
   (MODULE_MAJOR_VERSION*10000+MODULE_MINOR_VERSION*100+MODULE_PATCH_LEVEL)
