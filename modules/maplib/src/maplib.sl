set_import_module_path (".:" + get_import_module_path ());
#ifeval (_slang_version < 20000)
if (current_namespace () != "")
  import ("maplib", current_namespace ());
else
#endif
  import ("maplib");

#ifexists provide
provide ("maplib");
#endif
