% -*- mode: SLang; mode: fold -*-
() = evalfile ("inc.sl");
msg ("testing fs_comm.... ");

private define test_items()
{
   variable
     chars = ['a', 'b', 'c', 'd', 'e'],
     ints = [1:5],
     longs = [0L:5L],
     floats = typecast (urand(5), Float_Type),
     doubles = urand(5),
     strings = ["this", "is", "an", "array", "of", "strings"],
     nulls = [NULL, NULL, NULL],
     doubles_3d = urand(3,5,2);

   variable _struct = struct
     {
        cc = chars, ii = ints, ll = longs, ff = floats, dd = doubles, ss = strings,
        nulls = nulls, doubles_3d = doubles_3d
     };

   variable _list =
     {
        chars, ints, longs, floats, doubles, strings, nulls, doubles_3d
     };

   variable _assoc = Assoc_Type[];
   _assoc["chars"] = chars;
   _assoc["ints"] = ints;
   _assoc["longs"] = longs;
   _assoc["floats"] = floats;
   _assoc["doubles"] = doubles;
   _assoc["strings"] = strings;
   _assoc["nulls"] = nulls;
   _assoc["doubles_3d"] = doubles_3d;

   return _list, _struct, _assoc;
}

private variable _list, _struct, _assoc;
(_list, _struct, _assoc) = test_items();

private variable TEST_ITEMS = 1000;

private define send_test_items (s)
{
   send_msg (s, TEST_ITEMS);

   variable buf = __fs_initsend ();
   variable x;
   foreach x (_list)
     {
        __fs_pack (buf, x);
     }
   __fs_pack (buf, _list);
   __fs_pack (buf, _struct);
   __fs_pack (buf, _assoc);

   __fs_pack (buf, {_list, _struct, _assoc});
   __fs_pack (buf, struct {_list=_list, _struct=_struct, _assoc=_assoc});

   variable aa = Assoc_Type[];
   aa["_list"] = _list;
   aa["_struct"] = _struct;
   aa["_assoc"] = _assoc;
   __fs_pack (buf, aa);

   __fs_pack (buf, [_struct, _struct, _struct]);
   __fs_pack (buf, [_list, _list, _list]);
   __fs_pack (buf, [_assoc, _assoc, _assoc]);

   __fs_pack (buf, _list, _struct, _assoc);

   return __fs_send_buffer (s.fp, buf);
}

private define sender (s)
{
   loop (3)
     {
        () = send_test_items (s);
     }

   return 0;
}

private define check_list (_l)
{
   variable i, n = length(_l);
   _for i (0, n-1, 1)
     {
        if (any(_l[i] != _list[i]))
          throw ApplicationError;
     }
}

private define check_struct (_s)
{
   variable fn, fnames = get_struct_field_names (_struct);
   foreach fn (fnames)
     {
        if (any (get_struct_field (_s, fn) != get_struct_field (_struct, fn)))
          throw ApplicationError;
     }
}

private define check_assoc (_a)
{
   variable kn, knames = assoc_get_keys (_assoc);
   foreach kn (knames)
     {
        if (any (_a[kn] != _assoc[kn]))
          throw ApplicationError;
     }
}

private define handler (s, msg)
{
   variable buf;

   switch (msg.type)
     {
      case TEST_ITEMS:
        buf = __fs_recv_buffer (s.fp);
        variable l, x;
        foreach l (_list)
          {
             x = __fs_unpack (buf);
             if (any (x != l))
               throw ApplicationError;
          }

        variable _l = __fs_unpack (buf);
        check_list (_l);

        variable _s = __fs_unpack (buf);
        check_struct(_s);

        variable _a = __fs_unpack (buf);
        check_assoc (_a);

        variable ll = __fs_unpack (buf);
        check_list (ll[0]);
        check_struct (ll[1]);
        check_assoc (ll[2]);

        variable ss = __fs_unpack (buf);
        check_list (ss._list);
        check_struct (ss._struct);
        check_assoc (ss._assoc);

        variable aa = __fs_unpack (buf);
        check_list (aa["_list"]);
        check_struct (aa["_struct"]);
        check_assoc (aa["_assoc"]);

        variable a_s = __fs_unpack (buf);
        foreach x (a_s)
          {
             check_struct (x);
          }

        variable a_l = __fs_unpack (buf, 3);
        foreach x (a_l)
          {
             check_list (x);
          }

        variable a_a = __fs_unpack (buf, 3);
        foreach x (a_a)
          {
             check_assoc (x);
          }

        ll = __fs_unpack (buf, 3);
        check_list (ll[0]);
        check_struct (ll[1]);
        check_assoc (ll[2]);
     }
}

define isis_main ()
{
   variable slaves = new_slave_list ();
   loop (3)
     {
        variable s = fork_slave (&sender);
        append_slave (slaves, s);
     }
   manage_slaves (slaves, &handler);
   msg ("ok\n");
}
