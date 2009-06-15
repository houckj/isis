() = evalfile ("inc.sl");
msg ("testing ignore/notice_values.... ");

private define check_datasets (indices, list)
{
   variable i;
   foreach i (indices)
     {
        variable x = get_data_info (i);
        if (any (list != x.notice))
          throw ApplicationError;
     }
}

define isis_main()
{
   variable lo, hi, n=10, c = ones(n)*10;
   (lo,hi) = linear_grid (1,n,n);
   loop (5) () = define_counts (lo, hi, c, sqrt(c));

   variable a = [2.5, 5.0], b=[4.5, 9.0];

   variable datasets = [1:5];
   variable outside, inside = where ((a[0] <= lo and hi <= b[0])
                                     or (a[1] <= lo and hi <= b[1]), &outside);
   variable list = Integer_Type[n];

   notice_values (datasets, a[0], b[0], a[1], b[1] ; unit="a");
   list[inside] = 1; list[outside] = 0;
   check_datasets (datasets, list);

   notice_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_sum=51);
   list[*] = 0;
   check_datasets (datasets, list);

   notice_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_sum=50);
   list[inside] = 1; list[outside] = 0;
   check_datasets (datasets, list);

   notice_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_val=11);
   list[*] = 0;
   check_datasets (datasets, list);

   notice_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_val=10);
   list[inside] = 1; list[outside] = 0;
   check_datasets (datasets, list);

   notice_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_sum=50, min_val=10);
   list[inside] = 1; list[outside] = 0;
   check_datasets (datasets, list);

   notice (datasets);
   ignore_values (datasets, a[0], b[0], a[1], b[1] ; unit="a");
   list[inside] = 0; list[outside] = 1;
   check_datasets (datasets, list);

   notice (datasets);
   ignore_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_sum=51);
   list[*] = 1;
   check_datasets (datasets, list);

   notice (datasets);
   ignore_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_sum=50);
   list[inside] = 0; list[outside] = 1;
   check_datasets (datasets, list);

   notice (datasets);
   ignore_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_val=11);
   list[*] = 1;
   check_datasets (datasets, list);

   notice (datasets);
   ignore_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_val=10);
   list[inside] = 0; list[outside] = 1;
   check_datasets (datasets, list);

   notice (datasets);
   ignore_values (datasets, a[0], b[0], a[1], b[1] ; unit="a", min_sum=50, min_val=10);
   list[inside] = 0; list[outside] = 1;
   check_datasets (datasets, list);

   msg ("ok\n");
}

