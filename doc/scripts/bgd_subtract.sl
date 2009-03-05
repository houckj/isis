public define backsub (s_idx, b_idx)
{
   variable s, b;
   variable s_backscale, s_exposure, b_backscale, b_exposure;
   variable scale;
   
   s = get_data_counts (s_idx);
   b = get_data_counts (b_idx);
   
   s_backscale = get_data_backscale (s_idx);
   b_backscale = get_data_backscale (b_idx);
   
   s_exposure = get_data_exposure (s_idx);      
   b_exposure = get_data_exposure (b_idx);      
   
   scale = (s_backscale * s_exposure) / (b_backscale * b_exposure);
   s.value -= b.value * scale;
   s.err = sqrt (s.err^2 + (b.err * scale)^2);
   
   put_data_counts (s_idx, s.bin_lo, s.bin_hi, s.value, s.err);
}
