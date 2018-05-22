function [wx,wy] = rotate_nt2xy_vec(wt,wn,n)

  wx =  n(2)*wt+n(1)*wn;
  wy = -n(1)*wt+n(2)*wn;
