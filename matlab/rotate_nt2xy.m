function [wxx,wxy,wyy] = rotate_nt2xy(wtt,wnt,wnn,n)

  wxx = 2*n(1)*n(2)*wnt+n(2)^2*wtt+n(1)^2*wnn;
  wxy = -n(1)*n(2)*(wtt-wnn)+(n(2)^2-n(1)^2)*wnt;
  wyy = -2*n(1)*n(2)*wnt+n(1)^2*wtt+n(2)^2*wnn;