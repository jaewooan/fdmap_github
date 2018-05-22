function fdmap_write_body_prestress(name,mode,sxx,sxy,syy,szz,sxz,syz,endian)

% write elastic properties in file 

  if nargin<9, endian = 'n'; end
      
  prec = 'real*8';
  
  [fid,m] = fopen(name,'w',endian);
  if fid==-1,disp(m),return,end
  switch mode
    case 2
      fwrite(fid,sxx,prec);
      fwrite(fid,sxy,prec);
      fwrite(fid,syy,prec);
      fwrite(fid,szz,prec);
    case 3
      fwrite(fid,sxz,prec);
      fwrite(fid,syz,prec);
  end
  fclose(fid);