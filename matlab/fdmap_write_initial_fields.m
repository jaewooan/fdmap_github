function fdmap_write_initial_fields(name,mode,vx,vy,vz,sxx,sxy,syy,szz,sxz,syz,endian)

% write initial fields in file 

  if nargin<12, endian = 'n'; end
      
  prec = 'real*8';
  
  [fid,m] = fopen(name,'w',endian);
  if fid==-1,disp(m),return,end
  switch mode
    case 2
      fwrite(fid,vx,prec);
      fwrite(fid,vy,prec);
      fwrite(fid,sxx,prec);
      fwrite(fid,sxy,prec);
      fwrite(fid,syy,prec);
      fwrite(fid,szz,prec);
    case 3
      fwrite(fid,vz,prec);
      fwrite(fid,sxz,prec);
      fwrite(fid,syz,prec);
  end
  fclose(fid);
