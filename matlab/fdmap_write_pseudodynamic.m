function fdmap_write_pseudodynamic(name,Dmax,trup,Vpeak)

    % ensure all input vectors are same length

    N = length(Dmax);
    assert(length(trup) ==N,'all vectors must be of length N')
    assert(length(Vpeak)==N,'all vectors must be of length N')
    
    % write file

    prec = 'real*8'; endian = 'n';

    fid = fopen(name,'w',endian);
    fwrite(fid,Dmax,prec);
    fwrite(fid,trup,prec);
    fwrite(fid,Vpeak,prec);
    fclose(fid);
