function fdmap_write_ratestate_friction(name,a,b,V0,f0,L,fw,Vw)

    % ensure all input vectors are same length

    N = length(a);
    assert(length(b) ==N,'all friction vectors must be of length N')
    assert(length(V0)==N,'all friction vectors must be of length N')
    assert(length(f0)==N,'all friction vectors must be of length N')
    assert(length(L) ==N,'all friction vectors must be of length N')
    assert(length(fw)==N,'all friction vectors must be of length N')
    assert(length(Vw)==N,'all friction vectors must be of length N')
    
    % write file

    prec = 'real*8'; endian = 'n';

    fid = fopen(name,'w',endian);
    fwrite(fid,a,prec);
    fwrite(fid,b,prec);
    fwrite(fid,V0,prec);
    fwrite(fid,f0,prec);
    fwrite(fid,L,prec);
    fwrite(fid,fw,prec);
    fwrite(fid,Vw,prec);
    fclose(fid);
