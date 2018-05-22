function fdmap_write_fault_stress(name,S0,N0)
    % just make everything a row vector
    S0 = shiftdim(S0)';
    N0 = shiftdim(N0)';

    % number of elements
    nx = max([length(N0),length(S0)]);

    % expand out data that isn't varying
    if(length(S0) == 1)
        S0 = S0*ones(1,nx);
    end
    if(length(N0) == 1)
        N0 = N0*ones(1,nx);
    end
    % write file

    prec = 'real*8'; endian = 'l';

    fid = fopen(name,'w',endian);
    fwrite(fid,S0,prec);
    fwrite(fid,N0,prec);
    fclose(fid);
