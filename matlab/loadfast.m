function data = loadfast(pb,name,nt,st)
% LOADFAST quickly load the data
%   data = LOADFAST(problem, name) loads the field name from problem (created
%   with init)
%
%   data = LOADFAST(problem, name, nt) loads only the first nt timesteps
%
%   data = LOADFAST(problem, name, [], st) load the timestep(s) indicated by st
%   where st can be vector in any order, i.e. st = [1, 20, 3, 56] loads timesteps 1,
%   20, 3, and 56 in that order (can still use nt field but is max(st) > nt will
%   throw an error)

% load data

  eval(['data = pb.',name,';'])

  if isempty(data)
    disp('field not defined')
    return
  end

  % read base time vector

  [f m] = fopen([pb.f00,'_',data.name,'.t'],'rb',pb.endian);
  if f==-1
    %if ~isempty(m);disp(m);return,end
    nt = 1;
  else
    if nargin>=3 && ~isempty(nt)
      eval(['t = fread(f,[1 nt],',pb.pre,');']);
    else
      eval(['t = fread(f,[1 inf],',pb.pre,');']);
      nt = length(t);
    end
    fclose(f);

    disp(['Total time steps = ',num2str(nt)])
    if nargin>=4 && ~isempty(st)
        if max(st) > nt
            error(['seeking to ',num2str(max(st)),' passes end at ',num2str(nt)]);
        end
    else
        st = 1:nt;
    end

    nt = length(st);

    data.t = t(st);
    data.nt = nt;
    data.tmin = min(data.t);
    data.tmax = max(data.t);
  end
  st = [0;shiftdim(st)];

  % set dimensions

  dim = 2;
  if data.nx==1;dim=dim-1;end
  if data.ny==1;dim=dim-1;end

  % open file

  filename = [pb.f00,'_',data.name,'.dat'];
  [f m] = fopen(filename,'rb',pb.endian);
  if ~isempty(m);disp([filename ':' m]);end

  if f==-1;return;end

  % read data

  switch dim

   case 0

    readdata = ['data.',data.field,' = fread(f,[1 nt],',pb.pre,');'];

   case 1

    if data.nx~=1
      eval(['data.',data.field,' = zeros(data.nx,nt);']);
      readdata = ['data.',data.field,'(:,n) = fread(f,[data.nx 1],',pb.pre,');'];
    elseif data.ny~=1
      eval(['data.',data.field,' = zeros(data.ny,nt);']);
      readdata = ['data.',data.field,'(:,n) = fread(f,[data.ny 1],',pb.pre,');'];
    end

   case 2

    eval(['data.',data.field,' = zeros(data.nx,data.ny,nt);']);
    readdata = ['data.',data.field,'(:,:,n) = fread(f,[data.nx data.ny],',pb.pre,');'];

  end

  if(strcmp(pb.pre,'''real*4'''))
    pre = 4;
  elseif(strcmp(pb.pre,'''real*8'''))
    pre = 8;
  else
    error(['You need to add the precision bytes for: ',pb.pre])
  end
  seekdata = ['if((st(n+1)-st(n)-1)~=0),fseek(f,data.nx*data.ny*pre*(st(n+1)-st(n)-1),0);end;'];

  disp(['loading ',num2str(nt),' time steps'])

  switch dim
   case {1,2}
    for n=1:nt
      eval(seekdata)
      eval(readdata)
    end
   otherwise
    eval(readdata)
  end

  fclose(f);

  % only return field array for coordinates

  switch data.field
   case {'x','y','X','Y'}
    eval(['data = data.',data.field,';'])
  end
