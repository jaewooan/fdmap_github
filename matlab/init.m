function pb=init(name,datadir,endian,psav)
% initialize variables to read in data from a problem
  
  % name
    
  pb.name = name;  
  
  % data directory
  
  if datadir(end)~='/',datadir = [datadir '/'];end
  pb.datadir = datadir;
  
  % endian
  
  if nargin<3
    pb.endian = 'n';
  else
    pb.endian = endian;
  end
  
  % precision of data to be read in

  if nargin<4
    psav = 4;
  end

  switch psav
   case 4
    pb.pre = '''real*4''';
   case 8
    pb.pre = '''real*8''';
  end
  
  % base file name

  pb.f00 = [pb.datadir,deblank(pb.name)];
  
  % load basic problem data
  
  [f m] = fopen([pb.f00,'.m',],'r');
  if ~isempty(m);disp(m);end
  if f==-1;return;end
  while 1
    txt = fgetl(f);
    if txt==-1; break; end
    eval(['pb.',txt])
  end
  fclose(f);
