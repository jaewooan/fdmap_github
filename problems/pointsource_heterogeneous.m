% script to generate heterogeneous material properties for point source test,
% Ossian's 2018 SCEC project on free surface BC

% set up or load grid

grid_option = 1;
switch grid_option
  case 1 % load from homogeneous simulation
    A = init('point_r1','~/fdmap/data');
    x = loadfast(A,'B_x'); y = loadfast(A,'B_y');
  case 2 % set up grid by hand
    %x = ...; y = ...;
end

% set material properties
[nx,ny] = size(y); rho = zeros(nx,ny); cs = zeros(nx,ny); cp = zeros(nx,ny);
for j=1:ny
  for i=1:nx
    if y(i,j)<-7
      rho(i,j) = 2.8; cs(i,j) = 3.464; cp(i,j) = 6;
    elseif y(i,j)<-3
      rho(i,j) = 2.3; cs(i,j) = 2.7; cp(i,j) = 5;
    elseif y(i,j)<-1.5
      rho(i,j) = 1.5; cs(i,j) = 1.0; cp(i,j) = 4;
    else
      rho(i,j) = 1.0; cs(i,j) = 0.8; cp(i,j) = 3.6;
    end
  end
end

name = 'point_heterogeneous';
fdmap_write_elastic(name,rho,cs,cp);
