function [interf, DEM, noise, coherence, watermask] = siminterf(Bperp,fracdim,water,maxheight,numlines,numpixels,doplot,donoise,dorefpha)
% SIMINTERF -- SIMulate phase of INTERFerogram (not amplitude).
%
%   Function to simulate an (unwrapped) interferogram by 'radarcoding'
%   a terrain model (DEM). 'Radarcoding' means converting the terrain
%   coordinates to the radar system [azimuth lines, range pixels]. 
%   The terrain model either is a geometric body (i), a fractal DEM (ii),
%   or an input matrix (iii). The range pixel spacing of the DEM is 80 meters.
%
%   The algorithm used radarcodes the DEM per line by computing the slant
%   range r to master and slave satellite for each point of the DEM, yielding
%   a function phase(r). This function is interpolated to a regular grid
%   (radar range pixel coordinate system). Zero Doppler data is assumed,
%   thus the azimuth coordinates are equal in both systems. The phase of 
%   the 'flat Earth' is computed using a far field approximation, and
%   subtracted by default. To compute r, the position of master and slave
%   satellite are fixed by using a certain Bperp and baseline orientation.
%   The range, incidence angle to the first terrain point, and the wavelength
%   are also fixed. The orbits are assumed not to converge (easy to simulate.)
% 
%   Input parameters:
%     INT = SIMINTERF by itself simulates a DEM of a cone and the
%     corresponding interferogram. It is verbose and makes plots.
%
%     SIMINTERF(Bperp) uses the specified perpendicular baseline. A larger
%     baseline corresponds with a smaller height ambiguity (more fringes).
%     A crude approximation of the height ambiguity is ha = 10000/Bperp.
%
%     SIMINTERF(Bperp,D) where D is the fractal dimension of the simulated
%     DEM. Smaller D implies smoother surface. Earth topography can be simulated
%     by a fractal dimension of approximately 2.3.
%       D == DEM     Use this input matrix as DEM; 
%                     Input arguments: water, height, lines, pixels
%                     are not assigned if this option is used.
%             0      Simulate a cone;
%            -1      Simulate a ramp;
%            -2      Simulate pyramid;
%             other  Use this fractal dimension to simulate a DEM.
%
%     SIMINTERF(Bperp,D,WATER) uses the specified percentage to 
%     create water areas of height 0 (approximately). The phase of 
%     these areas is uniform noise, the coherence is gaussian below 0.2.
%
%     SIMINTERF(Bperp,D,water,HEIGHT) simulates a DEM with the
%     specified maximum HEIGHT. (The minimum height always equals 0.)
%
%     SIMINTERF(Bperp,D,water,height,LINES) creates an interferogram
%     with the specified number of azimuth lines.
%
%     SIMINTERF(Bperp,D,water,height,lines,PIXS) creates an
%     interferogram with the specified number of range pixels.
%
%     SIMINTERF(Bperp,D,water,height,lines,pixs,DOPLOT) where
%     DOPLOT=0 prevents the generation of plots. Figures 1:6 are used.
%
%     SIMINTERF(Bperp,D,water,height,lines,pixs,doplot,DONOISE)
%     where DONOISE = 0 disables noise generation.
%
%     SIMINTERF(Bperp,D,W,H,lines,pixs,doplot,donoise,DOREFPHA)
%     where DOREFPHA is 0 if reference phase does not have to be
%     subtracted.
%
%   Defaults:
%     Bperp    = 150   (meters)
%     D        = 0     (i.e. a cone)
%     WATER    = 20    (percentage)
%     HEIGHT   = 700   (meters)
%     LINES    = 512   (number of azimuth lines)
%     PIXELS   = 512   (number of range bins in range)
%     DOPLOT   = 1     (do plot results)
%     DONOISE  = 1     (do add noise based on terrain slopes)
%     DOREFPHA = 1     (do subtract reference phase of 'flat Earth')
%
%   Additional output:
%     [INT,DEM] = SIMINTERF(...) optionally outputs the simulated DEM (in 
%     terrain coordinates).
%     [INT,DEM,NOISE] = SIMINTERF(...) optionally outputs the noise matrix
%     that is added to the INTerferogram based on the geometrical decorellation.
%     [INT,DEM,NOISE,COHERENCE] = SIMINTERF(...) optionally outputs the
%     coherence map computed from the terrain slopes (geometric decorrelation).
%     [INT,DEM,NOISE,COHERENCE,WATERMASK] = SIMINTERF(...) optionally 
%     outputs the waterarea index vector as returned by FIND.
%
%   EXAMPLES:
%   To simulate an interferogram with a baseline of 100 meter,
%   rough terrain, and 30% water area:
%     Bperp = 100; D=2.4; water=30;
%     siminterf(Bperp,D,water);
%
%   To fastly simulate some interferograms to test this script:
%     Bperp=100;  D=2.1;    water=0; H=500;
%     nlines=100; npix=200; doplot=1; donoise=1; doflat=1;
%     siminterf(Bperp,D,water,H,nlines,npix,doplot,donoise,doflat);
%
%   To radarcode an input DEM:
%     Bperp=50; X=256.*peaks(256);
%     siminterf(Bperp,X);
%
%   To view the unwrapped interferogram, and obtain the coherence map:
%     Bperp=200;  D=2.3;    water=30; H=100;
%     nlines=256; npix=256; doplot=0; donoise=1; doflat=1;
%     [INT,DEM,NOISE,COH] = ...
%       siminterf(Bperp,D,water,H,nlines,npix,doplot,donoise,doflat);
%     figure(1);
%       subplot(221), imagesc(DEM); axis 'image'; title 'DEM'; colorbar;
%       subplot(222), imagesc(INT); axis 'image'; title 'INT'; colorbar;
%       subplot(223), imagesc(wrap(INT)); axis 'image'; title 'W(INT)'; colorbar;
%       subplot(224), imagesc(COH); axis 'image'; title 'COH'; colorbar;
%
%   To make the interferogram complex, use e.g.,  
%     cint   = complex(cos(interf),sin(interf));
%   or:
%     ampl   = ones(size(interferogram));
%     cint   = ampl.*complex(cos(interf),sin(interf));
%
%   See also FRACTOPO, ANAFRACDEMO, ANGLE, COMPLEX, WRAP, CONE, PYRAMID, HEIGHTAMB
%
%   The DEOS FRACTAL and INSAR toolboxes are used
%   (www.geo.tudelft.nl/doris/)
%

% Created by Bert Kampes 05-Oct-2000
% Tested by Erik Steenbergen
%
%   TODO:
%     -more geometrical figure if fracdim=-1,-2, etc., 
%     -demo for people new to insar
%


%%% better switch more, but how, get(0)?
more off

%%% Set defaults for general parameters.
%%% Handle input; could be done smarter...
% (Bperp,dem,water,maxheight,numlines,numpixels,doplot,donoise,dorefpha)
%   1     2   3      4         5         6         7     8       9
if (nargin<9) dorefpha  = 1;    end;
if (nargin<8) donoise   = 1;    end;
if (nargin<7) doplot    = 1;    end;
if (nargin<6) numpixels = 512;  end;
if (nargin<5) numlines  = 512;  end;
if (nargin<4) maxheight = 700; end;
if (nargin<3) water     = 20;   end;
if (nargin<2) fracdim   = 0;    end;%  cone
if (nargin<1) Bperp     = 150;  end;

%%% Set output parameters.
% [interferogram, DEM, noise, coherence, watermask]
%     1            2    3       4          5
if (nargout>=3) donoise = 1;  end;
%if (nargout<1) siminterf(); end;



%%% Check if input fracdim (D) was a matrix (use it as DEM).
if (prod(size(fracdim))>1) 
  DEM                  = 999;			% 999 identifier use input matrix
  [DEM,fracdim]        = deal(fracdim,DEM);	% swap
  [numlines,numpixels] = size(DEM);
  maxheight            = max(DEM(:));
  if (nargin<3)  water = 0;   end;
end;


%%% Fixed variables (do not change w/o reason).
alpha           = deg2rad(10.);%	[rad] baseline orientation
lambda          = 0.05666;%		[m]   wavelength
theta           = deg2rad(19.);%	[rad] looking angle to first pixel
R1              = 830000.;%		[m]   range to first point
pi4divlam       = (-4.*pi)./lambda;


%%% Provide some info.
disp(['Lines (azimuth):        ', num2str(numlines)]);
disp(['Pixels (range bins):    ', num2str(numpixels)]);
disp(['Perpendicular baseline: ', num2str(Bperp), ' m']);
disp(['Height ambiguity:       ', num2str(round(heightamb(Bperp))), ' m']);
disp(['Fractal dimension DEM:  ', num2str(fracdim)]);	
disp(['Water area:             ', num2str(water), '%']);	
disp(['Minimum height DEM:     ', num2str(0), ' m']);
disp(['Maximum height DEM:     ', num2str(maxheight), ' m']);
disp(['Make plots:             ', num2str(doplot)]);	
disp(['Compute noise:          ', num2str(donoise)]);	
disp(' ');


%%% Simulation fractal DEM or cone.
switch fracdim
  case  999
    disp('Method DEM:             Using input matrix.');
  case  0
    disp('Method DEM:             Creating cone');
    DEM   = cone(numlines,numpixels);
  case -1
    disp('Method DEM:             Creating ramp');
    DEM   = ones(numlines,1) * lying(linspace(0,maxheight,numpixels));
  case -2
    disp('Method DEM:             Creating pyramid');
    DEM   = pyramid(numlines,numpixels);
  otherwise
    disp('Method DEM:             Creating fractal');
    %%% Create DEM of dimensions (numlines)x(numpixels)
    %%% See script fractopo for improvements & water areas.
    %%% and 3D visualization.
    beta = 7 - 2*fracdim;	% valid for 2D area [pr.corr.RH]
    DEM = fracsurf(numlines,beta,'n',1000);
end


%%% Add water and rescale to [0:maxheight].
if ( fracdim ~= 999 )		% input matrix
  disp(['Rescaling DEM:          [0:',num2str(maxheight),']']);
  DEM = DEM(1:numlines,1:numpixels);
  minDEM         = min(DEM(:));
  maxDEM         = max(DEM(:));
  waterlevel     = minDEM + 0.01*water.*(maxDEM-minDEM);
  DEM            = (DEM-waterlevel) .* (maxheight./(maxDEM-waterlevel));
  DEM(DEM<0)     = 0.;
end


%%% The actual radarcoding by interp1.
% use nearest, certainly for DEM fractal (?)
% due to problems with layover areas with cubic method.
%// BK 21-Nov-2000
method = 'nearest';
%method = 'linear';
%method = 'cubic';
disp(['Radarcoding DEM:        ', method ,' interpolation method']);
numpixelsdem = size(DEM,2);%			in (oversampled) DEM
dx           = 80;%    				[m] DEM resolution
%dx           = scenewidth/(numpixelsdem-1);%	[m] DEM resolution
scenewidth   = dx.*(numpixelsdem-1);%		[m]   ground range
x0           = sin(theta) .* R1;%		x coord. first DEM point
sat1_x       = 0.;%				x coord. of master satellite
sat1_y       = cos(theta) .* R1 + DEM(1,1);%	y coord. (height)
Margin       = maxheight;% [m] be on save side
maxrange     = sqrt((x0+(numpixelsdem-1)*dx).^2+sat1_y.^2)-Margin;
R1extra      = R1+Margin;
totalrange   = maxrange-R1extra;
rangebinsize = totalrange/numpixels;
rangegrid    = R1extra:rangebinsize:maxrange-rangebinsize;

%%% Give some more info.
disp(['Wavelength:             ', num2str(round(lambda*1000)./10), ' cm']);
disp(['Resolution of DEM:      ', num2str(round(dx)), ' m']);
disp(['Width of scene:         ', num2str(round(scenewidth/1e3)), ' km']);
disp(['Resolution of interf:   ', num2str(round(rangebinsize)), ' m']);
disp(['Satellite height:       ', num2str(round(sat1_y./1e3)), ' km']);
disp(['Looking angle:          ', num2str(rad2deg(theta)), ' deg']);
disp(['First range bin at:     ', num2str(round(rangegrid(1))./1e3), ' km']);
disp(['Last range bin at       ', num2str(round(rangegrid(numpixels))./1e3), ' km']);


%%% Compute range diff to slave satellite
%%% Assume range to first bin of DEM = R1
%%% use local coord. sytem to compute range to master
%satH = cos(theta) .* R1 + DEM(1,1);
B      = Bperp ./ cos(theta-alpha);
sat2_x = B .* cos(alpha);
sat2_y = B .* sin(alpha) + sat1_y;
x      = x0:dx:x0+dx*(numpixelsdem-1);%	x coord. w.r.t. sat1
x2sqr  = (x - sat2_x).^2;
xsqr   = x.^2;


%%% Compute range for whole line in 1 go.
% and compute interferometric phase, flat earth corrected.
% refphase may be wrong, but trend should be ok this way.
%%%range2master = sqrt(sat1_y.^2 + xsqr);
%%%range2slave  = sqrt(sat2_y.^2 + x2sqr);
%%%refphase     = pi4divlam .* (range2slave-range2master);
if ( dorefpha == 0 )
%%%  %refphase = zeros(size(refphase));
%%%  refphase = 0.;
  disp('Reference phase:        not subtracted')
else
  disp('Reference phase:        subtracted')
end


%%% Compute range for all lines, same refphase
oldperc = -1;
for linenum=1:size(DEM,1)
  newperc = floor((100*linenum)./numlines);
  if (newperc ~= oldperc)
    oldperc = newperc;
    %comment out on x86 PCs where '\r' does not work...
    fprintf(1,'\rRadarcoding DEM:        %3.0f%%', newperc);
  end
  y            = sat1_y-DEM(linenum,:);
  range2master = sqrt(y.^2+xsqr);
  y2           = sat2_y-DEM(linenum,:);
  range2slave  = sqrt(y2.^2+x2sqr);
  phase = pi4divlam .* (range2slave-range2master);

  if ( dorefpha ~= 0 )% 			subtract it
    tantheta = x./y2;
    deltax   = DEM(linenum,:) ./ tantheta;%	far field approx.
    x2_0     = x - deltax;
    %refpharangemaster = range2master; (approx...)
    refpharangemaster = sqrt(sat1_y.^2 +  x2_0.^2);
    refpharangeslave  = sqrt(sat2_y.^2 + (x2_0-sat2_x).^2);
    refphase = pi4divlam .* (refpharangeslave-refpharangemaster);

    phase = phase - refphase;
%not ok.
%    phase = pi4divlam .* (refpharange2slave-range2slave);
%  else
%    phase = pi4divlam .* (range2slave-range2master);
  end

  %%% Interpolate p to grid rangebins
  %%% range is not always increasing due to foreshortning
  [range2master,sortindex] = sort(range2master);
  phase                    = phase(sortindex);         
  interf_nonoise(linenum,:) = interp1(range2master,phase,rangegrid,method);

  if ( donoise ~= 0 )
    %%% Method must be nearest! (why?)
    slopeDEM = atan2((diff(DEM(linenum,:),1,2)),dx);% 1 kleiner...
    slopeDEM = [slopeDEM,0];
    slopeDEM = slopeDEM(sortindex);
    slope(linenum,:) = interp1(range2master,slopeDEM,rangegrid,'nearest');
  end
end
disp(' ');
disp(' ');

%%% Get watermask for radarcoded matrices.
watermask = [];
if ( fracdim ~= 999 )		% input matrix
  watermask = find(interf_nonoise<=10*eps);
end

%%% Adding noise
if ( donoise ~= 0 )
  disp('Adding to phase:        geometric decorrelation noise')
  [noise,coherence]  = simnoise(slope,Bperp);% rad
else
  noise     = zeros(size(interf_nonoise));
  coherence = zeros(size(interf_nonoise));
end
interferogram = interf_nonoise + noise;
%%% Tidy up a little.
clear slope;


%%% Adjusting the coherence of watermask
%%% Coherence of watermask is calculated wrong using the formula :
%%% Bcritical = lambda*(Bw/c)*R1*tan(theta-slope)
%%% since watermask doen't have any slope, Bcritical is constant
%%% watermask is given random coherence between 0.0 and 0.2
disp('Water area:             Coherence [0.0,0.2]');
coherence(watermask)=rand(size(watermask))*0.2;


if ( donoise ~= 0 )
  disp('Water area:             Phase uniform noise.');
  interferogram(watermask)=rand(size(watermask))*2*pi-pi;
else
  disp('Water area:             Phase = 0.');
  interferogram(watermask)=0;
end


%%% Plot some.
if ( doplot ~= 0 )
disp(' ');
disp('Plotting:               DEM (3D)');
%j       = jet(64);
%topomap = [(j(1,:));j(64:-1:15,:)];
%phasemap = ph(256);
j       = flipud(jet(256));
topomap = [0,0,0.5625;j(1:200,:)];
phasemap = deos(256);
figure(1);
  clf;
  %h=surf(DEM);
  %set (h,'edgecolor','none')
  mesh(DEM);
  colormap(topomap);
  l=line(round([-0.1*numpixels -0.1*numpixels]), ...
         round([-0.3*numlines 1.2*numlines]) ,...
         round([1.2*maxheight 1.2*maxheight]));
  set(l,'LineWidth',3);
  %set(l,'Linestyle','--');
  set(l,'color',[1 0 0]);
  set(l,'markersize',6);
  set(l,'marker','<');
  text(round(-0.25*numpixels),round(numlines),round(1.2*maxheight),'orbit')
  axis ('tight');
  set(gca,'Zlim',[0 1.6*maxheight]);
  %
  %view(20,65);
  %view(-15,40);
  %view(10,50);
  view(25,60);
  t=title (['DEM (height [0:',num2str(maxheight),'], ',num2str(water),'% water)']);
  set(t,'fontweight','bold');
  xlabel ('range pixels');
  ylabel ('azimuth lines');
  zlabel ('height');
  c=colorbar;
  set(get(c,'title'),'string','[meter]');


disp('Plotting:               DEM (2D)');
figure(2);
  clf;
  imagesc(DEM);
  t=title ('Digital Elevation Model');
  set(t,'fontweight','bold');
  xlabel ('range pixels');
  ylabel ('azimuth lines');
  zlabel ('height');
  axis image;
  colormap(topomap);
  c=colorbar;
  set(get(c,'title'),'string','[meter]');

disp('Plotting:               Unwrapped radarcoded DEM w/o noise');
figure(3);
  clf;
  mesh(interf_nonoise);
  axis ('tight');
  set(gca,'Zlim',[0 1.6*max(interf_nonoise(:))]);
  colormap(phasemap);
  view(25,60);
  %g = surf(interf_nonoise);
  %set (g,'edgecolor','none')
  %t=title ('Radarcoded DEM (B\perp=',num2str(Bperp),')');
  t=title (['Radarcoded DEM (Bperp=',num2str(Bperp),')']);
  set(t,'fontweight','bold');
  xlabel ('range pixels');
  ylabel ('azimuth lines');
  zlabel ('height');
  %axis image;
  colormap(topomap);
  c=colorbar;
  set(get(c,'title'),'string','[rad]');

disp('Plotting:               Unwrapped radarcoded DEM with noise');
figure(4);
  clf;
  imagesc(interferogram);
  colormap(phasemap);
  axis image;
  c=colorbar;
  set(get(c,'title'),'string','[rad]');
  t=title (['Phase of radarcoded DEM (Bperp=',num2str(Bperp),')']);
  set(t,'fontweight','bold');
  xlabel ('range pixels');
  ylabel ('azimuth lines');

disp('Plotting:               Wrapped radarcoded DEM with noise');
figure(5);
  clf;
  imagesc(wrap(interferogram));
  colormap(phasemap);
  axis image;
  c=colorbar;
  set(get(c,'title'),'string','[rad]');
  %t=title (['Wrapped phase (B\perp =',num2str(Bperp),')']);
  t=title (['Wrapped phase (Bperp =',num2str(Bperp),')']);
  set(t,'fontweight','bold');
  xlabel ('range pixels');
  ylabel ('azimuth lines');

disp('Plotting:               Coherence map');
figure(6);
  clf;
  imagesc(coherence,[0 1]);
  colormap(gray);
  t=title ('Coherence');
  set(t,'fontweight','bold');
  xlabel ('range pixels');
  ylabel ('azimuth lines');
  axis image;
  c=colorbar;
  set(get(c,'title'),'string','[-]');

% make figures active
if (exist('tip','file'))
  tip;
else
  for ii=6:-1:1
    figure(ii);
  end;
end;
end

%%% EOF
if (nargout ~= 0) 
  interf=[];
  [interf, interferogram] = deal(interferogram,interf);
end
  
more on;

