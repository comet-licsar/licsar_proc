function [] = manuwbk_cb(action)
%function [] = MANUWBK_CB(ACTION)
%
% Perform action for MANUWBK, call back functions for GUI.
%
% ACTION can be either:
%   'viewregio'    - imagesc regios
%   'viewphase'    - imagesc phase
%   'viewampli'    - imagesc ampli
%   'profileH'     - plot a horizontal profile, estimate difference
%   'profileV'     - plot a vertical profile, estimate difference
%   'profileI'     - plot a improfile
%   'climcolorbar' - click on colorbar for new CLim
%   'climregion'   - click on region for new CLim
%
%   'getdimbk'     - prompt for numlines/width
%   'setnumcycles' - prompt for number of cycles (k2pi) to add.  (logged)
%   'addcycles'    - add number of cycles to all pixels in
%                    same region of phase matrix.                (logged)
%
%   'deleteregion' - delete a region.                            (logged)
%   'joinregion'   - join 2 regions.                             (logged)
%   'splitregion'  - split a region in 2 by polynomial.          (logged)
%
%   'help'         - give help
%
% This function contains a writemfile subfunction.
% (if global variable MFILENAME not empty, opens file a lot of times...)
% saving to a file is not logged at the moment.
% we should make filenames also GLOBALS and add action: 'savefile' etc.
%
% MANUWBK_CB depends on global variables from MANUWBK
%
% $Revision: 1.5 $  $Date: 2001/09/28 14:24:32 $
%

%// BK 19-Jul-2000

% these should exist from manuwbk
global PHASE AMPLI REGIO
global LINES
global NumRegions NumCycles
global PhaseFig RegioFig AmpliFig
%global MFILENAME


%switch action
switch lower(action)
%-------------------------------------------
  case {'addcycles'}
    regionnumber = 0;
    while (regionnumber < 1)
      disp('Please click on region with mouse.');
      % note: not ok, ginput returns xy, regionnumber = regions(floor(ginput(1)))
      [y x]=ginput(1);
      regionnumber = REGIO(floor(x),floor(y));
    end
    xxx        = find(REGIO==regionnumber);
    PHASE(xxx) = PHASE(xxx) + NumCycles*2*pi;
    
    %%%log if required... 
    reg=num2str(double(regionnumber));
    writemfile(['%%% Adding cycles to region:']);
    writemfile(['XX=find(REGIO==',reg,');']);
    writemfile(['PHASE(XX)=PHASE(XX)+2*pi*',num2str(NumCycles)]);
    writemfile(['']);

%-------------------------------------------
  case 'splitregion'
    % SPLITREGION  split regions for manuwbk
    %   better make this scripts, working on global workspace data for speed?
    %   [regions, numreg]=splitregion (regionmatrix, numregions)
    %   for manuwbk callback, joins to regions by clicking with mouse...
    % logging not yet ok...
    %
    %%% Select regions A
    regionA = 0;
    while (regionA < 1)
      disp('Click on region in which you would like to create a new region.');
      disp('Right mouse button cancels.');
      [y x button]=ginput(1);
      if (button==3) disp('Split canceled.'); return; end;
      regionA = REGIO(round(x),round(y));
      if (regionA < 1) warning('You did not click in any region (#0).'); end;
    end;
    %
    %%% Polygon new region.
    disp('Select polygon region as with roipoly.');
    disp('type: ''help roiploy'' for more info).');
    [BW,xi,yi] = roipoly;%			log xi, yi
    %%% Check enough points TODO
    %if (length(xi<=2))
    %  disp('Split canceled (less than 2 points).');
    %  return;
    %end
    %    
    %%% Create new region.
    xxx=find(BW~=0 & REGIO==regionA);%		indices roi.and.in region
    REGIO(xxx)=uint8(double(NumRegions)+1);
    NumRegions=uint8(double(NumRegions)+1);
    disp(['New region created with number: ',num2str(double(NumRegions))]);
    
    %%% Logging...
    writemfile(['%%% Splitting region in two:']);
    writemfile(['regionA = ',num2str(double(regionA))]);
    writemfile(['%%% xi, yi polynomial points selected with roipoly.']);
    writemfile(['xi = [',num2str(lying(xi)),']']);
    writemfile(['yi = [',num2str(lying(yi)),']']);
    writemfile(['BW = roipoly(PHASE,xi,yi);']);
    writemfile(['%%% Create new region.']);
    writemfile(['XX = find(BW~=0 & REGIO==regionA); %indices roi.and.in']);
    writemfile(['REGIO(XX)=uint8(double(NumRegions)+1)']);
    writemfile(['NumRegions=uint8(double(NumRegions)+1)']);
    writemfile(['']);


%-------------------------------------------
  case 'joinregion'
    % JOINREGIONS  for manuwbk
    % [regions, numreg]=joinregion (regionmatrix, numregions)
    % for manuwbk callback, joins to regions by clicking with mouse...
    %
    %%% Selection of regions A,B
    regionA = 0;
    regionB = 0;
    while (regionA < 1)
      disp('Click on first region to join with second one, right mouse button returns.');
      [y x button]=ginput(1);
      if (button==3) disp('Joining canceled.'); return; end;
      regionA = double(REGIO(round(x),round(y)));
      if (regionA < 1) warning('You did not click in any region (#0).'); end;
    end;
    while (regionB < 1)
      disp('Click on second region.');
      [y x button]=ginput(1);
      if (button==3) return; end;
      regionB = double(REGIO(round(x),round(y)));
      if (regionB < 1) warning('You did not click in any region (#0).'); end;
      if (regionB==regionA)
        regionB=0;
        disp('Try again, you selected same region as first time.');
        disp('Right mouse button returns.');
      end;
    end;

    %%% Join by giving same regionnumber
    %%% Keep regions up to date.
    regionnumbers=sort([regionA,regionB]);%				min,max
    disp(['Joined regions: ',num2str(regionnumbers)]);
    %
    regionnumbers=uint8(regionnumbers);
    REGIO(find(REGIO==regionnumbers(2)))=regionnumbers(1);%	set2smallest
    if (regionnumbers(2)~=NumRegions)
      REGIO(find(REGIO==NumRegions))=regionnumbers(2);%		update
    end;
    NumRegions = uint8(double(NumRegions)-1);

    %%% Log
    writemfile(['%%% Joining regions:']);
    writemfile(['regnum=uint8([',num2str(double(regionnumbers)),'])']);
    writemfile(['REGIO(find(REGIO==regnum(2)))=regnum(1);% set2smallest']);
    writemfile(['if (regnum(2)~=NumRegions)']);
    writemfile(['  REGIO(find(REGIO==NumRegions))=regnum(2);']);
    writemfile(['end']);
    writemfile(['NumRegions = uint8(double(NumRegions)-1);']);
    writemfile(['']);


%-------------------------------------------
  case 'deleteregion'
    % DELETEREGION for manuwbk
    % this function, ask for region to be deleted,
    % set this regio to zero, adjust regios, set phase to wrapped phase,
    % set amplitude to 0, exactly as in hgt format.
    %
    %%% Select regions A
    regionA = 0;
    while (regionA < 1)
      disp('Click on region that you would like to delete.');
      disp('Right mouse button cancels.');
      [y x button]=ginput(1);
      if (button==3) return; end;
      regionA = REGIO(round(x),round(y));
      if (regionA < 1) warning('You did not click in any region (#0).'); end;
    end;
    
    %%% Delete region.
    xxx=find(REGIO==regionA);%		indices in region
    AMPLI(xxx)=0;%				reset amplitude
    PHASE(xxx)=wrap(PHASE(xxx));%		wrap phase
    REGIO(xxx)=0;%				empty out region
    if (regionA~=NumRegions)%		update region matrix
      REGIO(find(REGIO==NumRegions))=regionA;
    end;
    NumRegions=uint8(double(NumRegions)-1);
    disp(['Deleted region: ',num2str(double(regionA))]);

    %%% Log
    writemfile(['%%% Deleting region:']);
    writemfile(['regionA = ',num2str(double(regionA))]);
    writemfile(['xxx=find(REGIO==regionA)']);
    writemfile(['AMPLI(xxx)=0']);
    writemfile(['PHASE(xxx)=wrap(PHASE(xxx))']);
    writemfile(['REGIO(xxx)=0']);
    writemfile(['if (regionA~=NumRegions)']);
    writemfile(['  REGIO(find(REGIO==NumRegions))=regionA']);
    writemfile(['end']);
    writemfile(['NumRegions=uint8(double(NumRegions)-1)']);
    writemfile(['']);

%-------------------------------------------
  case 'getdimbk'
    LINES=input('Enter number of lines of files on disk: ');
    disp('Thank you.');
    
    %%%log if required... 
    writemfile(['%%% getdimbk, not logged, interactive:']);
    writemfile(['']);

%-------------------------------------------
  case 'setnumcycles'
    NumCycles=input('Enter integer number of cycles to be added: ');
    disp('Thank you.');

    %%%log if required... 
    writemfile(['%%% new number of cycles entered:']);
    writemfile(['NumCycles = ',num2str(NumCycles)]);
    writemfile(['']);


%-------------------------------------------
  case 'viewregio'
    figure(RegioFig);
    imagesc(REGIO);
    title(['regions of unwrapped data (',num2str(double(NumRegions)),')']);
    xlabel ('width');
    ylabel ('lines');
    mymap=jet(2*(double(NumRegions)+1));% otherwise not distinct? no unique colors?
    colormap(mymap);
    colorbar;

%-------------------------------------------
  case 'viewphase'
    figure(PhaseFig);
    imagesc(PHASE);
    %set(PhaseFig,'FigureUnits','pixels');% werkt niet
    %set(0,'DefaultFigureUnits','pixels');
    if (exist('pixval'))
      set(PhaseFig,'units','pixels')
      pixval(PhaseFig,'on');% image toolbox
    end; %images toolbox
    title('unwrapped phase image');
    xlabel ('width');
    ylabel ('number of lines');
    colormap(jet);
    colorbar;

%-------------------------------------------
  case 'viewampli'
    figure(AmpliFig);
    q=AMPLI.^0.3;
    q=q./(mean(q(:))/150);
    q(find(q<16))=16;
    q(find(q>255))=255;
    imagesc(q);
    %imagesc(AMPLI.^.25);
    title  ('amplitude data (exp 0.3)');
    xlabel ('width');
    ylabel ('lines');
    colormap(gray);
    colorbar;

%-------------------------------------------
  case 'profileh'
    %%% Initialize
    disp('Please enter two points.');
    figure(PhaseFig);
    [x y] = ginput(2);
    x     = round(x);
    y     = round(y);
    step  = 1;
    if (x(1)>x(2)) step=-1; end;
    %%% Handle mouse input
    profiel = PHASE(y(1),x(1):step:x(2));
    reg     = REGIO(y(1),x(1):step:x(2));
    figure;
    plot(profiel,'k');
    hold on
    plot(reg,'r--');
    % estimate num cycles
    xxx=find(reg~=0);
    regs = reg(xxx);
    reg1 = regs(1);
    reg2 = regs(length(regs));
    reg1 = find(reg==reg1);
    reg2 = find(reg==reg2);
    iareg1 = reg1(1);%					first index
    ibreg1 = reg1(length(reg1));%	last index
    iareg2 = reg2(1);%					first index
    ibreg2 = reg2(length(reg2));%	last index

    mask=zeros(size(profiel));
    mask(iareg1:ibreg1)=1;
    plot(profiel.*mask,'r');
    mask=zeros(size(profiel));
    mask(iareg2:ibreg2)=1;
    plot(profiel.*mask,'b');

    last = max(iareg1,ibreg1);
    first= min(iareg2,ibreg2);
    if (last>first)
      last = min(iareg1,ibreg1);
      first= max(iareg2,ibreg2);
    end
    d = abs(profiel(last)-profiel(first));
    c = round(d/(2*pi));
    title(['profile; difference: ',num2str(d),' = ',num2str(c),' + ', ...
            num2str((d-c*2*pi)/(2*pi)),' cycles']);
    hold off

%-------------------------------------------
  case 'profilev'
    %%% Initialize
    disp('Please enter two points.');
    figure(PhaseFig);
    [x y] = ginput(2);
    x     = round(x);
    y     = round(y);
    step  = 1;
    if (y(1)>y(2)) step=-1; end;
    %
    %%% Handle mouse input.
    profiel = PHASE(y(1):step:y(2),x(1));
    reg     = REGIO(y(1):step:y(2),x(1));
    %%%
    figure;
    plot(profiel,'k');
    hold on
    plot(reg,'r--');
    % estimate num cycles
    xxx=find(reg~=0);
    regs = reg(xxx);
    reg1 = regs(1);
    reg2 = regs(length(regs));
    reg1 = find(reg==reg1);
    reg2 = find(reg==reg2);
    iareg1 = reg1(1);%					first index
    ibreg1 = reg1(length(reg1));%	last index
    iareg2 = reg2(1);%					first index
    ibreg2 = reg2(length(reg2));%	last index

    mask=zeros(size(profiel));
    mask(iareg1:ibreg1)=1;
    plot(profiel.*mask,'r');
    mask=zeros(size(profiel));
    mask(iareg2:ibreg2)=1;
    plot(profiel.*mask,'b');

    last = max(iareg1,ibreg1);
    first= min(iareg2,ibreg2);
    if (last>first)
      last = min(iareg1,ibreg1);
      first= max(iareg2,ibreg2);
    end
    d = abs(profiel(last)-profiel(first));
    c = round(d/(2*pi));
    title(['profile; difference: ',num2str(d),' = ',num2str(c),' + ', ...
                 num2str((d-c*2*pi)/(2*pi)),' cycles']);
    hold off

%-------------------------------------------
  case 'profilei'
    %%% Initialize
    disp('Please enter two points.');
    figure(PhaseFig);
    [x y]=ginput(2);
    x=round(x);
    y=round(y);
    figure;
    improfile(PHASE,x,y);%		check options, not in standard

%-------------------------------------------
  case 'help'
    help manuwbk;

%-------------------------------------------
  case 'climregion'
    %figure(RegioFig);
    figure(PhaseFig);
    regionA=0;
    while(regionA<1)
      disp('Click in region where you like to maximize clim.');
      [y x] = ginput(1);
      regionA = double(REGIO(round(x),round(y)));
      if (regionA < 1) warning('You did not click in any region (#0).'); end;
    end
    newclim=zeros(1,2);
    xxx=find(REGIO==regionA);
    newclim(1)=min(min(PHASE(xxx)));
    newclim(2)=max(max(PHASE(xxx)));
    disp(['CLim: ', num2str(newclim)]);
    %figure(PhaseFig);
    set(gca,'CLim',newclim);
    colorbar;

%-------------------------------------------
  case 'climcolorbar'
    figure(PhaseFig);
    disp('Click with left button an interval on the colorbar (phase).');
    [x newclim] = ginput(2);
    newclim=lying(sort(newclim));
    disp(['CLim: ', num2str(newclim)]);
    figure(PhaseFig);
    set(gca,'CLim',newclim);
    %imagesc(PHASE,newclim);
    colorbar;

%-------------------------------------------
  otherwise
    error(['Unknown action: ', action])
end



%-------------------------
%-------------------------
%-------------------------
function [] = writemfile(s)
%function [] = writemfile(char_string)
%writemfile subfunction assumes global fid of opened file or -1
%to make a log...
%for some stupid reason fid cannot be passed and file has
%to be opened/closed/opened etc?
%// BK 19-Jul-2000

global MFILENAME

%if (MFILEFID>0)%		checking if logging is requested, thus opened ...
%  fprintf(MFILEFID, '%s\n',s);
%end
if (~isempty(MFILENAME))
  %disp(['logging to: ', MFILENAME]);
  MFILEFID=fopen(MFILENAME,'at');%       append in text mode
  %if (MFILEFID<0)
  %  error('logmfile not opened ok, permissions?');
  %end;
  fprintf(MFILEFID, '%s;\n',s);
  fclose(MFILEFID);
end

