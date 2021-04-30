% MANUWBK -- manually unwrap tool (script)
%
% Interactively correct unwrapped regions.
%   To be used with treefbk unwrapping program.
%   (FORTRAN, unwrap regions by treef, saves regionsfile).
% The problem with unwraping is that it in general does not unwrap
%   the total area, but rather a number of isolated patches (regions).
%   There is a k2pi offset between these regions.
% This utility supports estimation of these offsets, and correction
%   of the unwrapped phase.
% 
% Alternatively, this tool can be used to support the poor man's
%   unwrapping algorithm; unwrap by digitizing fringes. To use this
%   option, convert your complex interferogram to the HGT format, and 
%   create a regionsfile with all ones (not unwrapped, 1 region, see below).
%   With this tool you now can split this region by following the fringes
%   and clicking with the mouse. After you have done this you can add
%   multiples of 2pi to the newly created regions.
%
% The hgtfile (see FREADHGT) contains the amplitude and phase in
%   a 2 channel line interleaved file in major row order (float32).
%   If a pixel is not unwrapped, then the amplitude equals 0, and the
%   wrapped is stored.
% The regionfile (short int) contains flags for unwrapped areas.
%   0 means that this pixel is not in any region, otherwise, a positive
%   number indicates the region number. The maximum number indicates the 
%   total number of regions (no gaps).
%
%
% MANUWBK works with global variables:
%     - PHASE, REGIO, AMPLI, LINES, WIDTH, NumCycles, etc.
%
%    You can change these variables, manipulate them. e.g. to display
%    the phase only in unwrapped areas, one could say:
%      PHASE(find(AMPLI==0))=NaN;%               click on view phase
%    Or make a histogram to detect outliers:
%      q=hist(PHASE(:),[-20:20]); plot(-20:20,q);
%    Or to get a better view, one could use the zoom option and
%      figure(1); caxis([0 10]); colorbar;
%
%
% MANUWBK works with variables in workspace (case sensitive):
%   (Which you should define yourself...)
%     - manuwbkdebug = 1       Debug/test: creates dummy HGT and REGIONS matrices
%     - LOGMFILE   = 'doit.m'  Echo input to mfile "doit.m" which, if executed,
%				repeats the processing (for non-interactive actions),
%                               thus reproducing the results.
%                               This avoids the need to save the result.
%     - hgtfile    = 'uinf.hgt' hgt file with phase/amplitude info of
%                                unwrapped interferogram.
%     - regionfile = 'reg.raw' int16 file with regions.
%     - lines      = 658;       number of lines in hgt/regionsfile (or prompted for).
%
% If you want to log to a mfile, note that not interactive things are not
%   repeated. The function split region is also not logged yet.
%   Perhaps it is now (BK?)
%
%
% USAGE:
%   First set numbers, then perform actions. (Make this into buttons later...)
%   E.g., first enter dimension of input files, then load them.
%   and first enter number of cycles to add, then perform addition.
%   If this utility proves useful, we will make interactive bars to set these things.
%
%
% FUNCTION DEFINED IN SUBMENUS:
%   ==Tools==
%      - Add multiple of 2pi to a region;
%      - Plot profiles to estimate offset between unwrapped regions;
%      - Load/save files, quit;
%   ==View==
%      - View unwrapped phase, amplitude of interferogram, or regions map;
%      - Zoom in to maximize contrast;
%   ==Regions==
%      - Join together, split up, or delete a region;
%        An arbitrary region can be deleted by first creating a new region with
%        the function 'split', and thereafter deleting this new region.
%   ==Info==
%      - About;
%      - Help;
%
%
% EXAMPLE to test this script:
%   manuwbkdebug=1; manuwbk
%
%
% See also: INSAR toolbox, IMPROFILE, UINT8, ...
%   MANUWBK requires IMAGES toolbox for some functions (improfile, roipoly, ?)
%

% $Revision: 1.14 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 03-Mar-2000


%%% Globals for manuwbk_cb callback functions
global PHASE AMPLI REGIO
global LINES 
global NumRegions NumCycles
global PhaseFig RegioFig AmpliFig
global MFILENAME


%%% Initialize
more off;
close all;
manuwbkversion = '1.0';

%manuwbkdebug = 1;	% check if declared in workspace
if (exist('manuwbkdebug','var'))%			Testing only
  if (exist('LOGMFILE','var'))
    warning ('LOGMFILE declared for debugging?');
  end;
  %%% Defaults
  disp('Debugging...');
  disp('generating random phase, amplitude, and region matrices.');
  hgtfile    = 'uint.hgt';
  regionfile = 'regions.raw';
  LINES      = 100;

  % make debug input data
  PHASE = abs(randn(LINES,2*LINES))*10;
  REGIO = zeros(LINES,2*LINES);
  % make regions
  REGIO(10:50,5:100) = 1;
  REGIO(10:50,125:175) = 3;
  REGIO(70:90,25:185) = 2;
  PHASE = PHASE .* REGIO;
  AMPLI = PHASE .^0.02;
  AMPLI = AMPLI .* REGIO;
  REGIO = uint8(REGIO);%		memory considerations

else
  %%% Handle input (current workspace).
  if (~exist('lines','var'))
    LINES=input('Enter number of lines in binary hgt ampli/phase file: ');
  else
    LINES=lines;
  end
  if (~exist('hgtfile','var'))
    hgtfile    = 'x';
    warning('name of hgtfile is not logged...');
  end
  if (~exist('regionfile','var'))
    regionfile = 'x';
    warning('name of regionfile is not logged...');
  end
  %
  disp('Start reading phase/ampl file (hgt format).');
  disp('Press key to continue');
  pause;
  [PHASE,AMPLI]  = freadhgt(hgtfile,LINES);
  %
  disp('Start reading region file (short int).');
  disp('Press key to continue');
  pause;
  REGIO = freadbk(regionfile,LINES,'int16');
  REGIO = uint8(REGIO);%			memory considerations
end


%%% Plot figures, make gui
RegioFig = figure(3);
AmpliFig = figure(2);
PhaseFig = figure(1);
set(RegioFig,'Name','Regions with unwrapped phase');
set(AmpliFig,'Name','Amplitude');
set(PhaseFig,'Name','Unwrapped phase');
manuwbk_cb('viewphase');
CLimOrig=get(gca,'CLim');

NumCycles  = 1;%			default value 2 add 2 region
NumRegions = max(max(REGIO));%		at start, number of regions (<=255)

%%% Check if logging is requested.
%MFILEFID   = -1;% 			this is assumed in call back sub function
if (exist('LOGMFILE','var'))
  MFILENAME = LOGMFILE;% 		this is assumed in call back sub function
  if (exist(LOGMFILE,'file'))
    warning (['existing file: ',LOGMFILE, ' will be deleted, or exit...']);
    warning ('press key to continue');
    pause
    %warning (['existing file: ',LOGMFILE, ' has been deleted, sorry...']);
  end;
  MFILEFID=fopen(LOGMFILE,'wt');%	text mode
  if (MFILEFID<0)
    error('logmfile not opened ok, permissions?');
  end;
  fprintf(MFILEFID, '%s\n',['%']);
  fprintf(MFILEFID, '%s\n',['% mfile generated by manuwbk.']);
  fprintf(MFILEFID, '%s\n',['% version: ',manuwbkversion]);
  fprintf(MFILEFID, '%s\n',['% date:    ',date]);
  fprintf(MFILEFID, '%s\n',['% for:     ',getenv('USER')]);
  fprintf(MFILEFID, '%s\n',['%']);
  fprintf(MFILEFID, '%s\n',['% Please send your comments to:']);
  fprintf(MFILEFID, '%s\n',['% Bert Kampes, kampes@geo.tudelft.nl']);
  fprintf(MFILEFID, '%s\n',['%']);
  fprintf(MFILEFID, '%s\n',['']);
  %
  % init. globals
  fprintf(MFILEFID, '%s\n',['%%% Initialize variables.']);
  fprintf(MFILEFID, '%s\n',['LINES         = ', num2str(LINES),';']);
  fprintf(MFILEFID, '%s\n',['hgtfile       = ''',hgtfile,'''']);
  fprintf(MFILEFID, '%s\n',['regionfile    = ''',regionfile,'''']);
  fprintf(MFILEFID, '%s\n',['[PHASE,AMPLI] = freadhgt(hgtfile,LINES);']);
  fprintf(MFILEFID, '%s\n',['REGIO         = freadbk(regionfile,LINES,''int16'');']);
  fprintf(MFILEFID, '%s\n',['REGIO         = uint8(REGIO);']);
  fprintf(MFILEFID, '%s\n',['NumRegions    = max(max(REGIO));']);
  fprintf(MFILEFID, '%s\n',['NumCycles     = ', num2str(NumCycles),';']);
  fprintf(MFILEFID, '\n\n');
  % rest is handled in manuwbk_cb functions...
  fclose(MFILEFID);
end

labels = str2mat( ...
  'Too&ls', ...
  '>&AddCycles', ...
  '>>&SetNumber', ...
  '>>&Add', ...
  '>&Profile', ...
  '>>&Horizontal', ...
  '>>&Vertical', ...
  '>>&Improfile', ...
  '>-------', ...
  '>Set ifile dimensions', ...
  '>Load new hgt (a,ph) file', ...
  '>Load new regions file', ...
  '>&Save in hgt file', ...
  '>Save regions file', ...
  '>&Quit^q', ...
  '&View', ...
  '>&Phase', ...
  '>&Amplitude', ...
  '>&Regions', ...
  '>&Spinmap', ...
  '>-------', ...
  '>SetCLimColorbar', ...
  '>SetCLimRegion', ...
  '>RestoreCLim', ...
  '&Regions', ...
  '>&Join', ...
  '>&Split', ...
  '>&Delete', ...
  '&Info^i' ...
);

calls = str2mat( ...
  '', ...
  '', ...
  'manuwbk_cb(''setnumcycles'');', ...
  'manuwbk_cb(''addcycles''); manuwbk_cb(''viewphase'');', ...
  '', ...
  'manuwbk_cb(''profileH'');', ...
  'manuwbk_cb(''profileV'');', ...
  'manuwbk_cb(''profileI'');', ...
  '', ...
  'manuwbk_cb(''getdimbk'');', ...
  '[PHASE,AMPLI]=freadhgt(hgtfile,LINES);manuwbk_cb(''viewphase'');', ...
  'REGIO=uint8(freadbk(''unknown'',LINES,''int16''));manuwbk_cb(''viewregio'');', ...
  'fwritehgt(AMPLI,PHASE);disp(''hgt file written.'')', ...
  'fwritebk(REGIO,''unknown'',''int16'');disp(''region file written.'')', ...
  'clear; close all force; break;', ...
  '', ...
  'manuwbk_cb(''viewphase'');', ...
  'manuwbk_cb(''viewampli'');', ...
  'manuwbk_cb(''viewregio'');', ...
  'figure(PhaseFig); spinmap(1,1);', ...
  '', ...
  'manuwbk_cb(''climcolorbar'');', ...
  'manuwbk_cb(''climregion'');', ...
  'figure(PhaseFig); set(gca,''CLim'',CLimOrig); colorbar;', ...
  '', ...
  'manuwbk_cb(''joinregion'');', ...
  'manuwbk_cb(''splitregion'');', ...
  'manuwbk_cb(''deleteregion'');manuwbk_cb(''viewphase'');', ...
  'manuwbk_cb(''help'');' ...
);

% To avoid multiple redraws, set(gcf,'sharecolors','off') and
%      set(gcf,'renderer','painters').
%  'figure(h); set(gcf,''sharecolors'',''off''); set(gcf,''renderer'',''painters'');', ...
% why doesn't setclim work, if i paste it to workspace it is ok...?

handles = makemenu(PhaseFig, labels, calls);

  
% function dragbar to add 2kpi for region
% function info/help to window
% verbose write m file to reproduce what you did...
% how to return to workspace? get(0)?

%if (MFILEFID>0)
%  fclose(MFILEFID);
%end

%%% EOF
more on
  
  
 
