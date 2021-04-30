function pra4 (handle);
% PRA4  print figure maximum enlarged to A4, while maintainig aspect ratio.
%    PRA4 prints the current figure to the default printer
%    with maximum dimensions, without distorting the figure.
%
%    PRA4(H) does the same for figure with handle H.
%
%    PRA4('MAX') prints the figure with maximum dimensions to fit
%    on a4 paper (no need to want this for printing to file).
%
%    PRA4(OPTIONS) enlarges the figure and executes PRINT with the
%    options in the string OPTIONS.
%
%    To print the current figure to a file named file.eps, in epsc format,
%    enlarged to fill an a4 paper, while maintaining the width/height ratio use:
%    PRA4('-depsc file.eps');
%
%    To print to a printer with the name colorprinter use (?):
%    PRA4('-f1 -depsc1 file.eps');
%    PRA4('-f1 -Pcolorprinter');
%    figure(1); PRA4('-Pcolorprinter');
%
%    The printing options can be set in file PRINTOPT.M in the directory
%    of environment variable MATLABPATH.
%    In my case [cmd,dev]=printopt returns:
%      cmd = 'lp -d5simx4de -oA4 -odpi6 -onb'
%      dev = '-dps2'
%
%    See also ORIENT, PRINT, PRINTOPT, GET, SET.
%

% $Revision: 1.8 $  $Date: 2001/09/28 14:24:32 $
% initial function
%  Rens Swart, 17-Mar-2000
% printing to printer instead of file, extended options, simplified.
%  Bert Kampes, 29-Mar-2000


%%% Handle input.
options='';
maintainaspectratio = 1;
if (nargin==0) handle = gcf; end;
if (nargin==1)
  if (~ishandle(handle))
    options = deblank(handle);
    q1      = findstr(options,'-f');
    if (isempty(q1))% no -f specified in print command
      handle = gcf;
    else
      q=options(q1:length(options));%    e.g., '-f 1 file.eps' or '-f1 file.eps'
      spaces = find(isspace(q));
      if (isempty(spaces))
	handle  = str2num(q(3:length(q)));
	options = options(1:q1-1);
      else
	for ii=1:length(spaces)
	  handle = str2num(q(3:spaces(ii)));
	  if (ishandle(handle))
	    q2 = q1+spaces(ii);
            options = options([1:q1-1,q2:length(options)]);
	    break;
	  end
	end
      end
    end
  end;
end;
if (~ishandle(handle)) error('pra4: no handle.'); end;
if (length(options)==3)
  if (lower(options)=='max')
    options='';
    maintainaspectratio=0;
  end;
end;

%%% set options 2b sure (and leave it like this at exit...).
set(handle,'papertype','a4letter');
set(handle,'paperunits','centimeters');
origpaperorientation=get(handle,'paperorientation');
orig=get(handle,'units');
set(handle,'units','centimeters');
figuresize=get(handle,'position');
%set(handle,'units','pixels');%			set back to avoid problems...
set(handle,'units',orig);%			set back to avoid problems...

%%% scale up figure for printer;
left     = figuresize(1);
right    = figuresize(2);
width    = figuresize(3);
height   = figuresize(4);
if (width>height)
  set(handle,'paperorientation','landscape');
else
  set(handle,'paperorientation','portrait');
end;
or=get(handle,'paperorientation');
a4papersize = get(handle,'papersize');
a4width  = a4papersize(1);
a4height = a4papersize(2);


offset = 0;%		keep free at edges [cm]
if (maintainaspectratio==1)
  a4width  = a4width -2*offset;%		keep off
  a4height = a4height-2*offset;
  scaleh   = a4width  ./ width;
  scalev   = a4height ./ height;
  scale    = min(scaleh,scalev);
  if (scale<1) scale=max(scaleh,scalev); end;%		prevent stupid situation.
  width    = width*scale; 
  height   = height*scale; 
  left     = offset + (a4width -width) ./2.;
  right    = offset + (a4height-height)./2.;
  maxpos = [left, right, width, height];
else
  maxpos = [offset, offset, a4width-offset, a4height-offset];
end;

%%% Some info for user.
disp(['  figure size:          ', num2str(figuresize(3)),' ', num2str(figuresize(4))]);
disp(['  printing size:        ', num2str(maxpos(3)),' ', num2str(maxpos(4))]);
disp(['  printing orientation: ', or]);
if (maintainaspectratio)
  disp('  keeping aspect ratio, no distortion of figure.');
else
  disp('  Distortion of figure.');
end;
disp(['  command: print -f', num2str(handle),' ', options]);



origpaperpos=get(handle,'paperposition');%
set(handle,'paperposition',maxpos);%			where it is all about
eval(['print -f', num2str(handle),' ', options]);%	or is it
set(handle,'paperposition',origpaperpos);%		set back
set(handle,'paperorientation',origpaperorientation);%	set back

%%% Note: figure object properties
%get(1)
%PaperPosition = [0.634517 6.34517 20.3046 15.2284]
%PaperSize = [20.984 29.6774]
%PaperType = a4letter


%%% EOF
