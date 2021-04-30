function tip(h,varargin)
% TIP -- Position figures to cover the whole screen.
%
%    Reset position and order of all figures so that they
%    cover the whole screen.
%
%    TIP(H) does the same for figure handles in H.
%    TIP(H1,H2,...) does the same for specified handles.
%    
%    Example:
%      TIP         sorts all figures of the main window.
%      TIP(10:100) only sorts the figure windows with handles between 10 and 100.
%    
%    See also TOP, TRAP
%    

%// $Revision: 1.4 $  $Date: 2001/09/28 14:24:33 $
%// Bert Kampes, 11-Dec-2000 insar toolbox, changed implementation

%%% Get the children of main win, sort it.
switch nargin
  case 0
    qch = sort(get(0,'children'));
  otherwise
    qch = lying(h(ishandle(h)));
    for q = 1:length(varargin)
      h   = varargin{q};
      h   = lying(h(ishandle(h)));
      qch = [qch,h];
    end
end
%qch = unique(qch);
%%%
qnumfigs = length(qch);
if (qnumfigs==0)
  warning('No figure windows found, exiting.');
  break;
elseif (qnumfigs>10)
  warning('Too many figure windows, exiting.');
  break;
end

%%% Prepare reshape display.
set(0,'DefaultFigureUnits','pixels')
qsc       = (get(0,'screensize'));
qscwidth  = qsc(3);
qscheight = qsc(4);
qnY       = 1;
if (qnumfigs>=4) qnY = 2; end;
qnX       = ceil(qnumfigs/qnY);
qsizeX    = floor(qscwidth/ qnX);
qsizeY    = floor(qscheight/qnY);
%
qfig = 0;
for qY = 1:qnY
  for qX = 1:qnX
    qpos    = [qsc(1)+(qX-1)*qsizeX, ...
	      qsc(2)+rem(qY*qsizeY,qscheight), ...
	      qsizeX, floor(0.85*qsizeY)];
    qfig    = qfig + 1;
    if (qfig>qnumfigs) break; end;
    qhandle = qch(qfig);
    if ishandle (qhandle)
      set(qhandle,'Position',qpos);
      figure(qhandle);
    end
  end
end

%%% EOF

