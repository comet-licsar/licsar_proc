function trap(h,varargin);
% TRAP -- Position figures as trap (escalator).
%
%    Reset position and order of all figures so that they
%    lie shifted with respect to eachother. 
%
%    TRAP(H) does the same for handles in H.
%    TRAP(H1,H2,...) does the same for specified handles.
%
%    Example:
%      TRAP         sorts all figures of the main window.
%      TRAP(10:100) only sorts the figure windows with handles between 10 and 100.
%    
%    See also TOP, TIP
%    

% Idea of Rens Swart * 22 maart 2000
%// Bert Kampes, 11-Dec-2000 insar toolbox, changed implementation
%// $Revision: 1.4 $  $Date: 2001/09/28 14:24:34 $

switch nargin
  case 0% Get the children of main win, sort it, change it.
    qch = flipud(sort(get(0,'children')));
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
if (isempty(qch))
  warning('No figure windows found, exiting.');
  break;
end
%
set(0,'DefaultFigureUnits','pixels')
qnwpos = get(0, 'DefaultFigurePosition');
for qnum = 1:length(qch)
  qhandle = qch(qnum);
  if ishandle (qhandle)
    set(qhandle,'Position',qnwpos);
    qnwpos(1:2) = qnwpos(1:2) + [15 -15];
    figure(qhandle);
  end
end


%%% EOF

