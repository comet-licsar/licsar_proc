function top(h,varargin);
% TOP -- Position figures on top of eachother.
%
%    Reset position and order of all figures so that they
%    are alligned on top of eachother. (use <TAB> key to view)
%
%    TOP(H) does the same for handles in H.
%    TOP(H1,H2,...) does the same for specified handles.
%
%    Example:
%      TOP         sorts all figures of the main window.
%      TOP(10:100) only sorts the figure windows with handles between 10 and 100.
%
%    See also TRAP, TIP 
%    

%// $Revision: 1.4 $  $Date: 2001/09/28 14:24:33 $
%// Bert Kampes, 11-Dec-2000 insar toolbox, changed implementation

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
    qnwpos2 = get(qhandle, 'Position');
    qnwpos2(1:2) = qnwpos(1:2);
    set(qhandle,'Position',qnwpos2);
    figure(qhandle);
  end
end

%%% EOF

