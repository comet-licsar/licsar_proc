function hh = titlenormal(h,string)
% TITLEBOLD -- Set title properties to normal face.
%    TITLEBOLD by itself set the title of the current figure to normal.
%    (or all children?)
%    TITLEBOLD(h) where h is a vector of handles, set all these to normal.
%    TITLEBOLD('text') where text is a character string, add this title
%    to the current figure and makes it normal.
%    TITLEBOLD(h,'text') where h is a vector of handles and text contains
%    the same number of title strings.
%    TITLEBOLD('text',h) does the same as TITLEBOLD(h,t).
%    g = TITLEBOLD('text') returns the handle(s) to the text object.
%    
%    See also TITLE, XLABEL, YLABEL, SET, GET
%    

%    WARNING: this function set the first child only.
%    thus first use title, then colorbar, then titlenormal.
%    not colorbar, title, titlenormal
%    

%// $Revision: 1.2 $  $Date: 2001/03/16 13:47:09 $
%// Bert Kampes, 18-Dec-2000

%%% Check input.
handles = [];
t       = [];
if (nargin==0)
  handles = gcf;
elseif (nargin==1)
  if (ishandle(h))
    handles=h;
  elseif (ischar(h))
    handles=gcf;
    t=h;
  end;
elseif (nargin==2)
  % get handles in 'handles'
  if (ishandle(h))
    handles = h;
  elseif (ishandle(string))
    handles = string;
  else
    handles = gcf;
    warning('handle not recognized, using gcf.');
  end;
  % get title in 't'
  if (ischar(h))
    t = h;
  elseif (ischar(string))
    t = string;
  else
    warning('text string not recognized, using no text.');
  end;
else 
  warning('i think number of input arguments is wrong, continuing.');
end
if (isempty(handles))
  error('No figure handles found.');
end

%%% Check input.
%if (length(handles)) ~= size(t,1)
%  error('size handles not equal to size titles');
%end


%%% Set to normal for input handles.
for ii = 1:length(handles)
  % assume do forall axes, (no colorbar? etc. / how to check for istitle?)
  % this also does it for title of colorbar!
  child = get(handles(ii),'children');
  child = child(1);%			hmmm...
  h=get(child,'title');
  if (~isempty(t)) set(h,'string',t); end;
  set(h,'fontweight','normal');
end

if (nargout==1) hh=h; end;

%%% EOF

