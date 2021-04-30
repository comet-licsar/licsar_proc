function enlargefig(factor)
% ENLARGEFIG -- make figure larger for printing/viewing.
%  
%   This function changes the 'position' and 'paperposition'
%   properties of the current figure by a factor 2.
%
%   ENLARGEFIG(F) uses F as factor.
%
%   See also PRA4
%

% $Revision: 1.2 $  $Date: 2001/09/28 14:24:32 $
% RH, 1999?
%// BK 07-May-2001

%not ok, bk: set(gcf,'paperposition',get(gcf,'paperposition')*factor)
if (nargin==0)
  factor = 2;% default
end

pos       = get(gcf,'position');
ppos      = get(gcf,'paperposition');
pos(1:2)  = pos(1:2).*0;
pos(3:4)  = pos(3:4).*factor;
ppos(3:4) = ppos(3:4).*factor;
set(gcf,                 ...
   'position',      pos, ...
   'paperposition', ppos);

