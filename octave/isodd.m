function res = isodd(in)
% ISODD  --  True for odd numbers.
%   ISODD(IN) returns 1 for the elements in IN that are odd,
%   0 otherwise.
%
%   To test if any/all elements in an array are odd, use any or all:
%     if(any(isodd(1:4))) disp('at least one element is odd'); end;
%     if(all(isodd(1:4))) disp('all elements are odd'); end;
%
%   Example:
%     a=[0 1 2 3 4];
%     isodd(a)% returns [0 1 0 1 0]
%
%   See also ISEVEN, ISINT, REM, MOD.
%

% $Revision: 1.2 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 31-Mar-2000

%res=0;
%if(rem(in,2))res=1;end;

res=(rem(in,2));

%%% EOF
