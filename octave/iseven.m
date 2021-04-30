function res = iseven(in)
% ISEVEN  --  True for even numbers.
%
%   ISEVEN(IN) returns 1 for the elements in IN that are even,
%   0 otherwise.
%
%   To test if any/all element is even in a array, use e.g., any or all:
%     if(any(iseven(1:4))) disp('at least one element is even'); end;
%     if(all(iseven(1:4))) disp('all elements are even'); end;
%
%   Example:
%     a=[1 2 3 4 5];
%     iseven(a)% returns [0 1 0 1 0]
%
%   See also ISODD, ISINT, REM, MOD.
%

% $Revision: 1.2 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 31-Mar-2000

% not ok, array may be mixed odd/even
%res=~isodd(in);

%res=0;
%if(~rem(in,2))res=1;end;

res=(~rem(in,2));

%%% EOF
