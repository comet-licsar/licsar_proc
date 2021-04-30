function res = ispow2(in)
% ISPOW2  True for powers of two.
%   ISPOW2(IN) returns 1 if all elements of IN are powers of two,
%   and 0 otherwise.
%   See also NEXTPOW2, POW2.
%

% $Revision: 1.2 $  $Date: 2000/03/16 09:38:31 $
% Bert Kampes, 15/03/00

res = 1;
if (prod(size(in))==1)
  if (2^nextpow2(in)~=in) res=0; end;

else% check all individual elements the clumsy way.
  for ii=1:size(in,1)
    for jj=1:size(in,2)
      if (2^nextpow2(in(ii,jj))~=in(ii,jj)) res=0; return; end;
    end;
  end;
end;

%%% EOF
