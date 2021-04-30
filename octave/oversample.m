function f2 = oversample(f,factorrow,factorcol);
% OVERSAMPLE  Oversample (harmonic interpolate) data.
%    OVERSAMPLE(F)  Oversample F harmonically by factor 2 in both
%    dimensions. If F is a vector only in 1 dimension.
%    The dimensions of F have to be powers of two.
%    Harmonical interpolation/oversampling is done by zero padding
%    in the spectral domain.
%
%    OVERSAMPLE(F,FACTOR)  Oversample F harmonically by integer FACTOR 
%    in both directions (if appropriate, vector in one direction only.).
%
%    OVERSAMPLE(F,FACTORROW,FACTORCOL)  Oversample F harmonically by
%    integer FACTORROW in row direction, and by factor FACTORCOL in
%    second, column direction.
%
%    OVERSAMPLE(F,F1,F2) returns array of dimensions SIZE(F).*[F1 F2]
%    with harmonically oversampled data.
%
%
%    Note: extrapolated at end...
%    Note: Implementation can be improved.
%
%    See also SIZE, FFT2, IFFT2,
%    seems to do the same as INTERPFT (matlab)?
%

% $Revision: 1.6 $ $Date: 2000/04/28 10:53:25 $
% Bert Kampes, 24-Feb-2000


%%% Handle input
if (nargin<1) error('OVERSAMPLE: nargin=0'); end;
if (nargin<2) factorrow = 2; end;%			default
if (nargin<3) factorcol = factorrow; end;%		default

if (~isint(factorrow)) error('OVERSAMPLE: only integer factors.'); end;
if (~isint(factorcol)) error('OVERSAMPLE: only integer factors.'); end;
if (factorrow<1) error('OVERSAMPLE: factorrow < 1.'); end;
if (factorcol<1) error('OVERSAMPLE: factorcol < 1.'); end;
if (factorrow*factorcol==1) error('OVERSAMPLE: factors<=1'); end;
[L P] = size(f);
if (L*P<=1) error ('OVERSAMPLE: data is scalar'); end;
if (~ispow2(L) & factorrow~=1) error('number of rows should be power of 2.'); end
if (~ispow2(P) & factorcol~=1) error('number of columns should be power of 2.'); end

realinput = 0;
lying     = 0;
standing  = 0;
if (isreal(f)) realinput = 1; end;
if (L==1)      lying     = 1; end;
if (P==1)      standing  = 1; end;

%%% Vector in one direction, perhaps better seperate routine.
if (standing)
  factorcol=1;
  if (factorrow==1) error('OVERSAMPLE: standing vector, factor==1'); end;
end;
if (lying)
  if (nargin==1) factorrow=1; end;%			factorcol default
  if (nargin==2) factorcol=factorrow; factorrow=1; end;
  if (factorcol==1) error('OVERSAMPLE: lying vector, factor==1'); end;
  if (factorrow~=1) error('OVERSAMPLE: lying vector, factorrow~=1'); end;
end


%%% Start oversampling
F  = fft2(f);
clear f;

%%% Last term has to be divided by 2 'cause of even L,P
L2 = L*factorrow;%			dimension of result
P2 = P*factorcol;
ll = floor(L/2);%			ok
pp = floor(P/2);
F2 = zeros(L2,P2);

%if (factorrow ~= 1)
if (factorrow ~= 1 & ~lying)
  F(ll+1,:) = F(ll+1,:)./2.;
end;
%if (factorcol ~= 1)
if (factorcol ~= 1 & ~standing)
  F(:,pp+1) = F(:,pp+1)./2.;
end;

%%% Actual zero padding at correct place
if (factorrow == 1)
  F2(       :,           1:pp+1) = F(    :,        1:pp+1);
  F2(       :,     P2-pp+1:P2  ) = F(    :,     pp+1:P   );
elseif (factorcol == 1)
  F2(      1:ll+1,        :)     = F(   1:ll+1,     :    );
  F2(L2-ll+1:L2,          :)     = F(ll+1:L,        :    );

else%	oversampling in both dimensions
  F2(      1:ll+1,       1:pp+1) = F(   1:ll+1,    1:pp+1);
  F2(      1:ll+1, P2-pp+1:P2  ) = F(   1:ll+1, pp+1:P   );
  F2(L2-ll+1:L2,         1:pp+1) = F(ll+1:L,       1:pp+1);
  F2(L2-ll+1:L2,   P2-pp+1:P2  ) = F(ll+1:L,    pp+1:P   );
end
clear F;


%%% Inverse transformation
f2 = ifft2(F2);
clear F2;


%%% Correct for ...
if (realinput==1) f2 = real(f2);    end;

%f2 = f2(1:L2-r+1,1:P2-r+1);
f2 = factorrow*factorcol*f2;%			correct for size

%%%EOF
