function f2 = cpxmultilook(f,rX,rY);
% CPXMULTILOOK  -- complex multilooking of interferogram.
%
%   CPXMULTILOOK(F,rX,rY) Multilook data in f by integer factors rX [rY=rX]
%
%   cpxmultilook(f,r) returns array of dimensions floor(size(f)./[rY rX])
%   with complex multilooked data.
%   defined by ~ f2(i,j) = sum(sum( f((i-1)*rY+1:i*rY, (j-1)*rY+1:j*rX) )) ./ (rX*rY);
%
%   Thus a box is jumping over the data, replacing the first rX,rY data
%   with the mean, etc. Commonly used in insar, but kinda strange.
%
%   Note that cpxmultilook yields different results for magnitude
%   then multilooking of magnitude image only.
%   (complex adds phasors, then computes length.)
%
%   To multilook real arrays one could use the (slow) form:
%     F = CPXMULTILOOK(complex(magnitude),5,5);
%
%   Examples:
%     To multilook an array with a noisy phase trend (10 fringes), magnitude noisy:
%     phase = 10.*2.*pi.*ramp(100,200) + 0.25*pi.*randn(100,200);% 10 fringes
%     mag   = ones(size(phase)) + 0.5.*randn(100,200);
%     cpxdata   = mag.*complex(cos(phase),sin(phase));
%     cpxdata55 = cpxmultilook(cpxdata,5,5);% 
%     figure(101);
%       subplot(1,2,1); imagesc(angle(cpxdata));   title('orig. phase');
%       subplot(1,2,2); imagesc(angle(cpxdata55)); title('multilooked phase');
%
%   See also .
%

% Bert Kampes, 28-Feb-2000
% $Revision: 1.6 $ $Date: 2001/09/28 14:24:31 $

% not true:
%   std::sum used, not sure if this is correct (not mean).

%%% Handle input
if (nargin~=2 & nargin~=3) helphelp; break; end;
if (nargin==2) rY=rX; end;
if (rY*rX==1)
  f2=f;
  return;
end;
isrealf = 0;
if (isreal(f)) isrealf = 1; end;
if (isrealf==1)
  %make complex data, use float values are wrapped ...
  disp('careful, input assumed to be wrapped float phase values.');
  %cosf = cos(f);
  %sinf = sin(f);
  %f    = complex(cosf,sinf);
  f    = complex(cos(f),sin(f));
end



%%% Get dimensions/multilook.
[L P]  = size(f);
f2Y    = floor(L/rY);
f2X    = floor(P/rX);


%%% It seems that with a loop it is fast for large multilook factors.
%%% For smaller, e.g., 2x2, the matrix multiplication is faster.
%%% In matrix multiplication one MUST use the sparse lib to increase speed.
way=2;% matrix is best... if sparse lib used.
if (way==1)
%%% Multilook.
f2     = zeros(f2Y,f2X);
startY = 1;
for ii=1:f2Y
  startX = 1;
  for jj=1:f2X
    f2(ii,jj) = sum(sum(f(startY:ii.*rY,startX:jj*rX)));
    startX = startX+rX;
  end
  startY = startY+rY;
end

else
%%% New way with matrix multiplication.
%%% If input matrix F(6,10) multilook 2,2 then
%%% premultiply by matrix A(F_y/r_y; F_y) (in Y direction, factor=2)
%%% A=[1 1 0 0 0 0;
%%%    0 0 1 1 0 0;
%%%    0 0 0 0 1 1]
%%% Then after multiply by B(F_x; F_x/r_x)
%%% B=[1 0 0 0 0;
%%%    1 0 0 0 0;
%%%    0 1 0 0 0;
%%%    0 1 0 0 0;
%%%    0 0 1 0 0;
%%%    0 0 1 0 0;
%%%    0 0 0 1 0;
%%%    0 0 0 1 0;
%%%    0 0 0 0 1;
%%%    0 0 0 0 1]
%%% Then multilooked = A*F*B; but how to create these A,B smartly?
%%% Like this:

%keyboard
%cputime
if (rY==1)
  A  = 1.;
else
  A  = repmat([ones(1,rY)/rY,zeros(1,L)],1,f2Y).';
  A  = A(1:length(A)-L);
  A  = sparse(reshape(A,L,f2Y).');
end
%cputime
if (rX==1)
  B  = 1.;
else
  B  = repmat([ones(1,rX)/rX,zeros(1,P)],1,f2X).';
  B  = B(1:length(B)-P);
  B  = sparse(reshape(B,P,f2X));
end;
%cputime
f2 = A*f*B;
%cputime
end%method for multilooking

%%% check both ways same result. (checked ok)
%max(max(abs(ff2-f2)))
%max(max(angle(ff2-f2)))


%%% Output.
if (isrealf==1)%		input was real
  f2 = angle(f2);
else%				input was complex
  if (way==1)
  f2 = f2./(rX*rY);%		take mean
  end
end


%%% EOF
