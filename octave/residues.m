function [nres, residumatrix] = residues(phasematrix,oldmethod);
% RESIDUES  --  Find residues for bert.
%
%   NRES = RESIDUES(MAT) returns number of residues NRES in MAT. If MAT is 
%   real it is assumed to contain wrapped phase values, if it is complex
%   then first the angle of MAT is computed.
%
%   [NRES, RMAT] = RESIDUES(MAT) optionally returns a matrix with
%   the locations of the residues {-2pi / 0 / +2pi}.
%
%   Residues are related to upperleft corner.
%
%   Residues are defined by .... see ghiglia&pritt 98, bamler&hartl 99 etc.
%
%   [NRES, RMAT] = RESIDUES(MAT,dummy) uses for loops to compute the residues.
%   Since matlab hates loops, it is nice to see the algorithm with diff is
%   much faster.
%
% See also ANGLE, WRAP, CMULTILOOK, DIFF.
%

% $Revision: 1.6 $  $Date: 2001/09/28 14:24:32 $
% Bert Kampes, 1/03/00

%%% Handle input
if (nargin==0 | nargin>2) helphelp; break; end;
if (~isreal(phasematrix)) phasematrix=angle(phasematrix); end;


if (nargin==2)% use oldmethod
nres = 0;
[L P] = size(phasematrix);
if (nargout==2)
  residumatrix=zeros(L,P);
  for ii=1:L-1
    for jj=1:P-1
      residumatrix(ii,jj) = ...
              wrap(phasematrix(ii+1,jj  ) - phasematrix(ii  ,jj  )) + ...
              wrap(phasematrix(ii+1,jj+1) - phasematrix(ii+1,jj  )) + ...
              wrap(phasematrix(ii  ,jj+1) - phasematrix(ii+1,jj+1)) + ...
              wrap(phasematrix(ii  ,jj  ) - phasematrix(ii  ,jj+1));
    end
  end
  % Count in residumatrix
  xxx = find(abs(residumatrix)<1);
  residumatrix(xxx)=0;
  nres = L*P-length(xxx);

else
  for ii=1:L-1
    for jj=1:P-1
      residue = wrap(phasematrix(ii+1,jj  ) - phasematrix(ii  ,jj  )) + ...
                wrap(phasematrix(ii+1,jj+1) - phasematrix(ii+1,jj  )) + ...
                wrap(phasematrix(ii  ,jj+1) - phasematrix(ii+1,jj+1)) + ...
                wrap(phasematrix(ii  ,jj  ) - phasematrix(ii  ,jj+1));
      if (abs(residue) > 1)
        nres = nres + 1;
      end 
    end
  end
end


%%% new method using diffs...
else

  %-------------------------------
  % this looks ok, but still one for loop...
  %-------------------------------
  %residumatrix=zeros(L,P);
  %  for ii=1:L-1
  %    phaseline = phasematrix(ii:ii+1,:);
  %    diffrow   = wrap(diff(phaseline,1,1));
  %    diffcol   = wrap(diff(phaseline,1,2));
  %    % note wrap(x) = -wrap(-x);
  %%%rubbish...
  %%%diffcol(2,:)     = -diffcol(2,:)
  %%%diffrow(2:2:L-1) = -diffrow(2:2:L-1)
  %diff1 = diff(diffrow,1,2)
  %diff2 = diff(diffcol,1,1)
  %diff2-diff1
  %residumatrix(ii,1:P-1) = diff2-diff1;
  %  end

  %-------------------------------
  %%% Compute residue matrix the fast (matlab) way.
  %note bk: wrap(x) = -wrap(-x);
  [L P] = size(phasematrix);
  diffrows = wrap(diff(phasematrix,1,1));
  diffcols = wrap(diff(phasematrix,1,2));
  clear phasematrix;

  %residumatrix=zeros(L,P);
  %diffrowscols = diff(diffrows,1,2);
  %diffcolsrows = diff(diffcols,1,1);
  %residumatrix(1:L-1,1:P-1) = diffcolsrows-diffrowscols;
  residumatrix = diff(diffcols,1,1) - diff(diffrows,1,2);
  %-------------------------------

  %%% Count in residumatrix.
  xxx = find(abs(residumatrix)<1);% non residues
  %nres = L*P-length(xxx);
  nres = (L-1)*(P-1)-length(xxx);
  if (nargout==2)% correct size
    residumatrix(xxx) = 0;%		not necessary (?).
    residumatrix(L,:) = 0;
    residumatrix(:,P) = 0;
  end;
end;

%%% EOF
