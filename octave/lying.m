function out = lying(in);
% LYING  lying vector (vectorize)
%   LYING(MAT) returns a row vector. If MAT is a matrix it is 
%   vectorized as MAT(:) would.
%   See also LYING, PUNCT, COLON
%

% $Revision: 1.3 $  $Date: 2001/03/16 13:46:44 $
% Bert Kampes, 1/03/00

out = in(:).';

%%% EOF
