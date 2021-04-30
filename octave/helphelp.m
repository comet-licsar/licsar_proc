function [] = helphelp(callingfile);
%  HELPHELP  --  general utility to give help if wrong.
% 
%  HELPHELP returns options menu for user.
%  HELPHELP('filename');
%  Example: in function that has wrong input arguments, call helphelp"
%    if (wronginput) helphelp; break; end;
%  or use the equivalent:
%    if (wronginput) helphelp(mfilename); break; end;
%
%  This utility sets more on, gives help of the callingfile, and
%  returns to callingfuntion without returning a value.
%

%// BK 10-Sep-2001
%// $Revision: 1.1 $  $Date: 2001/09/28 14:24:43 $

%%% Obtain calling file name(s).
disp(['---------------------']);
disp(['Input error detected:']);
disp(['---------------------']);
if (nargin==0)
  k = dbstack;
  switch length(k)
    case 1
      disp(['* Function "',mfilename, ...
	    '" called from mainwin, giving help on "',mfilename,'".']);
      callingfile = k(1).name;
    case 2
      callingfile = k(2).name;
    otherwise
      callingfile = k(2).name;
      for ii=length(k):-1:2
        disp(['* ', blanks(length(k)-ii+1), k(ii).name, ' --->']);
      end
  end
end
directory = [];
mfile     = callingfile;
[directory, mfile, extention, v] = fileparts(callingfile);
%for ii=length(callingfile):-1:1
%  if (callingfile(ii)=='/')%  UNIX dir seperation...
%    directory = callingfile(1:ii);
%    mfile = callingfile(ii+1:length(callingfile));
%    break;
%  end
%end
%disp(['Input error in file: ', mfile, '  (',directory,'/)']);
%disp(['directory:           ', mfile]);
%disp([callingfile, ': Input error']);
%disp('Press key for help, "q" to break.');

%%%
if (isempty(directory))
  disp(['* File: ', mfile]);
else
  disp(['* File: ', mfile, '  [', directory, ']']);
end

%%%
q=input(['!! Press key for help on ',mfile,';  "q" to break. '],'s');
if (length(q)~=1 | q~='q')
  %pause;% has advantage no enter required.
  more on;
  help(callingfile);
end

%return;
%error('dd');
%break;

%%% EOF.

