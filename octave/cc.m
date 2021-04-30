% CC clears all variables and closes figures

% $Revision: 1.4 $  $Date: 2001/09/28 14:24:30 $
% Bert Kampes, 1/03/00

%%% Remove really everything.  Reset more, warning, etc.
echo off;
close all force;
clear all;
dbclear all;
warning debug;
%dbstop if error;
more on;
clc;

%%% Be nice.
fprintf(1,'\a\nCleared all. Good Luck!\n\n');

%%% EOF
