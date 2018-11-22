% pop_readedf() - load a EDF EEG file (pop out window if no arguments).
%
% Usage:
%   >> [dat] = pop_readedf( filename);
%
% Inputs:
%   filename       - EDF file name
% 
% Outputs:
%   dat            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 13 March 2002
%
% See also: eeglab(), readedf()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 13 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: pop_readedf.m,v $
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, command] = pop_readedf(filename); 
command = '';

if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.EDF', 'Choose an EDF file -- pop_readedf()'); 
	if filename == 0 return; end;
	filename = [filepath filename];
end;

% load datas
% ----------
EEG = eeg_emptyset;
[EEG.data, header] = readedf(filename);  

EEG.filename        = filename;
EEG.filepath        = '';
EEG.setname 		= 'EDF file';
EEG.nbchan          = header.channels;
EEG.srate           = header.samplerate;
EEG.xmin            = 0; 

EEG = eeg_checkset(EEG);
command = sprintf('EEG = pop_readedf(''%s'');', filename); 

return;
