function stringOut = spacesToUnderscores(stringIn);
%
% spacesToUnderscores will replace all spaces with underscores within a string
%
% Cai Wingfield 1-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

if ~isempty(stringIn)
	stringOut = strrep(stringIn, ' ', '_');
end

end%function
