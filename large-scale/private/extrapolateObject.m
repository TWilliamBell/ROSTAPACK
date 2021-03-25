function obj = extrapolateObject(max_size,min_size,rollover_size)
%   extrapolateObject:
%       This object performs extrapolation on sequences of rank 1 or 
%       rank 2 matrices, where each matrix has unit norm.  The following
%       two types of sequences are supported:
%       
%           1) Only rank 1 matrices, either real or complex, and of unit
%              spectral norm.
% 
%           2) At most rank 2 matrices, all of unit Frobenius norm.  The
%              sequence can have both rank 1 and rank 2 matrices in it.
%
%   INPUT:
%       max_size                [positive integer >= 2]   
%           The number of sequence entries to use for an extrapolation
% 
%       min_size                [positive integer >= 2]
%           The minimum number of sequence entries necessary to attempt an
%           extrapolation of, when an early extrapolation is requested by
%           the caller.
% 
%       rollover_size           [nonnegative integer]
%           Number of the most recent accrued sequence entries that will be
%           reused for the next extrapolation attempt.  When this is 0,
%           performing an extrapolation will also clear all the stored
%           sequence entries.
%   
%   OUTPUT:
%       obj 
%           An "object", a struct containing the following functions for
%           performing extrapolations of sequences rank 1 or 2 matrices:
%
%       [U_ext,V_ext,k_used] = obj.updateAndExtrapolate(U,V,try_early) 
%           Adds U*V' to the internally maintained sequence of entries and,
%           if the extrapolation parameters are satisfied, returns the
%           extrapolation U_ext*V_ext', along with how many entries were
%           used to construct this extrapolation.  Note that if try_early
%           is true, an early extrapolation will be attempted, provided
%           that there are at least min_size entries.  If extrapolation was
%           not performed, U_ext and V_ext will each be returned as [] and
%           k_used will be zero.
%
%       .startNewSequence(U,V)
%           Clear the internal sequence of entries and start a new sequence
%           with U*V' as the first entry.
%
%   See also extrapolateVectorMPE.
% 
%
%   For more details, see [Mit14, Sections 4.2 and 6.3.5], 
%   [MO16, Section 6] and [GGMO17, Section 7.3].
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   extrapolateObject.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  ROSTAPACK: RObust STAbility PACKage                                  |
% |  Copyright (C) 2014-2018 Tim Mitchell                                 |
% |                                                                       |
% |  This file is part of ROSTAPACK                                       |
% |                                                                       |
% |  ROSTAPACK is free software: you can redistribute it and/or modify    |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  ROSTAPACK is distributed in the hope that it will be useful,         |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/agpl.html>.                             |
% =========================================================================

    U_seq   = cell(1,max_size);
    V_seq   = cell(1,max_size); 
    count   = 0;
  
    obj     = struct(   'update',               @updateAndExtrapolate,  ...
                        'startNewSequence',     @startNewSequence       );

    function [U_ext,V_ext,k_used] = updateAndExtrapolate(U,V,try_early)
        
        if nargin < 3
            try_early = false;
        end
    
        if count < max_size
            count           = count + 1;
            U_seq{count}    = U;
            V_seq{count}    = V;
        end
        
        [U_ext,V_ext,k_used] = extrapolate(try_early);
        
        % clear some of the older entries to make room for new ones
        if count == max_size
            % count may be smaller than rollover_size
            k_to_rollover           = min(rollover_size,count);
            start_index             = count - k_to_rollover + 1;
            indices                 = start_index:count;
            
            U_seq(1:k_to_rollover)  = U_seq(indices);
            V_seq(1:k_to_rollover)  = V_seq(indices);
            count                   = k_to_rollover;
        end
    end

    function startNewSequence(U,V)
        count       = 1;
        U_seq{1}    = U;
        V_seq{1}    = V;
    end

    % private helper function

    function [U_ext,V_ext,k_used] = extrapolate(try_early)
        U_ext   = [];
        V_ext   = [];
        k_used  = 0;
        if count == max_size || (try_early && count >= min_size)
            k_used  = count;
            try
                [U_ext,V_ext] = extrapolateUpToRank2Matrix(U_seq,V_seq);
            catch
                return
            end

            % if extrapolation is numerically invalid, return []
            if ~isFiniteValued(U_ext) || ~isFiniteValued(V_ext)
                U_ext   = [];
                V_ext   = [];
            end
        end
    end 
end
