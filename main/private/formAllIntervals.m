function ints = formAllIntervals(eHSP,opts)
%   formAllIntervals:
%       Given a set of complex values, this method takes the subset of
%       those that are (nearly) purely imaginary or (nearly) unimodular and
%       forms all nonoverlapping adjacent intervals using their frequencies
%       (imaginary part or angle, respectively) as endpoints.  The set of
%       intervals is sorted in increasing order by these endpoints and each
%       interval has positive length.
% 
%   INPUT:
%       eHSP                [vector of complex-valued scalars]
%           From the subset of these that are (nearly) purely imaginary or
%           unimodular, a vector of frequencies (imaginary parts or
%           angles), in increasing order and named freqs, will be
%           constructed and used to construct the intervals.
%               
%       opts                [struct of parameters]
%       .discrete_time      [logical]
%           False:  purely imaginary values, imaginary part frequencies
%           True:   unimodular values, angle frequencies, and the
%                   "wrap-around" interval will be included, i.e. the
%                   interval that includes pi, either by touching or
%                   crossing the negative x-axis.
%
%       .ham_symp_tol       [nonnegative real scalar]
%           Tolerance for determining whether a complex value is deemed
%           purely imaginary or unimodular.
%
%       .is_symmetric       [logical]
%           When this is true:
%   `       - all intervals with negative endpoints (both) are discarded
%           - the first interval will include 0, either as an endpoint or a
%             as a midpoint
%           - if opts.discrete_time is true, the last interval will include
%             pi, either as an endpoint or as a midpoint.
%  
%   OUTPUT:
%       ints                [cell array of intervals]
%           Each entry in the cell array is a two element vector giving the
%           endpoints of each interval.  The endpoints are in increasing
%           order and the cell array of intervals is also in increasing
%           order, also with respect to the endpoints.  If no intervals
%           could be constructed (e.g. not enough unique frequencies), then 
%           {} is returned.
%
%   See also splitAnInterval.
%             
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   formAllIntevals.m introduced in ROSTAPACK Version 2.0
%
% =========================================================================
% |  ROSTAPACK: RObust STAbility PACKage                                  |
% |  Copyright (C) 2014-2019 Tim Mitchell                                 |
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
    
    discrete        = opts.discrete_time;
    ham_symp_tol    = opts.ham_symp_tol;
    is_symmetric    = opts.is_symmetric;
    
    % First extract and sort the frequencies to use as interval endpoints
    if discrete
        % The (nearly) unimodular eigenvalues and their frequencies
        indx    = abs(abs(eHSP)-1) <= ham_symp_tol;
        freqs   = angle(eHSP);
    else
        % The (nearly) purely imaginary eigenvalues and their frequencies
        indx    = abs(real(eHSP)) <= ham_symp_tol;
        freqs   = imag(eHSP);
    end
    if is_symmetric
        % If problem is symmetric, only need the nonnegative frequencies
        indx    = indx & (freqs >= 0);
    end
    freqs       = freqs(indx);
    if isempty(freqs)
        ints    = {};
        return
    end
    % unique is slow so just sort - we will handle any duplicates below
    freqs       = sort(freqs,'ascend');
 
    % Second, form all the nonoverlapping adjacent intervals.
    % But, we first may need to add a frequency before and/or after, to
    % ensure 0 and possibly pi will be in the first and last intervals.
    f0      = [];
    fend    = [];
    if is_symmetric     
        % Interval that contains frequency 0 must be included, since it may
        % be a cross section 
        if freqs(1) > 0
            f0      = -freqs(1);
        end
        % For discrete-time problems, the interval that contains frequency
        % pi must all be included, since it too may be a cross section
        if discrete && freqs(end) < pi
            fend    = -freqs(end) + 2*pi;
        end       
    else    
        % Must include the wrap-around interval for discrete-time problems
        % but different than above because now pi may not be the midpoint
        if discrete % add wrap-around interval
            fend    = freqs(1) + 2*pi;
        end
    end
    freqs       = [f0; freqs; fend];
    % Need at least two frequencies to form intervals
    if numel(freqs) < 2
        ints    = {};
        return
    end
    ints        = makePairs(freqs);  
    % Finally, remove any interval of zero length
    idx         = cellfun(@(x) x(1) == x(2),ints);
    ints(idx)   = [];

end