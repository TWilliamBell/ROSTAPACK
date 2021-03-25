function [brkFn,rootFn,countFn] = makeInsideFns(A,B,C,D,E,epsilon,opts)
%   makeInsideFns:
%       For a given spectral value set, specified by matrices A,B,C,D,E and
%       epsilon >= 0, this method returns two function handles that each
%       compute f(x,y), the norm of the transfer function minus epsilon at
%       a point given by x and y (Cartesian or polar coordinates).  These
%       two functions also return its first and optionally second partial
%       derivatives with respect to x.  A point is:
%           - inside the spectral value set if f(x,y) > 0
%           - on the boundary if f(x,y) = 0
%           - outside the spectral value set if f(x,y) < 0
%       If B=C=I, D=0, and n=m=p, by default the functions will
%       equivalently compute the reciprocal of the smallest singular value
%       of zE-A, instead of ||(zE-A)^{-1}||, since the latter is more
%       expensive due to the inverse.  However, this can be overriden.
% 
%   INPUT:
%       A,B,C,D,E                   [matrices]
%           System matrices.  A must be provided while B,C,D,E can be given
%           explicitly or as [] for their shortcuts.
% 
%       epsilon                     [nonnegative real scalar]
%           Perturbation level of the spectral value set
%
%       opts                        [struct]
%       .discrete_time              [logical]
%           False:  x and y are Cartesian coordinates, z = x + 1i*y
%           True:   x and y are polar coordiantes, z = x*exp(1i*y)
%
%       .bracket_order              [1 or 2]
%           1:      brkFn will only compute the first derivative
%           2:      brkFn will compute both first and second derivatives
%
%       .root_order                 [real value in [1,2]]
%           < 2:    rootFn will only compute the first derivative
%           2:      rootFn will compute both first and second derivatives
%
%       .solve_type                 [0,1,2]
%           Sets how inverses of zE-A should be applied:
%           0: via upper triangular Hessenberg factorization 
%           1: via LU
%           2: via LU with permutations 
% 
%       .fast_search                [logical]
%           If this false, brkFn will only provide the function value while
%           rootFn will provide the function value and its first derivative
%           (but not its second), regardless of the settings of
%           opts.bracket_order and opts.root_order.
%
%       .force_two_norm             [logical]
%           When B=C=I, D=0, and n=m=p, force ||(zE-A)^{-1}|| to be used 
%           instead of the reciprocal of the smallest singular value of
%           zE-A.  (Since GESDD in LAPACK can have numerical difficulty for
%           very poorly scaled matrices.)
%   
%   OUTPUT:
%       brkFn                       [function handle]
%       rootFn                      [function handle]
%           where [f,g,h] = brkFn(x,y) and [f,g,h] = rootFn(x,y) where
%               f:  the norm of the transfer function minus epsilon
%               g:  the first derivative with respect to x
%               h:  the second derivative with respect to x 
%                   (nan if it was not asked to be computed)
%           at the point specified by x and y
% 
%       countFn                     [function]
%           Calling countFn() returns the total number of times either 
%           brkFn or rootFn has been called.
%
%   See also makeSigmaFn and makeSigmaMinFn.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   makeInsideFns.m introduced in ROSTAPACK Version 2.0
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

    brk_order           = opts.bracket_order;
    root_order          = opts.root_order;
    
    [ntf0,ntf1,ntf2]    = makeNtfFns(A,B,C,D,E,opts);
  
    if ~opts.fast_search 
        brkFn   = @(e,l) ntfRoot(ntf0,e,l);   
        rootFn  = @(e,l) ntfRoot(ntf1,e,l);   
    elseif max(brk_order,root_order) < 2
        brkFn   = @(e,l) ntfRoot(ntf1,e,l);
        rootFn  = brkFn;
    elseif min(brk_order,root_order) == 2
        brkFn   = @(e,l) ntfRoot(ntf2,e,l);
        rootFn  = brkFn;
    elseif root_order == 2
        brkFn   = @(e,l) ntfRoot(ntf1,e,l);
        rootFn  = @(e,l) ntfRoot(ntf2,e,l);
    else
        brkFn   = @(e,l) ntfRoot(ntf2,e,l);
        rootFn  = @(e,l) ntfRoot(ntf1,e,l);
    end
  
    target_gain = 1/epsilon;
    count       = 0;
    countFn     = @getCount;

    function [f,g,h] = ntfRoot(ntf,eta,loc)
        [f,g,h] = ntf(eta,loc,1);
        f       = f - target_gain;
        count   = count + 1;
    end

    function n_evals = getCount()
        n_evals = count;
    end
end