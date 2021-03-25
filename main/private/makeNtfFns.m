function [fn0,fn1,fn2] = makeNtfFns(A,B,C,D,E,opts)
%   makeNtfFns:
%       This returns three function handles for evaluating the norm of the
%       transfer function and its first and second derivative(s).
% 
%   INPUT:
%       A,B,C,D,E                   [matrices]
%           System matrices.  A must be provided while B,C,D,E can be given
%           explicitly or as [] for their shortcuts.
%
%       opts                        [struct]
%       .discrete_time              [logical]
%           False:  x and y are Cartesian coordinates, z = x + 1i*y
%           True:   x and y are polar coordiantes, z = x*exp(1i*y)
%
%       .solve_type                 [0,1,2]
%           Sets how inverses of zE-A should be applied:
%           0: via upper triangular Hessenberg factorization 
%           1: via LU
%           2: via LU with permutations 
% 
%       .force_two_norm             [logical]
%           When B=C=I, D=0, and n=m=p, force ||(zE-A)^{-1}|| to be used 
%           instead of the reciprocal of the smallest singular value of
%           zE-A.  (Since GESDD in LAPACK can have numerical difficulty for
%           very poorly scaled matrices.)
%   
%   OUTPUT:
%       [fn0,fn1,fn2]               [function handles]
%           fn0: only computes the norm of the transfer function
%           fn1: also computes the first derivative(s)
%           fn2: also computes the first and second derivative(s)
%
%   See also makeSigmaFn and makeSigmaMinFn.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   makeNtfFns.m introduced in ROSTAPACK Version 2.2
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

    polar       = opts.discrete_time;
    solve_type  = opts.solve_type; 
    [n,m,p]     = systemDimensions(A,B,C,D);
    
    if (m == n && (isempty(B) || isEye(B))) && ...
       (p == n && (isempty(C) || isEye(C))) && ...
       (isempty(D) || isZero(D))
        if opts.force_two_norm
            [fn0,fn1,fn2] = makeSigmaFn(A,B,C,D,E,polar,solve_type);
        else
            [fn0,fn1,fn2] = makeSigmaMinFn(A,E,polar);
        end
    else
        [fn0,fn1,fn2] = makeSigmaFn(A,B,C,D,E,polar,solve_type);
    end
end