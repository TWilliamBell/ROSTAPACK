function [solver,L,R] = hessSolver(A,E)
%   hessSolver:
%       An "object" for efficiently applying the inverse of zE-A for
%       different values of complex scalars z.  The following upper
%       triangular Hessenberg factorization is formed:
% 
%               (E - A)^{-1} = L(U - H)^{-1}R,
%
%       where U is upper triangular, H is upper Hessenberg, and L and R are
%       unitary matrices.  If E is the identity, then U is the identity and
%       L = R'.  The factorization requires O(n^3) work.
% 
%       However, since zU-H is always upper Hessenberg, (zE-A)^{-1} can now
%       be applied to a vector with only O(n^2) work for any complex scalar
%       z.  For instance, to compute (zE-A)^{-1}v, the computations would
%       proceed as follows:
%       
%       >> Rv = R*v;
%       >> X  = (z*U - H)\Rv;  % but quicker to use linsolve (zU-H is Hess)
%       >> b  = L*X;           % the result
% 
%       The ideas implemented in this routine are due to [1] (when E is the
%       identity) and [2] (otherwise).
%
%   INPUT:
%       A                           [matrix]
%           A square matrix
%       E                           [matrix, optional]
%           A square matrix of the same size as matrix A.  E is taken to be
%           the identity if it is either not provided at all of if it is
%           provided as a zero or empty matrix.
%
%   OUTPUT:         
%       solver
%           An hessSolver "object" with the following methods:
%
%           solver.factor(z)
%               Updates the factorization to zE-A
%
%           [U,H,L,R] = solver.get()
%               Returns the upper triangle Hessenberg factorization for E-A
%               such that L(E - A)R = U - H.
%               NOTE: L and R here are swapped versions of the L and R
%               matrices returned by hessSolver!
%              
%           tf = solver.hasInfNan()
%               Returns true if infs/nans are present in the upper triangle
%               Hessenberg factorization.
%
%           diff = solver.getNormDiff()
%               Checks the quality of the factorization by computing the
%               the relative difference of R'(U-H)L' with respect to E-A.
%           
%           X = solver.applyLeft(B)
%           X = solver.applyLeftCtrans(B)
%           X = solver.applyRight(B)
%           X = solver.applyRightCtrans(B)
%               Respectively computes 
%                   X = (zE-A)\B
%                   X = (zE-A)'\B
%                   X = B/(zE-A)
%                   X = B/(zE-A\)' 
% 
%           X = solver.applyLeftUH(B)
%           X = solver.applyLeftCtransUH(B)
%           X = solver.applyRightUH(B)
%           X = solver.applyRightCtransUH(B)
%               Respectively computes 
%                   X = (zU-H)\B
%                   X = (zU-H)'\B
%                   X = B/(zU-H)
%                   X = B/(zU-H\)' 
% 
%           NOTE 1: all eight of the above methods take advantage of the
%           upper triangular Hessenberg factorization using the current
%           value of z to do the solve in O(n^2) work.
%
%           NOTE 2: solver.factor(z) must have been already called at least
%           once before any of these eight functions can be called.
%
%           NOTE 3: If L and R can be pulled out so that they can be
%           preapplied just once, then UH variants will more efficient for
%           applying (zE-A)^{-1} more than once for a given z.  In this
%           case, the user is expected to manually construct the result
%           they want, i.e., (zE-A)\B, etc., from the output of the applyX
%           functions and the L and R matrices.
%   
%       L,R     
%           Matrices L and R such that (E - A)^{-1} = L(U - H)^{-1}R holds.
%
%   References:
%
%   [1] A. Laub
%       Efficient multivariable frequency response computations. 
%       IEEE Trans. Autom. Control, 26(2):407?408, April 1981.
% 
%   [2] P. Van Dooren and M. Verhaegen  
%       On the use of unitary state-space transformations.  
%       In Linear algebra and its role in systems theory (Brunswick, Maine, 
%       1984), volume 47 of Contemp. Math., pages 447?463.  Amer. Math. 
%       Soc., Providence, RI, 1985.
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%
%   hessSolver.m introduced in ROSTAPACK Version 2.2
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

    if nargin < 2 || isempty(E) || isZero(E) || isEye(E)
        E           = speye(length(A));
        U           = E;
        [L,H]       = hess(A);
        R           = L';          
    else
        [H,U,R,L]   = hess(A,E);
    end
    zUH             = [];
    opts.UHESS      = true;
    opts_CT.UHESS   = true;    
    opts_CT.TRANSA  = true;
    diff            = [];
     
    solver = struct(    'factor',               @factor,                ...
                        'get',                  @get,                   ...           
                        'getNormDiff',          @getNormDiff,           ...
                        'hasInfNan',            @hasInfNan,             ...
                        'applyLeft',            @applyLeft,             ...
                        'applyLeftCtrans',      @applyLeftCtrans,       ...
                        'applyRight',           @applyRight,            ...
                        'applyRightCtrans',     @applyRightCtrans,      ...
                        'applyLeftUH',          @applyLeftUH,           ...
                        'applyLeftCtransUH',    @applyLeftCtransUH,     ...
                        'applyRightUH',         @applyRightUH,          ...
                        'applyRightCtransUH',   @applyRightCtransUH     );

    function factor(z)
        zUH         = z*U - H;      
    end

    function [Uo,Ho,Lo,Ro] = get()
        Uo          = U;
        Ho          = H;
        % Since R*(E - A)*L = U - H, swap R and L matrices so 
        % Lo*(E - A)*Ro = Uo - Ho holds for the user.
        Lo          = R;   
        Ro          = L;  
    end

    function d = getNormDiff(varargin)
        if isempty(diff)
            EmA     = E - A;
            RUmHL   = R'*(U - H)*L';
            diff    = norm(EmA - RUmHL,'fro')/norm(EmA,'fro');
        end
        d           = diff;
    end

    function tf = hasInfNan()
        tf  = ~isFiniteValued(H) || ~isFiniteValued(U);
    end

    function X = applyLeftUH(B)
        X   = linsolve(zUH,B,opts);
    end

    function X = applyRightUH(B)
        X   = (linsolve(zUH,B',opts_CT))';
    end

    function X = applyLeftCtransUH(B)
        X   = linsolve(zUH,B,opts_CT);
    end

    function X = applyRightCtransUH(B)
        X   = (linsolve(zUH,B',opts))';
    end

    function X = applyLeft(B)
        X   = L*(applyLeftUH(R*B));
    end

    function X = applyRight(B)
        X   = (applyRightUH(B*L))*R;
    end

    function X = applyLeftCtrans(B)
        X   = R'*(applyLeftCtransUH(L'*B));
    end

    function X = applyRightCtrans(B)
        X   = (applyRightCtransUH(B*R'))*L';
    end

end