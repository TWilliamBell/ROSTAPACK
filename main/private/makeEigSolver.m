function solver = makeEigSolver(A,B,C,D,E,epsilon,discrete_time)
%   makeEigSolver:
%       For a given spectral value set, specified by matrices A,B,C,D,E and
%       epsilon >= 0, this method finds the farthest outward boundary point
%       along a line determined by y.
% 
%   INPUT:
%       A,B,C,D,E                   [matrices]
%           System matrices.  A must be provided while B,C,D,E can be given
%           explicitly or as [] for their shortcuts.
% 
%       epsilon                     [nonnegative real scalar]
%           Perturbation level of the spectral value set
%
%       discrete_time               [logical]
%           False:  find (purely imaginary) eigenvalues that correspond to
%                   boundary points along vertical lines
%           True:   find (unimodular) eigenvalues that correspond to
%                   boundary point along circles centered at the origin
%   
%   OUTPUT:
%       solver                      [function handle]
%           where d = solver(r) and
%               r                   [real-valued scalar]
%               d                   [column vector of complex values]
%           Returns all eigenvalues d of the pencil, where the:
%           - imaginary ones correspond to boundary points along the
%             vertical line given by x = r                  (continuous)
%           - unimodular ones correspond to boundary points along the 
%             circle centered at the origin with radius r   (discrete)
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   makeEigSolver.m introduced in ROSTAPACK Version 2.0
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
       
    % For fixed epsilon, it is most efficient to preconstruct as much of 
    % the pencil as possible and then do minor updates for each vertical
    % line or circle requested, rather constructing the pencil from scratch
    % each time.
    
    n       = length(A);
    E_ident = isempty(E) || isEye(E);
    if E_ident && ~issparse(E) 
        % Ensure E is a speye for efficiency
        E   = speye(n);
    end
   
    idx1    = 1:n;
    idx2    = n+1:2*n;
    
    if discrete_time
        [M,N]           = initDiscrete(A,B,C,D,epsilon);
        updatePencil    = @updateRadius;
        solver          = @getEigsMN;
        return
    end   
   
    % Continuous time pencils - two possibilities E = I or not
    updatePencil        = @updateShift;
    
    if E_ident
        % In this case, the pencil is just a regular eigenvalue problem
        % eig(H), which can be much cheaper than calling eig(H,I).
        [M,A_BRinvDhC]  = initContinuousM(A,B,C,D,epsilon);
        solver          = @getEigsM;
        return
    end
    
    [M,N,A_BRinvDhC]    = initContinuousMN(A,B,C,D,E,epsilon);
    solver              = @getEigsMN;
    
    % Private functions
       
    function d = getEigsM(shift)  
        updatePencil(shift);
        d               = eig(M);
        d(isinf(d))     = [];
    end

    function d = getEigsMN(shift_radius)        
        updatePencil(shift_radius);
        d               = eig(M,N);
        d(isinf(d))     = [];
    end

    function updateShift(shift)
        A_sE            = A_BRinvDhC - shift*E;
        M(idx1,idx1)    = A_sE;
        M(idx2,idx2)    = -A_sE';
    end

    function updateRadius(radius)
        rE              = radius*E;
        M(idx2,idx2)    = rE';
        N(idx1,idx1)    = rE;
    end
end

function [M,A_BRinvDhC] = initContinuousM(A,B,C,D,epsilon)
    [A_BRinvDhC,gBRinvBh,gChSinvC,n] = getComponents(A,B,C,D,epsilon);
    Z   = sparse(n,n);
    M   = toFull([ Z, -gBRinvBh; gChSinvC, Z]);
end

function [M,N,A_BRinvDhC] = initContinuousMN(A,B,C,D,E,epsilon)
    [A_BRinvDhC,gBRinvBh,gChSinvC,n] = getComponents(A,B,C,D,epsilon);
    Z   = sparse(n,n);
    M   = toFull([ Z, -gBRinvBh; gChSinvC, Z]);
    N   = toFull([ E, Z; Z, E' ]);
end

function [M,N] = initDiscrete(A,B,C,D,epsilon)
    [A_BRinvDhC,gBRinvBh,gChSinvC,n] = getComponents(A,B,C,D,epsilon);
    Z   = sparse(n,2*n);
    M   = toFull([ A_BRinvDhC, -gBRinvBh; Z ]);
    N   = toFull([ Z; -gChSinvC, A_BRinvDhC' ]);
end

function [A_BRinvDhC,gBRinvBh,gChSinvC,n] = getComponents(A,B,C,D,epsilon)
    [n,m,p]     = systemDimensions(A,B,C,D);
    gamma       = 1/epsilon;
    if isempty(D) || isZero(D)      
        % Doing matmuls with identities is more expensive than checking
        % whether or not they are identities and constructing the result
        Z       = sparse(n,n);
        if isempty(B) || isEye(B)
            BBh         = Z;
            indx        = 1:n+1:n*m;
            BBh(indx)   = 1;   
        else
            BBh         = B*B';
        end
        if isempty(C) || isEye(C)
            ChC         = Z;
            indx        = 1:n+1:p*n;
            ChC(indx)   = 1;
        else
            ChC         = C'*C;
        end         
        A_BRinvDhC  = A;
        gBRinvBh    = (-1/gamma)*BBh;
        gChSinvC    = (-1/gamma)*ChC;
    else
        if isempty(B)
            B       = speye(n,m);
        end
        if isempty(C)
            C       = speye(p,n);
        end
        gamma2      = gamma^2;
        R           = (D'*D) - gamma2*speye(m);
        S           = (D*D') - gamma2*speye(p);
        [L,U]       = lu(R);
        BRinv       = (B / U) / L;
        gBRinvBh    = gamma*(BRinv*B');
        gChSinvC    = gamma*(C'*(S \ C));
        A_BRinvDhC  = A - (BRinv*D'*C);
    end
end
