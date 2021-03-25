function [solve_fn,update_fn] = phiSystemsBackSlash(D_obj)
%   phiSystemsBackSlash:
%       This is a subroutine of phiSystemsSolver.
%
%       This routine solves the following pair of linear systems:
%       
%           Phi_m = (I - epsilon^2 D' * D)
%           Phi_p = (I - epsilon^2 D * D'),
%
%       which arise when approximating the complex-valued spectral value
%       set abscissa and/or radius.
%
%       This routine uses a simple backlash to do the solves but this is
%       generally not the most efficient way to solve these systems.
%       
%   INPUT:
%       D_obj
%           A matrixObject representation.  If D is not an explicit full or
%           sparse matrix, a full representation will be created.
%   
%   OUTPUT:
%       [b,c] = solve_fn(By,Cx,varargin)
%           Given right hand sides By and Cx, this returns:
%               b = Phi_m^(-1) By
%               c = Phi_p^(-1) Cx.
%
%       update_fn(new_epsilon)  
%           Calling this will update the value of epsilon with the user's
%           new value.  Note that this method does NOT check if the new
%           value violates epsilon's upper bound: epsilon < 1/norm(D).
%
%   See also phiSystemsFactorizeBoth, phiSystemsFactorizeLowRank,
%   phiSystemsFactorizeSmaller, phiSystemsPCG, phiSystemsSolver, and
%   phiSystemsSolverOptions.
%
%
%   For more details, see [GGO13, Section 3], [Mit14, Section 4.1] and 
%   [MO16, Section 6].
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   phiSystemsBackSlash.m introduced in ROSTAPACK Version 1.0
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

    epsilon_sq  = 0;
    [p,m]       = D_obj.getSize();      
    D           = D_obj.formFull();
    
    % Explicitly form Phi_p and Phi_m and use backslashes to solve
    DDh         = D_obj.apply(D');
    DhD         = D_obj.applyHermitian(D);
    Im          = speye(m);
    Ip          = speye(p);
    Phi_m       = Im;
    Phi_p       = Ip;
    
    update_fn   = @updateEpsilon;
    solve_fn    = @solveSystems;
    
    function updateEpsilon(epsilon) 
        new_epsilon_sq  = epsilon^2;
        if new_epsilon_sq ~= epsilon_sq
            epsilon_sq  = new_epsilon_sq;
            Phi_m       = Im - epsilon_sq*(DhD);
            Phi_p       = Ip - epsilon_sq*(DDh);
        end
    end
    
    function [b,c] = solveSystems(By,Cx,varargin)
        b = Phi_m\By;
        c = Phi_p\Cx;
    end
end
