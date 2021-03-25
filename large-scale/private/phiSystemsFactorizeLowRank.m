function [solve_fn,update_fn] = phiSystemsFactorizeLowRank(D_obj)
%   phiSystemsFactorizeLowRank:
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
%       This routine solves each system by forming two LU factorizations of
%       two (presumably small) k by k systems, where k is the shared
%       dimension of Dl and Dr and assumed to be smaller than either m or
%       p, where [p,m] = size(Dl*Dr').  Thus, by exploiting this low-rank
%       structure, this routine can solve both linear systems faster than
%       phiSystemsFactorizeBoth and can efficiently scale to
%       large-dimensional D matrices, where both m and p are large.
%       
%   INPUT:
%       D_obj
%           A matrixObject representation.  D must have been specified as a
%           low-rank outer product D = Dl*Dr'.
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
%   See also phiSystemsBackSlash, phiSystemsFactorizeBoth,
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
%   phiSystemsFactorizeLowRank.m introduced in ROSTAPACK Version 1.0
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
    [~,~,k]     = D_obj.getSize();
    [Dl,Dr]     = D_obj.getOuterProduct();
  
    % D is given as D = Dl*Dr', low rank-k factorization
    % Two rank-k LU factorizations for Phi_m and Phi_p
    Ik          = speye(k);  
    L_m         = Ik;
    U_m         = Ik;
    L_p         = Ik;
    U_p         = Ik;
    
    DlDl        = Dl'*Dl;
    DrDr        = Dr'*Dr;
    DlDrDr      = Dl*DrDr;
    DrDlDl      = Dr*DlDl;
    DlDlDrDr    = DlDl*DrDr;
    DrDrDlDl    = DlDlDrDr';
    
    update_fn   = @updateEpsilon;
    solve_fn    = @solveSystems;
     
    function updateEpsilon(epsilon) 
        new_epsilon_sq  = epsilon^2;
        if new_epsilon_sq ~= epsilon_sq
            epsilon_sq  = new_epsilon_sq;
            [L_m,U_m]   = lu(Ik - epsilon_sq*DrDrDlDl);
            [L_p,U_p]   = lu(Ik - epsilon_sq*DlDlDrDr);
        end
    end

    function [b,c] = solveSystems(By,Cx,varargin)
        b = By + epsilon_sq*DrDlDl*(U_m\(L_m\(Dr'*By)));
        c = Cx + epsilon_sq*DlDrDr*(U_p\(L_p\(Dl'*Cx)));
    end
end
