function [x,y] = getExtremalPoint(A,B,C,D,E,epsilon,x0,y,opts)
%   getExtremalPoint:
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
%       x0,y                        [real-valued scalar pair]
%           For fixed y and x0, do a search for an x along either:
%           - horizontal line x + 1i*y with x > x0          (continuous)
%           - rotated line x*exp(1i*y) with |x| > x0        (discrete)
%           to find the farthest boundary point.  In the continuous-time
%           case, this means to the right.  In the discrete-time case, this
%           means outward in either direction.
%    
%       opts                        [struct]
%       .discrete_time              [logical]
%           Select whether to do search along a horizontal line (false) or
%           a rotated line through the origin (true).
% 
%       .ham_symp_tol               [nonnegative real scalar]
%           Tolerance for determining whether a complex eigenvalue value is
%           deemed purely imaginary, i.e. one that corresponds to a
%           boundary point of the spectral value set along the search line.
%   
%   OUTPUT:
%       x,y                         [real-valued scalar pair]
%           Continuous time:
%           - The boundary point x + 1i*y with x > x0
%           Discrete time:  
%           - The boundary point x*exp(1i*y) with |x| > x0.  
%           - If this point is actually in the opposite direction of the 
%             angle given by y, then y is updated to y + pi to reflect
%             this.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   getExtremalPoint.m introduced in ROSTAPACK Version 2.0
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

    discrete_time   = opts.discrete_time;
    ham_symp_tol    = opts.ham_symp_tol;

    if discrete_time
        theta       = pi/2 - y;
        tau         = polar2cplx(1,theta);
        x           = 0;
        getPointsFn = @(r) polar2cplx(r,y);
        xValsFn     = @abs;
        imagOrderFn = @(d) abs(imag(d));
    else
        tau         = 1i;
        x           = -y;
        getPointsFn = @(x) x + 1i*y;
        xValsFn     = @real;
        imagOrderFn = @(d) imag(d);
    end    
   
    % Compute all the eigenvalues of the rotated matrix pencil
    eigSolver       = makeEigSolver(tau*A,tau*B,C,D,E,epsilon,false);
    eHSP            = eigSolver(x);
    
    % Get all the (nearly) imaginary eigenvalues
    indx            = abs(real(eHSP)) <= ham_symp_tol;  
    eHSP            = eHSP(indx);
    
    % Delete all eigenvalues that are either:
    %   - to the left of the current abscissa approximation (continuous) 
    %   - inside the current radius approximation (discrete)
    indx            = imagOrderFn(eHSP) <= x0;
    eHSP(indx)      = [];
    
    if isempty(eHSP)
        x           = -inf;
        return
    end
    
    % Sort them by largest imaginary part first 
    [~,indx]        = sort(imagOrderFn(eHSP),'descend');
    eHSP            = eHSP(indx);
    
    % Select outermost one on this ray
    eHSP            = eHSP(1);
    point           = getPointsFn(imag(eHSP));
    x               = xValsFn(point);
    
    % Extremal point was actually in the opposite direction on this line so
    % update the angle to reflect this.
    if discrete_time && imag(eHSP) < 0
        y           = wrapToPi(y + pi);
    end

end