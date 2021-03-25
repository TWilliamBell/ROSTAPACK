function [Un,Vn,beta,norm_U,norm_V] = normalizeTwoNormUV(U,V)
%   normalizeTwoNormUV:
%       Normalize U and V so that they each have unit norm (in the two
%       norm).
%
%   INPUT:
%       U,V    
%           Two matrices, real or complex, both nonzero.
%   
%   OUTPUT:
%       Un,Vn
%           Scaled versions of U and V such that Un and Vn have unit norm.
%           Furthermore, when U and V are both vectors, then Un*Vn' will
%           also have unit norm.
%
%       beta
%           1/(norm(U)*norm(V)).
% 
%       norm_U,norm_V
%           The norms of the original U and V.
%
%   See also normalizeFroNormUV.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   normalizeTwoNormUV.m introduced in ROSTAPACK Version 1.0
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

    norm_U  = norm(U);
    norm_V  = norm(V);
    Un      = U / norm_U;
    Vn      = V / norm_V;
    beta    = 1 / (norm_U * norm_V);
end
