function [Un,Vn,beta,norm_U,norm_V] = normalizeFroNormUV(U,V)
%   normalizeFroNormUV:
%       Normalize outer product U*V' so that it that has unit Frobenius
%       norm.  This routine use the trace inner product of the formulation
%       of the Frobenius norm in order to be more efficient than forming
%       U*V' explicitly and then computing norm(U*V','fro').
%
%   INPUT:
%       U,V    
%           Two matrices, real or complex, both nonzero.
%   
%   OUTPUT:
%       Un,Vn
%           Scaled versions of U and V such that Un*Vn' has unit Frobenius
%           norm.
%
%       beta
%           1/norm(Un*Vn','fro').
% 
%       norm_U,norm_V
%           The norms of the original U and V.
%
%   See also normalizeTwoNormUV.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   normalizeFroNormUV.m introduced in ROSTAPACK Version 1.0
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

    % first calculate the norm since beta = 1/norm(U*V','fro') is needed
    norm_fro_UV = frobeniusNormFast(U,V);
    beta        = 1/norm_fro_UV;

    % then renormalize U and V separately to try to make their entries
    % of similar magnitude
    norm_U      = norm(U);
    norm_V      = norm(V);
    U           = U / norm_U;
    V           = V / norm_V;
    % note that the Frobenius norm will likely have a different value 
    % than above since U and V have been renormalized separately
    scalar      = frobeniusNormFast(U,V)^0.5;
    Un          = U/scalar;
    Vn          = V/scalar;
end
