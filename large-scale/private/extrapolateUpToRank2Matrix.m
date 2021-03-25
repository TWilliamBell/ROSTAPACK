function [U_extrap,V_extrap] = extrapolateUpToRank2Matrix(U,V)
%   extrapolateUpToRank2Matrix:
%       This is a subroutine of extrapolateObject.
%   
%       Performs an implicit matrix extrapolation on a sequence of up to
%       rank-two matrices U_k*V_k' (without ever forming them explicitly),
%       where U_k and V_k are respectively entries in the cell arrays U and
%       V and each U_k*V_k' has unit Frobenius norm.  The user-supplied
%       sequence may contain any mixture of rank-one and rank-two matrices.
%
%       This returns the extrapolated matrix U_extrap*V_extrap', normalized
%       so that it has unit Frobenius norm.
%
%       Based on the dimensions of the provided matrix sequences, this
%       routine will automatically determine the most appropriate
%       subroutine:
%           - extrapolateColumnMatrix
%           - extrapolateRank1Matrix
%           - extrapolateRank2Matrix.
%
%       Note that when extrapolateRank2Matrix is employed, one must set the
%       random number generator seed explicitly in order to have repeatable
%       results, as it uses random initialization.
% 
%   See also extrapolateVectorMPE, extrapolateColumnMatrix,
%   extrapolateRank1Matrix, extrapolateRank2Matrix.  
%
% 
%   For more details, see [Mit14, Sections 4.2 and 6.3.5], 
%   [MO16, Section 6] and [GGMO17, Section 7.3].
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   extrapolateUpToRank2Matrix.m introduced in ROSTAPACK Version 1.0
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

    [m,r]   = size(U{1});
    p       = size(V{1},1);
    k       = length(U);
    
    % NOTE: getVectors() just extracts the first column if U or V is a cell
    % array containing two column wide matrices
    
    if p == 1           % V just contains unimodular scalars
        [U_extrap,V_extrap] = extrapolateColumnMatrix(getVectors(U,m,k),...
                                                      getVectors(V,p,k) );                                         
    elseif m == 1       % U just contains unimodular scalars
        [V_extrap,U_extrap] = extrapolateColumnMatrix(getVectors(V,p,k),...
                                                      getVectors(U,m,k) );   
    elseif r == 1 || allRank1(U) || allRank1(V)
        % implicit rank 1 extrapolation (using first columns of U and V)
        [U_extrap,V_extrap] = extrapolateRank1Matrix( getVectors(U,m,k),...
                                                      getVectors(V,p,k) );
    else 
        % implicit rank 2 extrapolation 
        [U_extrap,V_extrap] = extrapolateRank2Matrix(U,V);
    end
    
    % resize U_extrap and V_extrap so that they have the same number of 
    % columns as the input
    if r == 2 && size(U_extrap,2) == 1
        U_extrap = [U_extrap, zeros(m,1)];
        V_extrap = [V_extrap, zeros(p,1)];
    end        
end

function tf = allRank1(U)
    tf = true;
    for j = 1:length(U)
        U_j = U{j};
        [~,cols] = size(U_j);
        if cols > 1 && nnz(U_j(:,2)) > 0
            tf = false;
            return
        end
    end
    return
end

function vectors = getVectors(celldata,dim,k)
    vectors = zeros(dim,k);
    for j = 1:k
        vectors(:,j) = celldata{j}(:,1);
    end
end
