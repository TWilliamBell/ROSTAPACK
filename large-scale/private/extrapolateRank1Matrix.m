function [u,v] = extrapolateRank1Matrix(U,V)
%   extrapolateRank1Matrix:
%       This is a subroutine of extrapolateUpToRank2Matrix.
% 
%       This subroutine performs an implicit matrix extrapolation on a
%       sequence of rank one matrices u_k*v_k' (without ever forming them
%       explicitly), where u_k and v_k are respectively columns U and V and
%       each vector u_k and v_k has unit norm.
%
%       This returns the extrapolated matrix u*v', given as vectors u and
%       v, each normalized to have unit norm.
%
%       The matrix extrapolation is done efficiently by choosing a
%       particular and column and then performing two individual vector
%       extrapolations on this row and column of the matrix sequence.  A
%       normalized outerproduct u*v' is then reconstructed
%       from these two vector extrapolations.
% 
%   See also extrapolateVectorMPE and normalizeRank1Extrapolation.
% 
%
%   For more details, see [Mit14, Sections 4.2 and 6.3.5] and
%   [MO16, Section 6].
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   extrapolateRank1Matrix.m introduced in ROSTAPACK Version 1.0
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
    
    % since the perturbation is u*v', we need to conjugate the vectors in V
    V = conj(V);
     
    [row_u,row_index,u_nz] = findDensestRow(U);
    [row_v,col_index,v_nz] = findDensestRow(V);
    
    k = size(U,2);
        
    if u_nz < k || v_nz < k
        warning('Vector Extrapolation: not all entries are nonzero!');
    end
    
    % extrapolate a single row and column from dense matrix u*v'
    extrap_col  = extrapolateVectorMPE(U*diag(row_v));   
    extrap_row  = extrapolateVectorMPE(V*diag(row_u));
    
    % shared entry between the extrapolated row and col may not entirely 
    % agree so average them and replace the values with the average
    entry_col               = extrap_col(row_index);
    entry_row               = extrap_row(col_index);
    entry_avg               = 0.5 * (entry_col + entry_row);
    
    extrap_col(row_index)   = entry_avg;
    extrap_row(col_index)   = entry_avg;
     
    [u,v] = normalizeRank1Extrapolation(extrap_row, extrap_col, row_index);
   
end

function [row,indx,n_nz] = findDensestRow(A)    
    [n_nz,indx] = max(sum((A ~= 0),2));  
    row = A(indx,:);
end
