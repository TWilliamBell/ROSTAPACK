function [U_extrap,V_extrap] = extrapolateRank2Matrix(U,V)
%   extrapolateRank2Matrix:
%       This is a subroutine of extrapolateUpToRank2Matrix.
% 
%       This subroutine performs an implicit matrix extrapolation on a
%       sequence of up to rank-two matrices U_k*V_k' (without ever forming
%       them explicitly), where U_k and V_k are respectively entries in the
%       cell arrays U and V and each U_k*V_k' has unit Frobenius norm.  The
%       user-supplied sequence may contain any mixture of rank-one and
%       rank-two matrices.
%
%       This returns the extrapolated matrix U_extrap*V_extrap', normalized
%       so that it has unit Frobenius norm.
%
%       The matrix extrapolation is done efficiently by choosing up two
%       rows and columns and then performing four individual vector
%       extrapolations on these selections.  A normalized outerproduct
%       U_extrap*V_extrap' is then reconstructed from these four vector
%       extrapolations.
%
%       This method should not be called if the entire sequence consists of
%       rank-one matrices.  
%
%       Note that this method randomly selects which two rows and two
%       columns are selected.  One must set the random number generator
%       seed explicitly in order to have repeatable results.
% 
%   See also extrapolateVectorMPE, normalizeFroNormUV,  
%   normalizeRank1Extrapolation, and extrapolateRank1Matrix.
% 
%
%   For more details, see [Mit14, 6.3.5] and [GGMO17, Section 7.3].     
%       
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   extrapolateRank2Matrix.m introduced in ROSTAPACK Version 1.0
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

    m = size(U{1},1);
    p = size(V{1},1);
    k = length(U);

    % randomly choose two pairs of row and column indices
    row_indices     = chooseTwoIndices(m);
    col_indices     = chooseTwoIndices(p);
    % notational convenience 
    i1              = row_indices(1);
    i2              = row_indices(2);
    j1              = col_indices(1);
    j2              = col_indices(2);

    % Extract rows and cols to extrapolate (all as a column matrices)
    [cols1, cols2]  = extractCols(U,V,col_indices);
    [rows1, rows2]  = extractCols(V,U,row_indices);
     
    % perform four vector extrapolations (match my thesis notation)
    r1              = extrapolateVectorMPE(rows1);  
    r2              = extrapolateVectorMPE(rows2);  
    c1              = extrapolateVectorMPE(cols1);  
    c2              = extrapolateVectorMPE(cols2);  
    
    % the shared entries may not be equal so replace them by their averages
    c1(i1)          = mean(c1(i1) + r1(j1));
    c1(i2)          = mean(c1(i2) + r2(j1));
    c2(i1)          = mean(c2(i1) + r1(j2));
    c2(i2)          = mean(c2(i2) + r2(j2));
    r1(j1)          = c1(i1);
    r2(j1)          = c1(i2);
    r1(j2)          = c2(i1);
    r2(j2)          = c2(i2);
     
    % recover the rank-2 factorization
    u1_extrap       = (c1 + c2) / 2;
    u2_extrap       = (c1 - c2) / 2;
    
    v1_extrap       = recoverVExtrap(u1_extrap, u2_extrap);
    v2_extrap       = recoverVExtrap(u2_extrap, u1_extrap);
         
    U_extrap        = [u1_extrap, u2_extrap];
    V_extrap        = [v1_extrap, v2_extrap];
     
    [U_extrap,V_extrap] = normalizeFroNormUV(U_extrap, V_extrap);
    
    % construct the best rank-1 factorization using the available 
    % extrapolated pairs of rows and columns
    [u_extrap,v_extrap,rank1_error] = getBestRank1Extrapolation();
    
    % use the rank-1 extrapolation if it recovers the extrapolated pair of
    % rows and pair of columns better than the rank-2 factorization
    if rank1_error < getRank2Error()
        % need to add zero second columns to retain dimension 2 of output
        U_extrap = [u_extrap, zeros(m,1)];
        V_extrap = [v_extrap, zeros(p,1)];
    end
    
    % private helper functions
        
    function [cols1, cols2] = extractCols(U,V,indices)
        dim     = size(U{1},1);
        cols1   = zeros(dim,k);
        cols2   = zeros(dim,k);
        
        for j = 1:k
            [u1, u2]    = getColumns(U{j});
            [v1, v2]    = getColumns(V{j});
            cols1(:,j)  = v1(indices(1))*u1 + v2(indices(1))*u2;
            cols2(:,j)  = v1(indices(2))*u1 + v2(indices(2))*u2;
            Uk          = [u1 u2];
            Vk          = [v1 v2];
            Mk          = Uk*Vk.';
            if abs(Mk(:,indices(1)) - cols1(:,j)) + ...
                abs(Mk(:,indices(2)) - cols2(:,j)) > 10*eps
            
            end
        end
    end

    function v_extrap = recoverVExtrap(u1,u2)
        % recall that r1 and r2 are actually in column format
        % so there are no transposes on them
        v_extrap = ...
            ( u2(i2)*r1 - u2(i1)*r2 ) / ...
            ( u1(i1)*u2(i2) - u1(i2)*u2(i1) );
    end

    function [u,v,extrap_error] = getBestRank1Extrapolation()
        u   = cell(1,4);
        v   = cell(1,4);
        err = zeros(1,4);
        
        % compute all four possible rank-1 extrapolations that can be
        % recovered from the two pairs of extrapolated rows and columns
        [u{1},v{1},err(1)]      = getRank1Extrap(true,true);
        [u{2},v{2},err(2)]      = getRank1Extrap(true,false);
        [u{3},v{3},err(3)]      = getRank1Extrap(false,true);
        [u{4},v{4},err(4)]      = getRank1Extrap(false,false);
        
        % return the one which has the least error
        [extrap_error,index]    = min(err);
        u                       = u{index};
        v                       = v{index};
    end

    function [u,v,extrap_err] = getRank1Extrap(first_row,first_col)
        % get a rank-1 extrapolation and its error using only one 
        % of the rows and one of the columns
        
        if first_row
            r           = r1;
            r_index     = i1;
            ro          = r2;
            ro_index    = i2;
        else
            r           = r2;
            r_index     = i2;
            ro          = r1; 
            ro_index    = i1;
        end
        
        if first_col
            c           = c1;
            c_index     = j1;
            co          = c2;
            co_index    = j2; 
        else
            c           = c2;
            c_index     = j2;
            co          = c1;
            co_index    = j1; 
        end
        
        [u,v] = normalizeRank1Extrapolation(r,c,r_index);
        
        % error from chosen row and column
        r_err           = (r - u(r_index)*conj(v)).^2;
        c_err           = (c - conj(v(c_index))*u).^2;
      
        % error from other row and colum
        ro_err          = (ro - u(ro_index)*conj(v)).^2;
        co_err          = (co - conj(v(co_index))*u).^2;

        % remove the shared entries from the columns since these values are
        % already present in r_err and ro_err
        indx            = [r_index ro_index];
        co_err(indx)    = [];
        
        extrap_err = sum(r_err) + sum(ro_err) + sum(c_err) + sum(co_err);     
    end

    function err = getRank2Error()
        row_err = ([r1'; r2'] - U_extrap(row_indices,:) * V_extrap').^2;
        col_err = ([c1, c2] - U_extrap * (V_extrap(col_indices,:))').^2;
        
        % remove the shared entries from the columns since these values are
        % already present in row_err
        col_err(row_indices,:) = [];
        
        err = sum(row_err(:)) + sum(col_err(:));
    end
end

function [u1,u2] = getColumns(U)
    u1 = U(:,1);
    u2 = U(:,2);
end

function indices = chooseTwoIndices(k)
    indices = randi(k,1,2);
    while indices(1) == indices(2)
        indices(2) = randi(k);
    end
end
