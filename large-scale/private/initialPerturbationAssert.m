function pert = initialPerturbationAssert(sys,real_frobenius_norm,pert)
%   initialPerturbationAssert:
%       This routine does system and stability measure specific assertions
%       on the provided perturbation epsilon*U*V' and eigenvalue selection.
%       
%       Note this method does not assert properties which could have been
%       checked regardless of the system.  For examples, it is assumed that
%       pert.kth_eigenvalue is a positive integer, pert.epsilon is a
%       positive finite real value, and that pert.U and pert.V have at most
%       2 columns.  For non-specific asserts, see
%       initialPerturbationOptions.
%          
%       Finally, provided that all asserts succeed, this will output the
%       processed and verified perturbation and eigenvalue selection to be
%       used with subsequent routines.  Note that the processed unit-norm
%       perturbation U*V' may be slightly different than the user-supplied
%       U*V', since this returns the normalized version computed when
%       asserting that it is indeed unit norm.  Furthermore, when
%       real_frobenius_norm is true, U and V will each be padded with a
%       second column of zeros, if the provided U and V are column vectors.
%        
%   INPUT:
%       sys                     [required]
%           A stateSpaceABCD object.
% 
%       real_frobenius_norm     [required logical]
%           True if the spectral value sets are real-valued and bounded by
%           the Frobenius norm.  False for complex-valued spectral value
%           sets.
%
%       pert                    [required struct]
%           This struct provides the perturbation epsilon*U*V' and desired
%           eigenvalue selection.  It must have the following fields.
%               .kth_eigenvalue
%               .epsilon 
%               .U 
%               .V
% 
%   THROWS:
%       This method throws an error for the provided perturbation
%       epsilon*U*V' if:
%           - kth_eigenvalue > length(A) (from system sys)
%           - epsilon >= 1/norm(D) (only relevant for nonzero D matrices)
%           - U*V' is not dimensionally compatible with matrices B and C
%             (from system sys)
%           - U*V' does not have unit norm (in the chosen norm)
%           - U and V are not column vectors 
%             (only relevant when real_frobenius_norm is false)
%           - U and V have complex-valued entries
%             (only relevant when rela_frobenius_norm is true)
%   
%   OUTPUT:       
%       pert                    [struct]
%           Processed and verified version of the user-supplied eigenvalue
%           selection and initial perturbation epsilon*U*V', where U*V' has
%           unit norm:  
%               .kth_eigenvalue 
%               .epsilon 
%               .U
%               .V
%           
%   See also initialPerturbationOptions and stateSpaceABCD.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   initialPerturbationAssert.m introduced in ROSTAPACK Version 1.0
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

    NORM_MATCH  = 10*eps;
    
    [n,m,p]     = sys.getDimensions();
    
    k           = pert.kth_eigenvalue;
    epsilon     = pert.epsilon;
    U           = pert.U;
    V           = pert.V;
    
    assert(k <= n, dimensionErrorMsg(n));  
    
    if ~isempty(epsilon)
        limit   = epsilonUpperBound(sys);
        assert(epsilon < limit, epsilonErrorMsg(limit));
    end

    if isempty(U) % v must also be empty since it has already been checked
        return
    end
    
    [U_rows,U_cols] = size(U);
    [V_rows,V_cols] = size(V);
    
    assert(U_rows == m, uvDimensionErrorMsg('U',m));
    assert(V_rows == p, uvDimensionErrorMsg('V',p));

    if real_frobenius_norm
        assert(isRealValued(U),rfsvsRealValuedErrorMsg('U'));
        assert(isRealValued(V),rfsvsRealValuedErrorMsg('V'));
        if U_cols == 1
            U = [U zeros(m,1)];
            V = [V zeros(p,1)];
        end
        [pert.U,pert.V,beta] = normalizeFroNormUV(U,V);   
        norm_UV = 1/beta;  
        assert(abs(norm_UV - 1) <  NORM_MATCH, rfsvsNormErrorMsg());
    else
        assert(U_cols < 2,csvsColumnsErrorMsg('U')); 
        assert(V_cols < 2,csvsColumnsErrorMsg('V'));
        norm_U = norm(U);
        norm_V = norm(V);
        assert(abs(norm_U - 1) < NORM_MATCH, csvsNormErrorMsg('U'));
        assert(abs(norm_V - 1) < NORM_MATCH, csvsNormErrorMsg('V'));
        pert.U = U/norm_U;
        pert.V = V/norm_V;
    end   
end

function m = dimensionErrorMsg(dim)
    s = 'For this system, .initial_perturbation.kth_eigenvalue must be at most %d!';
    m = sprintf(s,dim);
end

function m = uvDimensionErrorMsg(name,dim)
    s = 'For this system, .initial_perturbation.%s must have %d rows!';
    m = sprintf(s,name,dim);
end

function m = epsilonErrorMsg(limit)
    s = 'For this system, .initial_perturbation.epsilon must be less than %g!';
    m = sprintf(s,limit);
end

function m = csvsNormErrorMsg(name)
    h = 'When real_frobenius_norm is false, .initial_perturbation.';
    m = 'must have have unit norm!';
    m = sprintf('%s%s %s',h,name,m);
end

function m = csvsColumnsErrorMsg(name)
    h = 'When real_frobenius_norm is false, .initial_perturbation.';
    m = 'must have at most 1 column!';
    m = sprintf('%s%s %s',h,name,m);
end

function m = rfsvsNormErrorMsg()
m = [
'When real_frobenius_norm is true, .initial_perturbation.U and .V '     ...
'must be normalized so that norm(u*v.'',''fro'') is equal to one.'      ];
end

function m = rfsvsRealValuedErrorMsg(name)
    h = 'When real_frobenius_norm is true, .initial_perturbation.';
    m = 'must only contain real values!';
    m = sprintf('%s%s %s',h,name,m);
end
