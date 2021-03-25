function opts = initialPerturbationOptions(user_opts)
%   initialPerturbationOptions:
%       Validate user options struct for specifying an (optional) initial
%       perturbation (full or partially given) and eigenvalue selection.
%       If user_opts is [] or not provided, returned opts will be the
%       default parameters.
%       
%   USAGE:
%       opts = initialPerturbationOptions();
%       opts = initialPerturbationOptions(user_opts);
%
%   INPUT:
%       user_opts   Struct of settable parameters.  No fields are required,
%                   irrelevant fields are ignored, and user_opts may be 
%                   given as [].
%   
%   OUTPUT:
%       opts        Struct of all validated parameters.  If a field is
%                   provided in user_opts, then the user's value is checked
%                   to whether or not it is a valid value, and if so, it is
%                   set in opts.  Otherwise, an error is thrown.  If a
%                   field is not provided in user_opts, opts will contain
%                   the field with the default value.
%
%   PARAMETERS
% 
%   .kth_eigenvalue             [positive integer | {1}]
%       Selects the kth rightmost (continuous-time) or outermost
%       (discrete_time) eigenvalue that the algorithm should be initialized
%       at.  If .epsilon, .U and .V are all supplied by the user, this
%       selects the kth eigenvalue of the matrix M(epsilon*U*V'), where
%
%           M(Delta) = A + B Delta (I - D Delta)^{-1} C
% 
%       is the perturbed system matrix.  Otherwise, this selects the kth
%       eigenvalue of the system matrix A to use to compute (either
%       completely or partially) an initial perturbation.
%   
%   .epsilon                    [value in (0,1/norm(D)) | {[]}]
%       Initial guess for epsilon.  
%   
%   .U                          [U from outer product U*V', {[]}]
%       Initial guess for an initial unit-norm perturbation U*V'.
%
%   .V                          [V from outer product U*V', {[]}]
%       Initial guess for an initial unit-norm perturbation U*V'.
%
%   One can specify either epsilon or U*V', both, or neither.  Unprovided
%   values will be computed automatically.
%
%   Note that U*V' must be dimensionally compatible with the provided
%   system.  If opts.real_frobenius_norm is true, U and V must additionally
%   have at most 2 columns and each must be real valued.  Otherwise, U and
%   V must each be column vectors but they may be either real or complex
%   valued.  U*V' must also have unit norm, in the appropriate norm.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   initialPerturbationOptions.m introduced in ROSTAPACK Version 1.0
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
    
    persistent validator;
   
    if isempty(validator)
        validator = optionValidator(                                    ...
            'ROSTAPACK', getDefaults(), [], 'initial_perturbation'      );
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getValidatedOpts();
        return
    end
   
    validator.setUserOpts(user_opts);
  
    try 
        validator.setIntegerPositive('kth_eigenvalue'); 
        if validator.isSpecified('epsilon')
            validator.setRealPositive('epsilon');
            validator.setFiniteValued('epsilon');
        end
        has_U = validator.isSpecified('U');
        has_V = validator.isSpecified('V');

        if ~has_U && ~has_V
            opts = validator.getValidatedOpts();
            return
        end

        validator.assert(has_U,'Initial V provided without initial U!');
        validator.assert(has_V,'Initial U provided without initial V!');
        U_cols =  size(user_opts.U,2);
        V_cols =  size(user_opts.V,2);
        validator.assert(   U_cols < 3,                                 ...
                            'Initial U must have at most 2 columns!'    );
        validator.assert(   V_cols < 3,                                 ...
                            'Initial V must have at most 2 columns!'    );
        validator.assert(   U_cols == V_cols,                           ...
                'Initial U and V must have the same number of columns!' );

        validator.setNonZero('U');
        validator.setNonZero('V');
        validator.setFiniteValued('U');
        validator.setFiniteValued('V');
    catch err
        err.throwAsCaller();
    end
    
    opts = validator.getValidatedOpts();
end

function default_opts = getDefaults()
    default_opts = struct(  'kth_eigenvalue',   1,                      ...
                            'epsilon',          [],                     ...
                            'U',                [],                     ...
                            'V',                []                      );
end
