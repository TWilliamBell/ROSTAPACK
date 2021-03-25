function opts = eigenSolverOptions(user_opts)
%   eigenSolverOptions:
%       Validate user options struct of parameters for eigenvalue solves.
%       If user_opts is [] or not provided, returned opts will be the
%       default parameters.
%       
%   USAGE:
%       opts = eigenSolverOptions();
%       opts = eigenSolverOptions(user_opts);
%
%   INPUT:
%       user_opts   Struct of settable parameters.  No fields are required,
%                   irrelevant fields are ignored, and user_opts may be 
%                   given as [].
%   
%   OUTPUT:
%       opts        Struct of all user-tunable parameters.  If a field is
%                   provided in user_opts, then the user's value is checked
%                   to whether or not it is a valid value, and if so, it is
%                   set in opts.  Otherwise, an error is thrown.  If a
%                   field is not provided in user_opts, opts will contain
%                   the field with the default value.
%
%   PARAMETERS
% 
%   .maxit                      [positive integer | {300}]
%       The maximum number of eigs/eigsPlus iterations to allow per solve.
%   
%   .tol                        [positive finite real value | {1e-12}]
%       Default convergence tolerance value for eigs/eigsPlus.
%   
%   .request_k_eigenvalues      [positive integer | {8}]
%       Default number of eigenvalues to request for eigs/eigsPlus.
% 
%   .krylov_dim                 [positive integer | {[]}]
%       Size of the Krylov subspace eigs/eigsPlus should use.  When this is
%       set to [], eigs/eigsPlus will automatically determine the Krylov
%       subspace dimension based on the problem parameters.  For more
%       details, see help eigs.
%
%   .v0_recycling_type          [integer in {0,1,2} | {2}]
%       This determines how the initial vector for eigs/eigsPlus is set for
%       each computation:
%       - 0: initial vector is randomly generated
%       - 1: initial vector is eigenvector corresponding to the desired / 
%            "selected" eigenvalue from the previous eigenvalue computation
%       - 2: initial vector is the average of all eigenvectors computed
%            in the previosu eigenvalue computation
% 
%   .use_default_eigs           [logical | {false}]
%       Set to this true if you wish to fallback to using MATLAB's default
%       eigs in lieu of eigsPlus.  The only likely reason to do fallback 
%       to the default eigs is when the latest version of MATLAB breaks
%       comptability with eigsPlus and eigsPlus has not yet been updated to
%       fix this.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   eigenSolverOptions.m introduced in ROSTAPACK Version 1.0
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
            'ROSTAPACK', getDefaults(), [], 'eig_solver_opts'           );
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
   
    validator.setUserOpts(user_opts);
  
    try 
        validator.setIntegerPositive('maxit');
        validator.setRealInIntervalOO('tol',0,inf);
        validator.setIntegerPositive('request_k_eigenvalues');
        
        if validator.isSpecified('krylov_dim')
            validator.setIntegerPositive('krylov_dim');
        end
        
        validator.setIntegerInRange('v0_recycling_type',0,2);
        validator.setLogical('use_default_eigs');
    catch err
        err.throwAsCaller();
    end
    
    opts = validator.getValidatedOpts();
end

function default_opts = getDefaults()
    default_opts = struct(  'maxit',                    300,        ...
                            'tol',                      1e-12,      ...
                            'request_k_eigenvalues',    8,          ...
                            'krylov_dim',               [],         ...
                            'v0_recycling_type',        2,          ...
                            'use_default_eigs',         false       );
end
