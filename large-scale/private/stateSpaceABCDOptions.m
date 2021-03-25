function opts = stateSpaceABCDOptions(user_opts)
%   stateSpaceABCDOptions:
%       Validate user options struct of parameters for the stateSpaceABCD
%       object.  If user_opts is [] or not provided, returned opts will be
%       the default parameters.
%       
%   USAGE:
%       opts = stateSpaceABCDOptions();
%       opts = stateSpaceABCDOptions(user_opts);
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
%   .discrete_time              [logical | {false}]
%       By default, systems are considered to be continuous time and thus
%       the root function for stability is the spectral abscissa.  If the
%       system is instead discrete time, then set .discrete_time to true; 
%       the root function for stability will then be the spectral radius
%       minus one.  
%
%   .sparse_mode                [integer in {0,1,2} | 2]
%       Controls whether dense or sparse eigenvalue solvers will be used.
%           0 - force dense eigenvalue solves
%           1 - force sparse eigenvalue solves
%           2 - automatically select dense or sparse eigenvalue solves
%               based on how the system matrix A is provided.  If A is
%               given as a sparse matrix, an outer product, or as a
%               function handle, sparse eigenvalues solves will be used.
%
%   .eig_solver_opts            [struct | {eigSolverOptions()}]
%       This is a struct of the parameters for sparse eigenvalue solves.
%       If dense eigensolves are being used, this struct of parameters is
%       ignored.  For more details, see eigenSolverOptions.
%
%   .count_multiplies           [logical | {false}]
%       When this is set to true, stateSpaceABCD will keep track of the
%       total number of matrix multiplications for A,A',B,B',C,C',D,D'.
%       Note enabling this adds some overhead.
%
%   See also stateSpaceABCD.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   stateSpaceABCDOptions.m introduced in ROSTAPACK Version 1.0
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
        [defaults,sub_validators] = getDefaults();
        validator = optionValidator(                                    ...
            'ROSTAPACK',defaults,sub_validators,'stateSpaceABCD'        );
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
   
    validator.setUserOpts(user_opts);
  
    validator.setLogical('discrete_time');  
    validator.setIntegerInRange('sparse_mode',0,2);
    validator.setStructWithValidation('eig_solver_opts');
    validator.setLogical('count_multiplies');

    opts = validator.getValidatedOpts();
end

function [default_opts,sub_validators] = getDefaults()
    default_opts = struct(  'discrete_time',    false,                  ...
                            'sparse_mode',      2,                      ...
                            'eig_solver_opts',  [],                     ...
                            'count_multiplies', false                   );
    sub_validators = struct('eig_solver_opts',  @eigenSolverOptions);     
end
