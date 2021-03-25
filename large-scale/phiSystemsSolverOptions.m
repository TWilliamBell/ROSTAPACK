function opts = phiSystemsSolverOptions(user_opts)
%   phiSystemsSolverOptions:
%       This routine creates and validates a struct of user options related
%       to solving the following pair of linear systems:
%       
%           Phi_m = (I - epsilon^2 D' * D)
%           Phi_p = (I - epsilon^2 D * D'),
%
%       which arise when approximating the complex-valued spectral value
%       set abscissa and/or radius.  If user_opts is [] or not provided,
%       returned opts will be the default parameters.
%
%       There is generally no reason to change the default options, unless
%       solving these systems becomes difficult.
%
%       In general, these linear systems can be solved via direct and
%       efficient methods and the most efficient method will be selected
%       automatically.  If D is large (in both dimensions), then pcg can be
%       used to iteratively solve the systems without expensive dense
%       matrix operations.
%       
%       Note that suboptions are only relevant when:
%           opts.real_frobenius_norm 
%           opts.complex_ode_iteration 
%       are both false.
%       
%   USAGE:
%       opts = phiSystemsSolverOptions();
%       opts = phiSystemsSolverOptions(user_opts);
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
%                   set in opts. Otherwise, an error is thrown.  If a field
%                   is not provided in user_opts, opts will contain the
%                   field with the default value.
%
%   PARAMETERS
% 
%   .force_method               [value in {0,1,2,3,4,5} | {0}]
%       Force a particular method to be used:
%           0 - D is zero, nothing to solve, so can't/don't force a method.
%           1 - Use back slash
%           2 - Use two Cholesky factorizations of Phi_m and Phi_p
%           3 - Use one Cholesky factorization: smaller of Phi_m and Phi_p
%           4 - Exploit low-rank structure D when given as D = Dl*Dr'
%           5 - Solve the systems iteratively using pcg.
%       By default, the most efficient method amongst {0,3,4,5} will be
%       chosen automatically, based on the user-supplied representation of
%       D but the user may optionally override this choice if they desire
%       (provided the supplied D has a representation that is
%       compatible with one's choice of method).
% 
%    .maxit                      [positive integer or empty | {[]}]
%       The maximum number of pcg iterations to allow per solve.  When this
%       is set to [], pcg automatically determines its own maxit based off
%       of the dimension of the matrix it is given.  Note that this
%       parameter only has an effect when pcg is used.
%
%   .tol                        [positive finite real value | {1e-12}]
%       Default convergence tolerance value for pcg.
%   
%   .initial_vector_type        [integer in {0,1,2,3} | {3}]
%       If the linear solves are being done via pcg, this determines what 
%       initial vectors are used for pcg for the pair of linear systems:
%           0:  zero vectors
%           1:  randomly generated
%           2:  the solutions to the previous pair of solves
%           3:  optimally rotated versions of the solutions to the previous
%               pair of solves.  This is generally the most efficient
%               choice.
%
%
%   For more details, see [GGO13, Section 3], [Mit14, Section 4.1] and 
%   [MO16, Section 6].
%   
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   phiSystemsSolverOptions.m introduced in ROSTAPACK Version 1.0
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
            'ROSTAPACK', getDefaults(), [], 'phi_solver_opts'           );
    end
    
    if nargin < 1 || isempty(user_opts)
        opts = validator.getDefaultOpts();
        return
    end
    
    validator.setUserOpts(user_opts);
    
    try 
        validator.setIntegerInRange('force_method',0,5);
        if validator.isSpecified('maxit') 
            validator.setIntegerPositive('maxit');          
        end
        validator.setRealInIntervalOO('tol',0,inf);
        validator.setIntegerInRange('initial_vector_type',0,3); 
    catch err
        err.throwAsCaller();
    end
    
    opts = validator.getValidatedOpts();
end

function default_opts = getDefaults()
    default_opts = struct(  'force_method',         0,      ...
                            'maxit',                [],     ...
                            'tol',                  1e-12,  ...
                            'initial_vector_type',  3       );
end
