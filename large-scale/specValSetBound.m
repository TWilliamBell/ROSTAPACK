function [result,info] = specValSetBound(varargin)
%   specValSetBound:
%       Given a linear dynamical system, in state-space matrix form,
%       specValSetBound computes a guaranteed lower bound to the epsilon
%       spectral value set abscissa (or radius), where the spectral value
%       sets either allow:
% 
%       - both complex-valued and real-valued perturbations, bounded by the
%         spectral norm (or equivalently, the Frobenius norm)
%       - only real-valued perturbations, bounded by the Frobenius norm.
%
%       If the system is defined only by matrix A, or if B=C=I and D=0,
%       then specValSetBound computes a lower bound to the pseudospectral
%       abscissa|radius (the real structured pseudospectral abscissa|radius
%       if opts.real_frobenius_norm is set to true).
%       
%       In practice, specValSetBound often converges to the actual epsilon
%       spectral value set abscissa|radius, and even when it does not, it
%       still almost always converges to a lower bound that is at least
%       locally optimal and which is frequently a good approximation to the
%       epsilon spectral value set abscissa|radius.
%
%       The main computational cost of specValSetBound is a sequence of
%       sparse eigenvalue solves (using eigs/eigsPlus), which can allow it
%       to be efficient for large-scale systems.  The underlying algorithm
%       generally has a linear rate of convergence but extrapolation and
%       interpolation features are enabled by default to accelerate the
%       computation.
% 
%   USAGE:
%       [result,info] = specValSetBound(A,epsilon)
%       [result,info] = specValSetBound(A,epsilon,opts)
%       [result,info] = specValSetBound(A,B,C,D,epsilon)
%       [result,info] = specValSetBound(A,B,C,D,epsilon,opts)
%
%   INPUT:
%       System matrix A         [required]
%           Matrix A can be specified in any of the following formats:
%           - as an explicit matrix, dense or sparse: A
%           - as an explicit matrix, dense or sparse, enclosed in single 
%             element cell array: {A}
%           - as an outer product A = U*V': {U,V}
%           - as a single function handle for multiplying A and A':
%                   {applyA,rows,cols,is_real}
%             where:
%                   applyA(x,false) returns A*x
%                   applyA(x,true) returns A'*x 
%                   rows is the number of rows of matrix A
%                   cols is the number of columns of matrix A
%                   is_real:logical, true if A only contains real entries
%           - as two separate function handles for multiplying A and A':
%                   {applyA,applyAh,rows,cols,is_real}
%             where:
%                   applyA(x,false) returns A*x
%                   applyAh(x,true) returns A'*x 
%                   rows is the number of rows of matrix A
%                   cols is the number of columns of matrix A
%                   is_real: logical, true if A only contains real entries
% 
%           Note that the function handles should be able to multiply by A
%           and A' when x is either a vector or matrix.
%
%       System matrices B,C,D   [optional]
%           Matrices B,C,D can also be provided in any of above formats,
%           as well as via the following shortcuts:
%           B: [] indicates B is an appropriately-sized sparse identity.
%           C: [] indicates C is an appropriately-sized sparse identity.
%           D: [] indicates D is an appropriately-sized zero matrix.
%           If B and D are both [], m is automatically set to n.
%           If C and D are both [], p is automatically set to n.
%           If B,C,D are not provided, B=C=I and D=0, all n by n matrices.
% 
%       Note: when explicitly providing B=C=I and D=0, B and C should be
%       given as sparse identities (using speye or via function handles)
%       and D should be given as a sparse zero matrix.  Otherwise,
%       computation and memory usage may be very inefficient.
%
%       epsilon                 [required: real value in (0,1/norm(D))]
%           Specifies the perturbation level of the spectral value set.
%   
%       opts                    [optional: struct of parameters]
%           An optional struct of settable parameters or [].
%           To see available parameters and their descriptions, type:
%           >> help specValSetBoundOptions
%                 
%   OUTPUT:
%       result      
%           Struct containing the computed lower bound to the spectral
%           value set abscissa|radius and related information:
%
%       .halt_status             
%           0:  Maximum number of iterations reached
%           1:  Converged: relative difference condition attained
%               opts.rel_diff_tol satisfied
%           2:  Relative step size condition attained, halting early
%               opts.rel_step_size_tol satisfied
%           4:  Line search failed to produce a monotonic step
%           5:  Magnitude of line search's derivative sufficiently small
%           6:  New perturbation was identical to current perturbation
% 
%           Note that .halt_status == 3 is missing intentionally.
%           
%       .converged
%           Logical: true if .halt_status is 1, 4, 5, or 6 as these are all
%           indications of convergence to a locally rightmost/outermost
%           point of the spectral value set (or at least a stationary
%           point).
%
%       .svsar_lb
%           The computed lower bound to the epsilon spectral value set
%           abscissa or radius.
% 
%       .epsilon
%           The user's selected perturbation level.
% 
%       .U and .V
%           The perturbation vectors/matrices U and V such that .svsar_lb 
%           is given by the spectral abscissa|radius of 
%               M = A + B*(epsilon*U*V')*(I - D(epsilon*U*V'))^(-1)*C.
%        
%       .f
%           The value of the root function determining instability of the
%           system, which is either the spectral abscissa (continuous-time)
%           or the spectral radius minus 1 (discrete-time).  In the latter
%           case, .svsar_lb == .f + 1;
% 
%       .z 
%           The eigenvalue of M found by specValSetBound that provides the
%           value of .svsar_lb
% 
%       .x  The right eigenvector corresponding to z.
% 
%       .y and .absyhx      [optionally included] 
%           If opts.left_eigenvector is set to true, .y will be the
%           corresponding left eigenvector for eigenvalue .z of matrix M,
%           .absyhx = abs(y'*x).  Note that enabling this may cause an
%           additional eigenvalue solve to be incurred.  If sparse
%           eigensolves are being used, there is a possibility that this
%           additional required solve can fail.  In this case, the left
%           eigenvector and abs(y'*x) will each be returned as [].
%
%       info
%           Struct containing metadata about the computation process to
%           produce the computed lower bound.  The following two fields are
%           always included:
%       
%       .iters
%           Number of iterations incurred before halting.
%   
%       .cost 
%           Struct containing the following stats related to the total cost 
%           of the computation:
%
%           .eig_solves_left    total number of left eigenvector solves
%           .eig_solves_right   total number of right eigenvector solves
%                      
%           The following are only provided when sparse computations are
%           used and eigsPlus is installed:
%
%           .eigs_iters_left    total number of eigs iterations for 
%                               computing all the left eigenvectors
%           .eigs_iters_right   total number of eigs iterations for 
%                               computing all the right eigenvectors
%
%           The remaining three subfields of .cost are only provided when
%           opts.real_frobenius_norm and opts.complex_ode_iteration are
%           both set to false and D is a nonzero matrix.  In this case,
%           each expansion step requires solving the pair of linear systems
%           described in phiSystemsSolverOptions.  Note that solving
%           this pair of systems is generally a negligible cost.
% 
%           .phi_pair_solves    total number of times these particular pair
%                               of linear systems was solved
%
%           If PCG was selected as the method of choice to solve the pair
%           of linear systems, the following two fields will also present:
%
%           .phi_m_iters       
%           .phi_p_iters
%
%           which provide the accumulated number of PCG iterations incurred
%           to solve each of the two systems in the pairs, respectively.
%
%       .multiplies     [only present if opts.count_multiplies is true]
%           A struct containing the total number of matrix-vector products
%           per matrix, i.e. for A, A', B, B', etc.
%
%       .initialization [only present if opts.record_level > 0]
%           The info output argument from initializePerturbation
%
%       .expansion      [only present if opts.record_level > 0]
%           The info output argument from uvExpand
%   
%   CONSOLE:
%       COMMON COLUMN HEADERS
%       Iter
%           The current iteration number of the relevant routine/phase
%       Epsilon
%           Value   Current value of epsilon in Delta = epsilon*U*V'
%           Diff    Absolute difference between current and previous values
%       Root Fn   
%           Value   Current value of r(z), where r(z) is either 
%                       r(z) = real(z)      (continuous time)
%                       r(z) = abs(z) - 1   (discrete time)
%                   and z is (typically) a globally rightmost/outermost
%                   eigenvalue of
%                       M(Delta) = A + B*Delta*(I - D*Delta)^{-1}*C
%                   for Delta = epsilon*U*V' with U*V' of unit norm.  Note
%                   that the minimal norm of Delta such that r(z) is zero 
%                   is the complex|real stability radius.
%           Diff    Absolute difference between current and previous values
%       
%       Note that:
%       1)  The value of epsilon (and the difference with its previous 
%           value) will only be printed when Iter is 0 or if epsilon has
%           been changed.  Otherwise, these two columns will remain blank.
%       2)  The values of epsilon and/or r(z) will be printed in ()'s for
%           rejected steps.  Rejected steps are only printed for
%           sufficiently high values of opts.print_level.
%
%       INITIALIZATION SPECIFIC COLUMN HEADERS
%       Matrix   
%           Type 
%               A       z is an eigenvalue of M(0) = A     
%               M0      z is an eigenvalue of M(Delta), Delta nonzero 
%           Eig. #
%               k       z is the kth eigenvalue of either A or M0, where
%                       the eigenvalues are in decreasing order of real
%                       part (continuous time) or modulus (discrete time)
%
%       Note that the r(z) for an eigenvalue of M0 may be less than its
%       value for an eigenvalue of A.
%
%       EXPANSION SPECIFIC COLUMN HEADERS
%       Step Type
%            -      initial value
%           UV      regular UV' update step
%           LS      rejected line search evaluation 
%           E       extrapolated step
%           I       interpolated step
%       Line Search
%           #       UV: total number of line search evaluations 
%                   LS: total number of line search evaluations so far
%                   (not applicable for E, I, or the initial value)
%           t       UV: final step length in (0,1] (typically accepted)
%                   LS: current trial step length
%                   I:  optimal length of step determined by interpolation
%                   (not applicable for E or the initial value)
%
%   See also specValSetBoundOptions and getStabRadBound.
%
%   
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   specValSetBound.m introduced in ROSTAPACK Version 1.0
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

    try 
        [A,B,C,D,epsilon,opts] = processInputs(varargin{:});
    catch err
        err.throwAsCaller(); 
    end
    info_requested = nargout > 1;
    if ~info_requested && nargin > 5
        % for efficiency, make sure recording is completely disabled when
        % the caller does not ask for the corresponding output argument
        opts.record_level = 0;
    end

    try 
        % validate all the options before doing any computations
        opts = specValSetBoundOptions(opts);
    catch err
        err.throwAsCaller();
    end
    
    [print_level,printer,printBox] = configurePrinting(opts);
    if print_level
        fprintf('\n\n\n');
        printBox(false,'COPYRIGHT NOTICE',copyrightNotice(mfilename()));    
        if opts.print_info_msg && opts.discrete_time
            printBox(   opts.print_use_orange,                          ...
                        'PRINT INFO MESSAGE',printInfoMsg()             );
        end
    end
    
    % Add .epsilon field, but set to empty.  initialPerturbationAssert
    % requires that the .epsilon field exist but it will not assert its
    % value if it is not set.  We will set and check it separately, in
    % order to control the output of the error messages.
    epsUV0          = opts.initial_perturbation;
    epsUV0.epsilon  = [];
    real_fro_norm   = opts.real_frobenius_norm;
    
    % Tell stateSpaceABCD that it can forgo verifying its opts, since they
    % have already been verified by getStabRadBoundOptions.
    opts.verified   = true;
    
    try 
        sys         = stateSpaceABCD(A,B,C,D,opts); 
        limit       = epsilonUpperBound(sys);
        assert( epsilon < limit,                                        ...
                'For this system, epsilon must be less than %g!',limit  );
        epsUV0      = initialPerturbationAssert(sys,real_fro_norm,epsUV0);
    catch err
        err.throwAsCaller();
    end
    
    epsUV0.epsilon  = epsilon;
    sys_c           = [];
    
    if opts.record_level > 0 
        info_init   = cell(1);
    else
        info_init   = {};
        phases      = {};
    end
    
    if real_fro_norm 
        assert(sys.isReal(),rsrErrorMsg());
        [~,m,p] = sys.getDimensions();
        assert(max(m,p) > 1,rsrErrorDimenionsMsg());
        sys_c = systemRealFrobeniusNorm(sys);
    elseif opts.complex_ode_iteration
        sys_c = systemComplexTwoNormODE(sys);
    else
        sys_c = systemComplexTwoNorm(sys,opts.phi_solver_opts);
    end
    
    if print_level
        fprintf('\n\n\n');
    end
         
    % compute initial perturbation or (partially) use user-supplied delta                                                       
    [~,info_init{:}] = initializePerturbation(sys_c,initPertOpts());  
                 
    expand_opts                 = opts.expansion_opts;
    expand_opts.record_UV       = opts.record_UV;
    expand_opts.record_level    = max(opts.record_level - 1,0);
    expand_opts.print_level     = print_level;
    if print_level 
        expand_opts.printer     = printer.expand;
    end
   
    [~,~,UV_info] = uvExpand(sys_c,expand_opts);
    
    prepareExit(UV_info.halt_status);
    
    function prepareExit(status)
        [f,z]       = sys_c.getf();
        epsilon     = sys_c.getEpsilon();
        [U,V]       = sys_c.getUV();
        if opts.left_eigenvector
            [x,y,absyhx]    = sys_c.getEigenvectors();
            evecs           = {'x',x,'y',y,'absyhx',absyhx};
        else
            evecs   = {};
        end
        svsar       = f;
        if opts.discrete_time 
            svsar   = svsar + 1;
        end
        
        result  = struct(   'halt_status',  status,             ...
                            'converged',    UV_info.converged,  ...  
                            'svsar_lb',     svsar,              ...
                            'epsilon',      epsilon,            ...
                            'U',            U,                  ...
                            'V',            V,                  ...
                            'f',            f,                  ...
                            'z',            z,                  ...
                            evecs{:}                            );
        if info_requested
            if opts.record_level > 0 
                phases      = { 'initialization',   info_init,      ...
                                'expansion',        UV_info         };
            end
            if opts.count_multiplies
                multiplies  = {'multiplies', sys.getCounts()};
            else
                multiplies  = {};
            end
            info = struct(  'iters',            UV_info.stats.iters,...
                            'cost',             sys_c.getTotals(),  ...
                            multiplies{:},                          ...                         
                            phases{:}                               ); 
        end
    end
    
    function init_opts = initPertOpts()
        init_opts     = opts;
        
        init_opts.check_stability = false;
        
        % initialPerturbationAssert may renormalize u and v
        init_opts.initial_perturbation        = epsUV0;
        
        init_opts.record_UV       = opts.record_UV;
        init_opts.record_level    = min(max(opts.record_level - 1,0),1);
        
        if ~isempty(printer)
            init_opts.printer     = printer.init;
        end
    end
end

function [A,B,C,D,epsilon,opts] = processInputs(varargin)
    switch nargin 
        case {2,3}
            A           = varargin{1};
            [B,C,D]     = deal([]);
            vars        = varargin(2:end); 
        case {5,6}
            [A,B,C,D]   = deal(varargin{1:4});
            vars        = varargin(5:end);
        otherwise
            error('specValSetBound:invalidUsage',invalidUsageMsg()); 
    end
    epsilon         = vars{1};
    assert( isARealNumber(epsilon) && epsilon > 0 && ~isinf(epsilon),   ...
            'Epsilon must be a positive real number!'                   );
    if length(vars) > 1
        opts        = vars{2};
        assert(isstruct(opts),'opts must be a struct!');
    else
        opts        = [];
    end
end

function [print_level,printer,printBoxFn] = configurePrinting(opts)
        
    print_level = opts.print_level;
    printer     = [];
    printBoxFn  = [];
    
    % nothing to do if printing is disabled
    if print_level < 1
        return
    end
    
    % Check if printer is already configured and if so, return it
    if isfield(opts,'printer')  
        printer = opts.printer;
        return
    end
    
    % Otherwise we need to build a new printer
    print_ascii     = opts.print_ascii;
                       
    col_specs       = getColumnSpecs(opts);
    cols            = printerColumnFormatters(col_specs);
    printer.init    = printerInit(print_ascii,cols);
    printer.expand  = printerExpand(print_ascii,cols);
    
    printBoxFn      = @(use_orange,title,lines) printMessageBox(        ...
                        print_ascii, use_orange, 2, title, [], lines    );  
end

function c = getColumnSpecs(opts)
    expansion_opts  = opts.expansion_opts;
    maxit           = expansion_opts.maxit;
    ls_maxit        = expansion_opts.line_search_opts.maxit;
    kth_eigenvalue  = opts.initial_perturbation.kth_eigenvalue;
    c = struct( 'maxit',            maxit,          ...
                'step_type',        true,           ...
                'ls_maxit',         ls_maxit,       ...
                'kth_eigenvalue',   kth_eigenvalue  );
end

function m = invalidUsageMsg()
m = [                                                                   ...
'Invalid usage: specValSetBound may be called in the following ways:\n' ...
'  [result,info] = specValSetBound(A,epsilon)\n'                        ...
'  [result,info] = specValSetBound(A,epsilon,opts)\n'                   ...
'  [result,info] = specValSetBound(A,B,C,D,epsilon)\n'                  ...
'  [result,info] = specValSetBound(A,B,C,D,epsilon,opts)\n'             ];
end

function m = rsrErrorMsg()
m = [                                                                   ...
'System cannot contain complex matrices if opts.real_frobenius_norm is '...
'set to true!'                                                          ];
end

function m = rsrErrorDimenionsMsg()
m = [                                                                   ...
'Error: max(m,p) must be at least 2 if opts.real_frobenius_norm is '    ...
'set to true!'                                                          ];   
end

function s = printInfoMsg()
s = {                                                                           ...
'NOTE: specValSetBound() prints out the history of the root function for '      ...
'determining whether or not a system is stable (under perturbation), with '     ...
'nonnegative values indicating instability.  In the continuous-time case, '     ...
'the root function is simply the spectral value set abscissa but, for '         ...
'discrete-time systems, the root function is the spectral value set radius '    ...
'MINUS ONE.  Thus, for this discrete-time system, all printed values with '     ...
'be one less than the approximation to the spectral value set radius '          ...
'computed at each iteration.'                                                   ...
''                                                                              ...
'To disable this notice, set opts.print_info_msg = false.'                      ...
};
end
