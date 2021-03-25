function [result,info] = getStabRadBound(varargin)       
%   getStabRadBound:
%       Given a linear dynamical system, in state-space matrix form,
%       getStabRadBound computes a guaranteed upper bound to the complex
%       stability radius (equivalently a lower bound to the H-infinity
%       norm).  If the system is defined only by matrix A (B=C=I and D=0),
%       then complex stability radius is known as the distance to
%       instability.  Alternatively, getStabRadBound can instead compute a
%       guaranteed upper bound to the real stability radius, using the
%       Frobenius norm (see opts.real_frobenius_norm).
% 
%       Whether the continuous- or discrete-time stability radius is
%       approximated is selected via opts.discrete_time.
%       
%       In practice, the bounds computed by getStabRadBound are often tight
%       (i.e. the bounds are the complex|real stability radius values).
%       If not, the bounds still are almost always locally optimal and
%       frequently good approximatiomns to the actual stability radius
%       values.
%
%       The main computational cost of getStabRadBound is a sequence of
%       sparse eigenvalue solves (using eigs/eigsPlus), which can allow it
%       to be efficient for large-scale systems.  This routine implements 
%       the Hybrid Expansion-Contraction (HEC) algorithm [MO16], which is
%       generally quadratically convergent.
% 
%   USAGE:
%       [result,info] = getStabRadBound(A)
%       [result,info] = getStabRadBound(A,opts)
%       [result,info] = getStabRadBound(A,B,C,D)
%       [result,info] = getStabRadBound(A,B,C,D,opts)
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
%       opts                    [optional: struct of parameters]
%           An optional struct of settable parameters or [].
%           To see available parameters and their descriptions, type:
%           >> help getStabRadBoundOptions
%                 
%   OUTPUT:
%       result      
%           Struct containing the computed upper bound to the complex|real
%            stability radius and related information:
%
%       .halt_status  
%          -1:  A was not stable
%           0:  Failed to find an initial upper bound. 
%           1:  HEC converged to tolerances.
%           2:  HEC stagnated as it can no longer make any contraction 
% `             or expansion progress but the tolerances have not been
%               satisfied.  This is generally an indication that the
%               tolerances are too tight for the particular problem.
%           3:  Maximum number of iterations reached.
%           
%       .epsilon
%           The computed upper bound to the complex|real stability radius.
% 
%       .gain
%           The reciprocal of .epsilon.  When approximating the complex
%           stability radius, .gain is a lower bound to the H-infinity
%           norm.  
% 
%       .U and .V
%           The final perturbation vectors/matrices U and V such that the
%           for perturbation Delta = .epsilon*.U*.V', the
%           rightmost/outermost eigenvalue z of
%               M = A + B*Delta*(I - D*Delta)^(-1)*C
%           lies on the stability boundary (imaginary axis or unit circle).
%     
%       .f
%           The value of the root function determining instability of the
%           system, which is either real(z) (continuous-time) or abs(z)-1
%           (discrete-time), evaluated at the eigenvalue z described
%           immedatiately above.
%
%       .z 
%           The eigenvalue z described above.
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
%           produce the computed upper bound.  The following four fields
%           are always included:
%       
%       .iters_upperbound
%           Number of iterations to find an initial upper bound.
%
%       .iters_frequency
%           Number of iterations trying to improve the initial frequency
%           guess once an initial upper bound has been obtained.  
%           
%       .iters_hec
%           Number of HEC iterations incurred before halting.
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
%       .upperbound     [only present if opts.record_level > 0]
%           The info output argument from findUpperBound
% 
%       .hec            [only present if opts.record_level > 0]
%           The info output argument from hybridExpansionContraction
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
%       UPPER BOUND SPECIFIC COLUMN HEADERS
%       Phase
%           Type 
%                -      initial value
%               EPS     increasing epsilon phase (done via epsilonExpand)
%               UV      updating U*V' phase (done via uvExpand)
%           Iters  
%               Number of steps taken by epsilonExpand or uvExpand
%           TC
%               Termination code returned by epsilonExpand or uvExpand  
%
%       HYBRID EXPANSION-CONTRACTION SPECIFIC COLUMN HEADERS
%       Phase
%           Type 
%                -      initial value
%               CON     increasing epsilon phase (done via epsilonContract)
%               EXP     updating U*V' phase (done via uvExpand)
%           Iters  
%               Number of steps taken by epsilonContract or uvExpand
%           TC
%               Termination code returned by epsilonContract or uvExpand  
% 
%       Note that when opts.print_level > 1 (the default is 1), the
%       printing is also enabled for the subroutines epsilonExpand, 
%       epsilonContract, and uvExpand, thus showing progress of the
%       subroutines in addition to the overall algorithm.
%
%   See also getStabRadBoundOptions and specValSetBound. 
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   getStabRadBound.m introduced in ROSTAPACK Version 1.0 
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
        [A,B,C,D,opts] = processInputs(varargin{:});
    catch err
        err.throwAsCaller(); 
    end
    info_requested = nargout > 1;
    if ~info_requested && nargin > 4
        % for efficiency, make sure recording is completely disabled when
        % the caller does not ask for the corresponding output argument
        opts.record_level = 0;
    end

    try 
        % validate all the options before doing any computations
        opts = getStabRadBoundOptions(opts);
    catch err
        err.throwAsCaller();
    end
    
    [print_level,opts.printers] = printersConfig(opts);
    
    if print_level
        fprintf('\n\n\n');
        printMessageBox(opts.print_ascii, false, 2, 'COPYRIGHT NOTICE', ...
                        [], copyrightNotice(mfilename())                );        
    end
     
    epsUV0          = opts.initial_perturbation;
    real_fro_norm   = opts.real_frobenius_norm;
    
    % Tell stateSpaceABCD that it can forgo verifying its opts, since they
    % have already been verified by getStabRadBoundOptions.
    opts.verified   = true;
    
    try 
        sys         = stateSpaceABCD(A,B,C,D,opts); 
        epsUV0      = initialPerturbationAssert(sys,real_fro_norm,epsUV0);
    catch err
        err.throwAsCaller();
    end
    
    sys_c               = [];
    iters               = 0;
    iters_ub.upperbound = 0;
    iters_ub.frequency  = 0;
    
    if opts.record_level > 0 
        info_init   = cell(1);
        info_ub     = cell(1);
        info_hec    = cell(1);
    else
        info_init   = {};
        info_ub     = {};
        info_hec    = {};
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
    [A_stable,info_init{:}] = initializePerturbation(sys_c,initPertOpts());  
       
    if ~A_stable
        prepareExit(-1);
        return
    end
    
    if print_level
        fprintf('\n\n\n');
    end
   
    % find upper bound for epsilon
    [found_upperbound,iters_ub,info_ub{:}] = findUpperBound(sys_c,opts);
    
    if ~found_upperbound
        prepareExit(0);
        return
    end
    
    if print_level
        fprintf('\n\n\n');
    end
    
    % convergent phase
    [halt_code,iters,info_hec{:}] = hybridExpansionContraction(sys_c,opts);
    
    prepareExit(halt_code);
    
    function prepareExit(status)
        epsilon     = sys_c.getEpsilon();
        [U,V]       = sys_c.getUV();
        [f,z]       = sys_c.getf();
        if status > 0 && opts.left_eigenvector
            [x,y,absyhx]    = sys_c.getEigenvectors();
            evecs           = {'x',x,'y',y,'absyhx',absyhx};
        else
            evecs   = {};
        end
        result  = struct(   'halt_status',  status,     ...
                            'epsilon',      epsilon,    ...
                            'gain',         1/epsilon,  ...
                            'U',            U,          ...
                            'V',            V,          ...
                            'f',            f,          ...
                            'z',            z,          ...
                            evecs{:}                    );
        if info_requested
            if opts.record_level > 0 
                phases      = { 'initialization',   info_init,      ...
                                'upperbound',       info_ub,        ...
                                'hec',              info_hec         };
            end
            if opts.count_multiplies
                multiplies  = {'multiplies', sys.getCounts()};
            else
                multiplies  = {};
            end
            info = struct(  'iters_upperbound', iters_ub.upperbound,...
                            'iters_frequency',  iters_ub.frequency, ...
                            'iters_hec',        iters,              ...
                            'cost',             sys_c.getTotals(),  ...
                            multiplies{:},                          ...
                            phases{:}                               ); 
        end
    end
    
    function init_opts = initPertOpts()
        init_opts   = opts;
        
        eps_lf      = opts.upperbound_opts.epsilon_limit_fraction;
        eps_sm      = opts.upperbound_opts.epsilon_step_multiplier;
        init_opts.epsilon_limit_fraction    = eps_lf;
        init_opts.epsilon_step_multiplier   = eps_sm;
        
        % initialPerturbationAssert may renormalize u and v
        init_opts.initial_perturbation      = epsUV0;
        
        init_opts.record_UV         = opts.record_UV;
        init_opts.record_level      = min(max(opts.record_level - 1,0),1);
        
        if init_opts.print_level
            init_opts.printer       = opts.printers.init;
        end
    end
end

function [A,B,C,D,opts] = processInputs(varargin)
    switch nargin 
        case {1,2}
            A           = varargin{1};
            [B,C,D]     = deal([]);
            vars        = varargin(2:end); 
        case {4,5}
            [A,B,C,D]   = deal(varargin{1:4});
            vars        = varargin(5:end);
        otherwise
            error('getStabRadBound:invalidUsage',invalidUsageMsg()); 
    end
    if isempty(vars)
        opts        = [];
    else
        opts        = vars{1};
        assert(isstruct(opts),'opts must be a struct!');
    end
end

function m = invalidUsageMsg()
m = [                                                                   ...
'Invalid usage: getStabRadBound may be called in the following ways:\n' ...
'  [result,info] = getStabRadBound(A)\n'                                ...
'  [result,info] = getStabRadBound(A,opts)\n'                           ...
'  [result,info] = getStabRadBound(A,B,C,D)\n'                          ...
'  [result,info] = getStabRadBound(A,B,C,D,opts)\n'                     ];
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
