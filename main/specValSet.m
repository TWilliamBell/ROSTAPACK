function [eta,loc,info] = specValSet(A,varargin)
%   specValSet:
%       Given a linear dynamical system, in state-space matrix form,
%
%           E xdot = Ax + Bu
%                  = Cx + Du,
% 
%       specValSet computes the epsilon spectral value set abscissa or
%       radius.  When B=C=E=I and D=0, with n=m=p, then specValSet computes 
%       the epsilon pseudospectral abscissa or radius.  E should be
%       invertible.
%       
%       Computationally, specValSet requires O(n^3) work and O(n^2) memory
%       per iteration.  Every iteration, a dense matrix pencil of order 2n
%       is formed and all of its eigenvalues are computed.  The underlying
%       algorithm has a quadratic rate of convergence.
% 
%   USAGE:
%       [eta,loc,info] = specValSet(A,epsilon)
%       [eta,loc,info] = specValSet(A,epsilon,opts)
%       [eta,loc,info] = specValSet(A,B,C,D,E,epsilon)
%       [eta,loc,info] = specValSet(A,B,C,D,E,epsilon,opts)
%
%   INPUT:
%       System matrix A         [required]
%           Matrix A must be provided as an explicit dense or sparse matrix
%       
%       System matrices B,C,D,E [optional]
%           Matrices B,C,D,E must be provided as explicit dense or sparse
%           matrices.  Setting any to [] acts as a shortcut for their
%           corresponding default values:
%           B: [] indicates B is a n by m (sparse) identity.
%           C: [] indicates C is a p by n (sparse) identity.
%           D: [] indicates D is the zero matrix.
%           E: [] indicates E is a n by n (sparse) identity.
%           If B and D are both [], m is automatically set to n.
%           If C and D are both [], p is automatically set to n.
% 
%           When opts.fast_search is enabled, if the A,B,C,D,E matrices are
%           indeed sparse and are provided in sparse format, there can be
%           additional efficiency benefits in the horizontal/radial
%           searches.  In the general case, the fast searches require
%           applying (zE-A)^{-1} for complex scalars z, which by default is
%           done via Hessenberg factorizations.  If A and E are large but
%           then this can optionally be done with LU factorizations.  When
%           B=C=I, D=0, and n=m=p, an alterate more efficient formulation
%           is used by default, where only zE-A is required (and not its
%           inverse), thus avoiding any factorizations and backsolves.
%           However, as this alterate form involves computing the smallest
%           singular value, numerical accuracy can be adversely affected
%           for some very poorly scaled matrices.  In this case, one can
%           set opts.force_two_norm to true to force the the computation
%           involving (zE-A)^{-1}.
% 
%       epsilon                 [required: real value in [0,1/norm(D))]
%           Specifies the perturbation level of the spectral value set.
%           This can be set to zero, to compute the spectral abscissa or
%           radius of (A,E), with or without respect to infinite and/or
%           uncontrollable/unobservable eigenvalues, depending on how opts
%           is set.
%   
%       opts                    [optional: struct of parameters]
%           An optional struct of settable parameters or [].
%           To see available parameters and their descriptions, type:
%           >> help specValSetOptions
%                 
%   OUTPUT:
%       eta      
%           The computed value of the epsilon spectral value set abscissa
%           or radius.
%
%       loc
%           The location of where the computed globally rightmost/outermost
%           point is, i.e. its imaginary part or angle.  The computed
%           rightmost/outermost point is:
%               abscissa case:  eta + 1i*loc
%               radius case:    eta*exp(1i*loc)
%
%       info
%           Struct containing metadata about the computation.
%       
%       .iters
%           Number of iterations incurred before halting.
%   
%       .eig_count
%           Number of times all eigenvalues of the matrix pencils of order
%           2n were computed.  When opts.fast_search is enabled, this is
%           equals to opts.iters, i.e. the number of vertical/circular
%           searches performed.  When opts.fast_search is disabled, this
%           this is the total number of all the vertical/circular and 
%           horizontal/radial searches.
%
%       .svd_count 
%           Number of times an SVD of a p times m matrix was computed.  In
%           the spectral value set case, this corresponds to computing the
%           norm of the transfer function.  If B=C=I, D=0, and n=m=p (and
%           opts.force_two_norm is false), then this corresponds to
%           computing the reciprocal of the minimum singular value of zE-A
%           for some complex scalar z.
%
%       .brk_steps
%           Only nonzero if opts.fast_search is enabled.  This is the total
%           number of steps incurred for finding upper bounds for all the
%           root problems encountered, in order to establish an initial
%           bracket for each root-finding problem. 
%
%       .root_steps
%           Only nonzero if opts.fast_search is enabled.  This is the total
%           number of steps incurred for finding roots over all problems,
%           but not including the steps necessary to first establish the
%           initial brackets.  
%
%       .full_steps
%           Only nonzero if opts.fast_search is enabled.  This is the total
%           subset of steps of info.root_steps that were not bisection
%           steps.
%       
%       .hr_problems
%           This is a vector containing the total number of possible
%           horizontal/radial searches per iteration.
%
%       .hr_solved
%           This is a vector containing the total number of
%           horizontal/radial searches per iteration that were actually
%           solved.  If opts.fast_search is disabled, then this will be
%           identical to info.hr_problems.
%       
%       .n_perts
%           Only nonzero if opts.fast_search is enabled.  This is the total
%           number of small perturbations possibly done after each
%           horizontal/radial search phase, in order to ensure the computed
%           root is just outside the spectral value set.  This can help
%           prevent incurring an additional expensive vertical/circular
%           search near convergence, which would be unnecessary.  This 
%           technique also does not affect numerical accuracy with respect
%           to the relative tolerance parameter provided in opts.
%
%       .rand_pushes
%           Only nonzero if opts.fast_search is enabled.  This is the total
%           number of times the new method of avoiding (nearly) singular
%           pencils and/or problematic interior searches succeeded in
%           increasing the current value of eta, in order to push past this
%           problematic areas.
%
%   See also specValSetOptions.
%
%   
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   specValSet.m introduced in ROSTAPACK Version 2.0.
%
% =========================================================================
% |  ROSTAPACK: RObust STAbility PACKage                                  |
% |  Copyright (C) 2014-2019 Tim Mitchell                                 |
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

    % Check all input arguments
    switch nargin 
        case {2,3}  
            % Only A matrix so set other matrices to [] for pseudospectra
            [B,C,D,E]   = deal([]);
            epsilon     = varargin{1};   
            if nargin == 2
                opts    = [];
            else
                opts    = varargin{end};
            end
        case {6,7}  
            % A,B,C,D,E provided
            [B,C,D,E]   = deal(varargin{1:4});
            epsilon     = varargin{5};       
            if nargin == 6
                opts    = [];
            else
                opts    = varargin{end};
            end
        otherwise
            error('specValSet: invalid number of input arguments');
    end
    % B,C,E will be replaced by [] if they are explicit identities.
    % D will be replaced with either [] or sparse(p,m) if it is zero.
    [A,B,C,D,E,D_zero] = processSystemMatrices(A,B,C,D,E);
    assert( isARealNumber(epsilon) && epsilon >= 0, ...
            'specValSet: epsilon must be >= 0'      );
    assert( isempty(opts) || isstruct(opts),        ...
            'specValSet: invalid opts argument'     );
    if ~D_zero && epsilon >= 1/norm(toFull(D))
        error('specValSet: epsilon must be less than 1/norm(D)');
    end
    
    % Grab options 
    opts            = specValSetOptions(opts);
    spec_opts       = getSpectrumOptions(opts);
    discrete        = opts.discrete_time;
    warm_start      = opts.warm_start;
    maxit           = opts.maxit;
    tol             = opts.tol;
    fast_search     = opts.fast_search;
    vertical_first  = opts.vertical_search_first;
    root_order      = opts.root_order; 
    rand_dirs       = opts.random_directions;
    safeguard_width = opts.safeguard_width;
    plotting        = opts.plotting;
    suppress_warns  = opts.suppress_warnings;
    
    % The internal computations, e.g. the LUs and backsolves used in the
    % norm of transfer function evaluations, can produce warnings,
    % particularl when done near eigenvalues.  However, as this code is
    % designed to cope with causes of these warnings, there is typically no
    % need to see them.  As such, by default, we disable all warnings and
    % then restore the warning state after the computation is finished.
    if suppress_warns
        w_state     = warning();
        warning('off');
    end
    
    info_requested  = nargout > 2;
    n_perts         = 0;    
    rand_pushes     = 0;    
    if info_requested 
        [countRootSteps,rootTotals]     = rootCounters();
        [countSearches,searchTotals]    = hrSearchCounters(maxit);
    end
    
    if discrete
        cplx2EtaLoc = @cplx2polar;
        toCplx      = @polar2cplx;
    else
        cplx2EtaLoc = @cplx2cart;
        toCplx      = @cart2cplx;
    end
    
    % If the problem is symmetric, we can restrict all searches to the
    % upper half-plane. Checks if all matrices are real and converts any to
    % bonafide real-valued matrix structures if any have an imaginary part
    % that is physically present but nonetheless zero.
    [is_real,A,B,C,D,E] = isRealSystem(A,B,C,D,E);  
    
    % Compute spectrum to find initial starting eigenvalue.
    [d0,~,d_co,d_u,inf_co,inf_u] = getEigenvalues(A,B,C,D,E,spec_opts);
    if inf_co || inf_u
        error('E must be invertible');
    end
    [eta,loc]       = cplx2EtaLoc(d0);
    if is_real
        [d0,loc]    = ensureInUpperHalf(d0,loc); 
    end
   
    if epsilon == 0 || isinf(eta) || isempty(d0)
        prepareExit();
    end
    
    if plotting
        arc_res     = 2*pi/1000;   % increase denominator for smoother arcs
        plotFn      = @(z,t,varargin) plot( real(z),imag(z),t,          ...
                                            'LineWidth',2, varargin{:}  );
        plotLarge   = @(z,t) plotFn(z,t,'MarkerSize',15);
        plotSmall   = @(z,t) plotFn(z,t,'MarkerSize',8);  
        figure
        hold on
        plotLarge(d_u,'g.'); % Plot uncontrollable/unobservable eigenvalues
        plotLarge(d_co,'k.');% Plot controllable and observable eigenvalues
        plotLarge(d0,'rsquare');    % Highlight starting point
    end
      
    sub_opts    = getSubOpts();
    eigSolver   = makeEigSolver(A,B,C,D,E,epsilon,discrete);
    [ntfFn,ntfRootFn,ntfCountFn] = makeInsideFns(A,B,C,D,E,epsilon,opts);
    
    if fast_search
        searchFn        = @fastSearch;
        initialSearchFn = @initialSearchFast;
        use_interp      = root_order > 1 && root_order < 2;   
    else
        searchFn        = @eigBasedSearch;
        initialSearchFn = @initialSearchEigBased;
    end
    
    % Do initial horizontal/radial search.  Computing the radius must begin
    % with an initial radial search while it is generally more efficient to
    % begin with a horizontal search when computing the abscissa.
    if discrete || ~vertical_first
        [eta,loc]       = initialSearchFn(); 
    end
        
    % Convergent phase of the criss-cross algorithm
    iters               = 0;
    while iters < maxit
        iters = iters + 1;
        [intervals,mps,fs,dfs,ddfs] = computeIntervals(loc);
                              
        % If there are no intervals, we can quit
        if isempty(intervals) 
            if fast_search && discrete
                [intervals,mps,fs,dfs,ddfs] = getRandomAngles();
                if isempty(intervals)
                    prepareExit();
                    return
                end
                rand_pushes = rand_pushes + 1;
            else
                prepareExit();
                return
            end
        end
        
        if plotting
            plotIntervals(intervals,mps);
        end
        
        [eta_new,loc]   = searchFn(mps,fs,dfs,ddfs);    
        eta_diff        = (eta_new - eta)/abs(eta);
        eta             = eta_new;         
        if eta_diff <= 0  % tol 
            prepareExit();   
            return
        end      
    end
    
    prepareExit();
    return
    
    % Private nested functions
      
    function prepareExit()
        if suppress_warns
            warning(w_state);
        end
        if plotting
            plotSmall(toCplx(eta,loc),'bo')
        end
        if ~info_requested
            return
        end
        
        [n_b,n_r,n_f]   = rootTotals();
        [totals,solves] = searchTotals();
        
        if fast_search
            n_hr_eig    = 0;
        else
            n_hr_eig    = sum(totals);
        end
        
        info = struct(  'iters',        iters,              ...
                        'eig_count',    iters + n_hr_eig,   ...
                        'svd_count',    ntfCountFn(),       ...
                        'brk_steps',    n_b,                ...
                        'root_steps',   n_r,                ...
                        'full_steps',   n_f,                ...
                        'hr_problems',  totals,             ...
                        'hr_solves',    solves,             ...
                        'n_perts',      n_perts,            ...
                        'rand_pushes',  rand_pushes         );
    end

    function [ints,mps,f,df,ddf] = computeIntervals(w_prev)
        
        % This finds all the cross sections along a given vertical line or
        % circular.  One expensive eigenvalue problem is solved.  As the
        % purely imaginary eigenvalues when using eig may not be computed
        % so accurately, a relatively large tolerance is used to whether an
        % eigenvalue is purely imaginary or not.
        
        eHSP        = eigSolver(eta);
        ints        = formAllIntervals(eHSP,sub_opts);
        if isempty(ints)
            mps     = [];
            f       = [];
            df      = [];
            ddf     = [];
            return
        end
        ints                = splitAnInterval(ints,w_prev,safeguard_width);
        mps                 = cellfun(@mean,ints);
        if is_real 
            indx            = mps < 0;
            mps(indx)       = [];
            ints(indx)      = [];
            if discrete
                indx        = mps > pi;
                mps(indx)   = [];
                ints(indx)  = [];
            end
        end
        [mps,f,df,ddf,indx] = evaluatePoints(eta,mps);
        ints                = ints(indx);
    end

    function [mps,f,df,ddf,indx] = evaluatePoints(eta,mps)
        % Evaluate each point: eta+1i*mps(j) or eta*exp(1i*mps(j)) 
        [f,df,ddf]  = arrayfun(@(w) ntfFn(eta,w),mps);
        
        % Only keep points that are inside the spectral value set
        indx        = f > 0;
        mps         = mps(indx);
        f           = f(indx);
        df          = df(indx);
        ddf         = ddf(indx);
    end

    function [e_max,l_max] = eigBasedSearch(mps,varargin)  
        
        % The original criss-criss method of doing all the horizontal
        % or radial searches every iteration, using expensive eigenvalue
        % problems.  This cannot only be much more costly but also less
        % accurate.
        e_max = eta;
        l_max = []; 
        total = length(mps);

        for j = 1:total
            l = mps(j);
            [e,l] = getExtremalPoint(A,B,C,D,E,epsilon,eta,l,sub_opts);  
            if is_real
                l = abs(l);
            end
            if e > e_max
                e_max = e;
                l_max = l;
            end
            if plotting
                % Note, in the radius case, if the extremal point that
                % happens to be found is actually in the opposite direction
                % of the given angle, the radial search will not be drawn
                % from the midpoint of this arc cross section but from its
                % opposite point, i.e. e*exp(1i*(l + pi)).
                % However, if the problem has real symmetry, the directions
                % will be drawn in the upper half plane.
                plotSearch(eta,e,l)
            end
        end
        
        if info_requested
            countSearches(total,total);
        end
    end

    function [r,f,delta] = getRoot(l,e0,f0,df0,step0)
        % Find a root > e0 along horizontal/radial line given by l
        [r,f,delta,n_b,n_r,n_f] = findRoot(@(x) ntfFn(x,l),         ...
                                            @(x) ntfRootFn(x,l),    ...
                                            e0,f0,df0,step0,        ...
                                            use_interp,tol          );
        if info_requested
            countRootSteps(n_b,n_r,n_f);
        end
    end

    function [e,l] = fastSearchGivenOrder(mps,f0,df0,step0) 
        
        % This performs one iteration of the horizontal/radial phase, using
        % the new root-finding techniques and intelligently ordering.  It
        % is assumed that the first search is from a point that is indeed
        % inside the spectral value set.  This must not be violated.
        solves      = 1;
        l           = mps(1);
        total       = length(mps);
        [e,f,delta] = getRoot(l,eta,f0,df0,step0);
        if plotting
            plotSearch(eta,e,l);
        end
        
        % Do the remaining searches, where each one is warm started by the
        % increasing value of e, which is the evolving update of eta.  As
        % such, these subsequent searches may no longer be relevant, i.e.
        % the (warm) starting point could be outside the set.
        for j = 2:total
            mp              = mps(j);
            [f,df,ddf]      = ntfFn(e,mp);
            if f > 0 
                e0          = e;
                step0       = rootStep(f,df,ddf);
                [e,f,delta] = getRoot(mp,e,f,df,step0);
                l           = mp;
                solves      = solves + 1;
                if plotting
                    plotSearch(e0,e,l);
                end
            end
        end
        
        % Finally, if the rightmost root encountered is inside the spectral
        % value set, slightly perturb it so that it is outside.  This can
        % help prevent an additional vertical/circular search from being
        % incurred unnecessarily on convergence.
        if f > 0
            delta  = abs(delta);
            if delta/abs(e) <= tol
                delta   = abs(e)*tol;
            end
            k           = 0;
            while f > 0
                k       = k + 1;
                f       = ntfRootFn(e + k*delta,l);
            end
            e           = e + k*delta;
            n_perts     = n_perts + k;
        end
        
        if info_requested
            countSearches(total,solves);
        end
    end

    function [e,l] = fastSearch(mps,f,df,ddf) 
        [mps,f0,df0,step0]  = prioritize(mps,f,df,ddf);
        [e,l]               = fastSearchGivenOrder(mps,f0,df0,step0); 
    end

    function [x,f0,df0,step0] = prioritize(x,f,df,ddf)
        % Sort in order of largest initial Newton/Halley step first
        steps           = arrayfun(@rootStep,f,df,ddf);
        [steps,indx]    = sort(steps,'descend');
        x               = x(indx);
        f               = f(indx);
        df              = df(indx);
        % ddf             = ddf(indx);
        f0              = f(1);
        df0             = df(1);
        step0           = steps(1);
    end

    function [e,l] = initialSearchFast()
        
        % The basic initialization for computing the spectral value set
        % abscissa is to find a root (boundary point of the spectral
        % value set) further to the right of an eigenvalue attaining the
        % spectral abscissa.  For the radius case, a root is found further
        % outward from an eigenvalue attaing the spectral radius.  
        % 
        % Note that the initial value of eta is either the spectra abscissa
        % or spectral radius of (A,E) (modulo infinite, uncontrollable,
        % unobservable eigenvalues, as preferred by the user).
        %
        % However, for real-valued problems, if these initial searches are
        % not along the real axis, complex-valued arithmetic will be
        % incurred when evaluating the norm of the transfer function; not
        % only is each evaluation then significantly more expensive, but
        % finding this initial root may also require many such evaluations
        % (recall that the norm of the transfer function is infinity at the
        % starting point, an eigenvalue of (A,E)).  Thus, if possible, we
        % first attempt to increase eta by doing searches along the real
        % axis, which only incurs real-valued arithmetic, and then do the
        % original search from the hopefully larger value of eta.  This can
        % greatly reduce the number of evaluations needed for the original
        % search and can thus significantly reduce the total computation
        % time.  For the abscissa case, only the positive real axis
        % direction is searched (i.e. y=0); in the radius case, both are
        % (i.e. angles 0 and pi).
        % 
        % The user may also provide a point to warm start the computation.
        % In order to be usable (and still guarantee global convergence),
        % the point must:
        %   (a) be inside the spectral value set and 
        %   (b) its real part (or modulus for the radius case) must be at 
        %       least eta.
        % The norm of the transfer function must be evaluated at this warm
        % start point to check these conditions.  If (a) does not hold, the
        % warm start is ignored.  If (a) holds but (b) does not, the
        % horizontal/radial search direction determined by the warm start
        % can at least be added to the initial searches; even though it may
        % not be helpful, there is little added cost to trying it.  If (a)
        % and (b) both hold, then eta is updated and the search direction
        % for the warm start is added to the initial searches.  The routine
        % must take care to still check the real-axis direction(s) first
        % (for efficiency) and to also ensure that the warm start direction
        % does not duplicate the spectral abscissa/radius direction or the
        % real-axis direction(s).
        % 
        % The following code implements the above algorithm and makes sure
        % that norm of the transfer function is not computed more than once 
        % at any given point.
        
        lx              = [];       % Stores the real-axis directions
        l               = loc;      % Other directions 
        
        % Note that either lx or l may, as the routine proceeds, end up as
        % empty.  As such, care must be taken when using these variables in
        % tests (e.g. ismember below) or when passing them to functions
        % that assume nonempty inputs.  Even though (eta,loc) is an
        % eigenvalue of (A,E) and is always included as an initial point,
        % it may be that it *numerically* evaluates to outside the set and
        % so it may be that both lx and l end up empty.
        
        if is_real
            if discrete
                lx      = [0 pi];
            else
                lx      = 0;
            end
            if ismember(l,lx)
                l       = [];
            end
        end
        
        if isempty(warm_start)
            % No warm start, so do normal initial search
            [lx,fx,dfx,ddfx]        = evaluatePoints(eta,lx);
            [l,f,df,ddf]            = evaluatePoints(eta,l);
        else
            % Warm start provided, need to check if it is usable
            [ew,lw]                 = cplx2EtaLoc(warm_start);
            if is_real 
                [warm_start,lw]     = ensureInUpperHalf(warm_start,lw);
            end
            
            [lw,fw,dfw,ddfw]        = evaluatePoints(ew,lw);
            if isempty(lw) 
                % Warm start is outside the spectral value set so ignore it
                [lx,fx,dfx,ddfx]    = evaluatePoints(eta,lx);
                [l,f,df,ddf]        = evaluatePoints(eta,l);
            else
                % Warm start is in the spectral value set so we can use it
                if plotting
                    plotSmall(warm_start,'bsquare');
                end
                if ew < eta   
                    % Warm start is worse than spec. abs/rad.  The most we
                    % can do is add this warm start direction to the
                    % initial search directions, if it isn't already there
                    if ~ismember(lw,[lx l])
                        l       = [ l lw ];
                    end
                    % Need to do evaluations of the lx and l directions
                    [lx,fx,dfx,ddfx]            = evaluatePoints(eta,lx); 
                    [l,f,df,ddf]                = evaluatePoints(eta,l);
                else
                    % Warm start is at least as good as spec. abs/rad so we
                    % should update eta and add the computation for the
                    % warm start point, i.e. we don't want to recompute it
                    eta         = ew;
                    if ismember(lw,lx) 
                        % Warm start direction is along the real axis
                        [l,f,df,ddf]            = evaluatePoints(eta,l);   
                        if numel(lx) > 1
                            % Only have one of two real-axis directions so
                            % evaluate the other and then merge
                            lx                  = lx(lx ~= lw); 
                            [lx,fx,dfx,ddfx]    = evaluatePoints(eta,lx);   
                            lx                  = [lx lw];
                            fx                  = [fx fw];
                            dfx                 = [dfx dfw];
                            ddfx                = [ddfx ddfw];                     
                        else
                            % Warm start is the real-axis direction
                            lx                  = lw;
                            fx                  = fw;
                            dfx                 = dfw;
                            ddfx                = ddfw;    
                        end     
                    else
                        % Warm start direction is NOT along the real axis
                        [lx,fx,dfx,ddfx]        = evaluatePoints(eta,lx); 
                        if ismember(lw,l)
                            % Warm start direction is same as spec. abs/rad
                            l                   = lw;
                            f                   = fw;
                            df                  = dfw;
                            ddf                 = ddfw;
                        else
                            % Warm start direction is different
                            [l,f,df,ddf]        = evaluatePoints(eta,l);
                            l                   = [l lw];
                            f                   = [f fw];
                            df                  = [df dfw];
                            ddf                 = [ddf ddfw]; 
                        end
                    end
                end     
            end
        end
        
        % Both can empty (e.g. d0 is uncontrollable/unobservable)
        
        % We need to process l before lx, since if lx is nonempty, the last
        % update to f,df,step will set it to the lx step for the initial
        % real-axis search, and we want the real-axis searches first.
        if ~isempty(l)
            [l,f,df,step]   = prioritize(l,f,df,ddf);
        end
        if ~isempty(lx)
            [lx,f,df,step]  = prioritize(lx,fx,dfx,ddfx);
        end
        l                   = [lx l];
        if isempty(l)
            e               = eta;
            l               = loc;
            return
        end
        [e,l]               = fastSearchGivenOrder(l,f,df,step);  
    end

    function [e,l] = initialSearchEigBased()
        e = eigBasedSearch(loc);
        if discrete 
            e = avoidInteriorOld(e,loc);
        end
        l = loc;
    end

    function e = avoidInteriorOld(e,l)
        % This implements the older technique of Mengi and Overton for
        % avoiding singular pencils / problematic interior searches in the
        % radius case, with a modification to ensure the perturbation is at
        % least large enough to make a relative change beyond machine
        % precision.
        fn          = @(e) ntfRootFn(e,l);
        [f,df]      = fn(e);
        delta       = newtonStep(f,df);
        % Make a tolerance a bit above machine precision (1e-15 for 64 bit)
        tol_mach    = 10^(ceil(log10(eps)));     
        if delta/e <= tol_mach
            delta   = e*tol_mach;
        end   
        delta       = abs(delta);
        k           = 0;
        while f > 0
            k       = k + 1;
            f       = fn(e + k*delta);
        end
        e           = e + k*delta;
    end

    function [ints,angles,f,df,ddf] = getRandomAngles()
        limit                   = (1 + ~is_real)*pi;
        angles                  = limit*rand(1,rand_dirs);
        [ints,angles]           = makeIntervals(angles,0,limit);
        [angles,f,df,ddf,indx]  = evaluatePoints(eta,angles);
        ints                    = ints(indx);
    end

    function opts = getSubOpts()
        opts.discrete_time  = discrete;  
        nEA                 = norm(toFull(A));
        if ~isempty(E)
            nEA             = cond(toFull(E))*nEA;
        end
        opts.ham_symp_tol   = sqrt(eps)*max(nEA,epsilon);
        opts.is_symmetric   = is_real;
    end
 
    function plotIntervals(intervals,freqs)
        % Plots the cross section for the vertical/circular searches
        for j = 1:length(intervals)
            cs      = intervals{j};
            plotSmall(toCplx(eta,cs),'ro');         % Endpoints
            plotSmall(toCplx(eta,freqs(j)),'rx');   % Midpoint
            if discrete
                k   = ceil((cs(2) - cs(1))/arc_res) + 2;
                cs  = linspace(cs(1),cs(2),k);
            end
            plotFn(toCplx(eta,cs),'r-');            % Joining segment/arc
        end
    end

    function plotSearch(e0,e,l)
        % Plots the horizontal/radial searches
        pts = toCplx([eta e],l);
        plotFn(pts,'k--');              % Horizontal/radial search line
        plotSmall(toCplx(e0,l),'ko');   % Actual start point of search
        plotSmall(pts(end),'ro')        % Endpoint of search
    end

end

function o = getSpectrumOptions(opts)
    o.discrete_time         = opts.discrete_time;
    o.ignore_infinite       = false;
    o.ignore_unperturbable  = opts.ignore_unperturbable;
    o.find_unperturbable    = opts.ignore_unperturbable || opts.plotting;
end

function [d,y] = ensureInUpperHalf(d,y)
    
    % Given complex scalar d, with either imaginary part of angle y, ensure
    % that d is in the upper half-plane and y is positive.
    
    if y < 0
        y  = abs(y);
        d  = conj(d);
    end
end

function [intervals,values] = makeIntervals(values,lb,ub)
 
    % Make intervals with midpoints in values, such that intervals have lb
    % at the smallest endpoint and ub as the largest endpoint.
    
    values      = sort(values,'ascend');
    midpoints   = 0.5*(values(1:end-1) + values(2:end));
    dividers    = [ lb midpoints ub ];
    intervals   = num2cell([dividers(1:end-1);dividers(2:end)]',2)';
end

function [updateFn,totalsFn] = rootCounters()

    % Maintain counts of the total number of steps taken to bracket roots
    % and then, given a bracket, to converge to a root.  The total number
    % of non-bisection steps in the root phases is also collected.
    
    n_brk   = 0;
    n_root  = 0;
    n_full  = 0;
    
    updateFn    = @updateCount;
    totalsFn    = @getTotals;
    
    function updateCount(n_b,n_r,n_f)
        n_brk   = n_brk + n_b;
        n_root  = n_root + n_r;
        n_full  = n_full + n_f;
    end

    function [n_b,n_r,n_f] = getTotals()
        n_b     = n_brk;
        n_r     = n_root;
        n_f     = n_full;
    end
end

function [updateFn,totalsFn] = hrSearchCounters(maxit)

    % Maintain counts of the number of horizontal and radial searches
    % per iteration and the number that were actually solved per iteration.

    index       = 0;
    totals      = zeros(maxit,1);
    solves      = zeros(maxit,1);
    
    updateFn    = @updateCount;
    totalsFn    = @getTotals;
    
    function updateCount(total,solved)
        index           = index + 1;
        totals(index)   = total;
        solves(index)   = solved;  
    end
    
    function [t,s] = getTotals()
        t   = totals(1:index);
        s   = solves(1:index);
    end

end