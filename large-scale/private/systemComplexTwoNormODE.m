function obj = systemComplexTwoNormODE(sys,varargin)
%   systemComplexTwoNormODE:
%       This object manages all spectral-value-set-specific computations,
%       keeps track of state, and supports snap shots.  The algorithmic
%       routines in this folder of ROSTAPACK are abstracted so that they
%       can operate over any kind of spectral value set (a particular norm,
%       constraints on the kind of perturbations allowed, e.g. real or
%       complex), provided that a corresponding system[Type] object is
%       implemented.
%
%       systemComplexTwoNormODE is an implementation for regular
%       complex-valued spectral value sets, using an ODE-based update step
%       to approximate the spectral value set abscissa and/or radius, which
%       is used as a subroutine in the Hybrid Expansion-Contraction (HEC).
%
%       To support a specific type of spectral value set, one must
%       implement at least the following four functions:
%
%           [df_svs,U1,V1] = computeDfSVSAndUVDzero() [nested]
%           [df_svs,U1,V1] = computeDfSVSAndUV() [nested]
%           psi = computePsiDzero(U1,V1) [nested]
%           psi = computePsi(U1,V1) [nested]
%
%       which provide the following values:
%
%           df_svs
%               The derivative of the spectral value set abscissa|radius
%               with respect to epsilon. 
%
%           U1,V1
%               Given a current perturbation U0V0', U1V1' is the
%               perturbation that pushes an eigenvalue rightward/outward.
%
%           psi
%               This is the derivative of the line search function with
%               respect to t at t = 0 on the unit norm interpolated path
%               between U0V0' and U1V1'.  
%
%       These four functions are nested so that all values are easily
%       accessible and the function signatures can be minimal in terms of
%       their number of input arguments. You may add your own nested or
%       local helper functions as desired.
%
%       You may also need to modify the following three functions
%       accordingly:
%
%           obj = makeSystemObject() [nested]
%           s = snapShot() [nested]
%           restoreSnapShot() [nested]
%
%       in order to:
%           - specify a particular normalization routine for U*V' 
%           - update snapshots for new variables if you need them
%           - add support for additional subsolver that is required.
%             (For an example, see phi_solver in systemComplexTwoNorm.)
%
%       NOTE: when implementing a new type, much of this file of code
%       should NOT change at all.  Inheritance could have been used to
%       avoid this duplication of shared code across different spectral
%       value set types, but there are a couple reasons for not doing this.
%       First, it is difficult to predict how much code will be able to be
%       shared for new spectral value set types.  For example, the routine
%       for computing the derivative necessary for the contraction phase of
%       HEC (computeDerivativeEpsilon) can be shared across complex and
%       real spectral value sets but it is difficult to predict if a
%       modified routine will be necessary for other types of spectral
%       value sets.  Second, it was decided to avoid using MATLAB's OOP and
%       use structs with handles of nested functions to make objects
%       instead.
%
%   See also systemComplexTwoNorm and systemRealFrobeniusNorm.
%
%
%   For more details, see [GGO13], [MO16], and [GGMO17, Appendix B], the
%   last of which gives the formula for the perturbation update. 
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   systemComplexTwoNormODE.m introduced in ROSTAPACK Version 1.0
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

    MAX_RANK    = 1;

    [n,m,p]     = sys.getDimensions();
    B           = sys.getMatrixObject('B');
    C           = sys.getMatrixObject('C');
    D           = sys.getMatrixObject('D');

    epsilon     = 0;
    epsilon_ub  = D.getNorm()^-1;
  
    U           = zeros(m,MAX_RANK);
    V           = zeros(p,MAX_RANK);
    I           = eye(MAX_RANK);
    
    % These update if new U and V are provided
    BU          = zeros(n,MAX_RANK);
    ChV         = zeros(n,MAX_RANK);
    DU          = zeros(p,MAX_RANK);
    VhD         = zeros(MAX_RANK,m);
    VhDU        = zeros(MAX_RANK);

    % These update if new epsilon or U,V are provided
    I_epsVhDU   = I;
    I_T         = I;
    BU_eI_T     = U;

    % These update any time an eigenvalue is requested
    f           = [];
    z           = [];
    x           = zeros(n,1);

    % These update when derivatives or a new perturbation is requested
    y           = zeros(n,1);
    absyhx      = [];
    Bhy         = zeros(m,1);
    Cx          = zeros(p,1);
    yhBU        = 0;
    VhCx        = 0;
    z_epsilon   = [];
    
    % Derivative of the line search w.r.t. to t at t = 0, between two
    % perturbations U0*V0' and U1*V1', along a path of interpolated but
    % still unit-norm matrices between them.
    df_UV_t0    = [];

    check_point = [];

    if D.isZero()
        dfsvs_and_UV_fn = @computeDfSVSAndUVDzero;
        compute_psi_fn  = @computePsiDzero;
    else
        dfsvs_and_UV_fn = @computeDfSVSAndUV;
        compute_psi_fn  = @computePsi;
    end  
    eig_solver  = sys.getEigSolver();
    root_fn     = sys.getRootFunction();

    obj = makeSystemObject();

    % IMPLEMENT:
    %   - [df_svs,U1,V1] = computeDfSVSAndUVDzero() [nested]
    %   - [df_svs,U1,V1] = computeDfSVSAndUV() [nested]
    %   - psi = computePsiDzero(U1,V1) [nested]
    %   - psi = computePsi(U1,V1) [nested]
    %     
    %   Add your own helper functions as desired.
    %
    % MODIFY AS NECESSARY:
    %   - obj = makeSystemObject()
    %   - s = snapShot()
    %   - restoreSnapShot()

    % This is a helper function for the two computeDfSVSAndUV functions
    function [Bhy,Cx] = getBhYCX(Bhy,Cx)
        % intentionally blank, since this is a pass-through variant for
        % consistency with systemRealFrobeniusNorm. 
    end

    function [df_svs,U1,V1] = computeDfSVSAndUVDzero()                                
        [U1,V1]         = getBhYCX(Bhy,Cx);   
        [U1,V1,beta]    = normalizeTwoNormUV(U1,V1);
        df_svs          = 1/(beta * absyhx);  
    end

    function [df_svs,U1,V1] = computeDfSVSAndUV()  
        [BhY,CX]        = getBhYCX(Bhy,Cx);  
        DhV             = D.applyHermitian(V);
        VhD             = DhV';                
        U1              = BhY + DhV*(epsilon*( I_epsVhDU' \ (U'*BhY)) );
        V1              = CX + DU*(epsilon*( I_epsVhDU \ (V'*CX)) );   
        [U1,V1,beta]    = normalizeTwoNormUV(U1,V1);
        df_svs          = 1/(beta * absyhx);  
    end  

    function psi = computePsiDzero(U1,V1)
        % When D == 0, Ubar == Uhat
        Ubar    = computeUhat(U1,V1,U,V);
        BUbar   = B.apply(Ubar);
        ChV1    = C.applyHermitian(V1);

        psi     = epsilon * ((y'*BUbar)*(ChV'*x) + (y'*BU)*(ChV1'*x));
        psi     = real(psi);
    end

    function psi = computePsi(U1,V1)
        Uhat    = computeUhat(U1,V1,U,V);
        V1hDU   = V1'*DU;
        VhDUhat = VhD*Uhat;

        That    = epsilon * ( V1hDU / I_epsVhDU );
        eps_I_T = epsilon * I_T; 

        Ubar    = Uhat + (U * eps_I_T * VhDUhat);                                       
        Vbar    = ( V1' + ( That * V' ) )';

        BUbar   = B.apply(Ubar);
        ChVbar  = C.applyHermitian(Vbar);

        psi     = (y'*BUbar) * eps_I_T * (ChV'*x) + ...
                  (y'*BU) * (epsilon*(I+eps_I_T*VhDU)) * (ChVbar'*x);       
        psi     = real(psi);
    end

    function obj = makeSystemObject()
        
        % We need to set four things specifically for this type of system:
        normalize_UV_fn         = @normalizeTwoNormUV;
        get_totals_fn           = @eig_solver.getTotals;
        update_phi_epsilon_fn   = @NOP;
        update_phi_v0           = @NOP;
        
        % Keep everything below the same.  This defines and sets all the 
        % of this object's methods.
        
        obj = struct( ...
            'isReal',               sys.isReal,                         ...
            'isDiscreteTime',       sys.isDiscreteTime,                 ...
            'normalizeUV',          normalize_UV_fn,                    ...
            'getEpsilonUpperBound', @() epsilon_ub,                     ...
            'getEpsilon',           @getEpsilon,                        ...
            'getUV',                @getUV,                             ...
            'getf',                 @getf,                              ...
            'getDfEpsilon',         @computeDerivativeEpsilon,          ...
            'getEigenvectors',      @getEigenvectors,                   ...
            'getTotals',            get_totals_fn,                      ...
            'computeNextUV',        @computeNextUV,                     ...
            'computeEigA',          @computeEigA,                       ...
            'computeDeltaFromA',    @computeDeltaFromA,                 ...
            'computeEigEpsilon',    @getEigEpsilon,                     ...
            'computeEigUV',         @getEigUV,                          ...
            'computeEigDelta',      @getEigDelta,                       ...
            'snapShot',             @snapShot,                          ...
            'restoreSnapShot',      @restoreSnapShot,                   ...
            'updateEigsV0',         @eig_solver.updateInitialVectors,   ...
            'updatePhiV0',          update_phi_v0,                      ...
            'updatePhiSolver',      update_phi_epsilon_fn               );
    end

    function s = snapShot()
        
        s.epsilon   = epsilon;
        s.U         = U;
        s.V         = V;

        s.BU        = BU;
        s.ChV       = ChV;
        s.DU        = DU;
        s.VhD       = VhD;
        s.VhDU      = VhDU;
        s.I_epsVhDU = I_epsVhDU;
        s.I_T       = I_T;
        s.BU_eI_T   = BU_eI_T;

        s.f         = f;
        s.z         = z;
        s.x         = x;

        s.y         = y;
        s.absyhx    = absyhx;
        s.Bhy       = Bhy;
        s.Cx        = Cx;
        s.yhBU      = yhBU;
        s.VhCx      = VhCx;
        s.z_epsilon = z_epsilon;

        s.df_UV_t0  = df_UV_t0;

        check_point = s;
    end

    function restoreSnapShot(s)
        if nargin == 0
            if ~isempty(check_point)
                s = check_point;
            else
                return
            end
        end
  
        epsilon     = s.epsilon;
        U           = s.U;
        V           = s.V;
        
        BU          = s.BU;
        ChV         = s.ChV;
        DU          = s.DU;
        VhD         = s.VhD;
        VhDU        = s.VhDU;
        I_epsVhDU   = s.I_epsVhDU;
        I_T         = s.I_T;
        BU_eI_T     = s.BU_eI_T;
        
        f           = s.f;
        z           = s.z;
        x           = s.x;

        y           = s.y;
        absyhx      = s.absyhx;
        Bhy         = s.Bhy;
        Cx          = s.Cx;
        yhBU        = s.yhBU;
        VhCx        = s.VhCx;
        z_epsilon   = s.z_epsilon;
        
        df_UV_t0    = s.df_UV_t0;
    end

    % THE FUNCTIONS BELOW ALL ARE NECESSARY BUT IT IS LIKELY THEY THAT 
    % DO NOT NEED TO BE MODIFIED
    
    function epsilon_out = getEpsilon()
        epsilon_out = epsilon;
    end

    function [U_out,V_out] = getUV()
        U_out = U;
        V_out = V;
    end

    function [f_out,z_out] = getf()
        f_out = f;
        z_out = z;
    end

    function [x_out,y_out,absyhx_out] = getEigenvectors()
        [x_out,y_out,absyhx_out] = deal([]);
        if isempty(z)
            return
        elseif ~isempty(absyhx)
            [x_out,y_out,absyhx_out] = deal(x,y,absyhx);
            return
        end
        if isZero(BU_eI_T) || isZero(ChV)
            getEvecsFn  = @() eig_solver.eigentripleA(z,x);
        else
            getEvecsFn  = @() eig_solver.eigentripleAUV(BU_eI_T,ChV,z,x);
        end
        try 
            [x_out,y_out,absyhx_out] = getEvecsFn();
        catch
            x_out = x;
        end
    end

    function df = computeDerivativeEpsilon()
        if isempty(absyhx)
            [x,y,absyhx]    = eig_solver.eigentripleOfAUV(BU_eI_T,ChV,z,x);
            Bhy             = B.applyHermitian(y);
            Cx              = C.apply(x);
            yhBU            = Bhy'*U;
            VhCx            = V'*Cx;
            z_epsilon       = (yhBU*(I_T^2)*VhCx) / absyhx;
        end
        df = real(z_epsilon);
    end

    function [U1,V1,df_UV,df_eps,df_svs] = computeNextUV()
        df_eps          = computeDerivativeEpsilon();
        [df_svs,U1,V1]  = dfsvs_and_UV_fn();
        psi             = compute_psi_fn(U1,V1);
        df_UV_t0        = real(psi)/absyhx;

        if df_UV_t0 < 0
            U1          = -U1;
            V1          = -V1;
            df_UV_t0    = -df_UV_t0;
        end
        df_UV           = df_UV_t0;
    end

    function [f1,z1,fk,zk,xk] = computeEigA(kth_eval)
        absyhx          = [];
        
        selection       = ternOp(kth_eval > 1, [1 kth_eval], 1);
        [z_vals,X_vals] = eig_solver.eigenpairsOfA(selection);
        f_all           = arrayfun(root_fn, z_vals);
        
        getKthEigFn     = @(k) deal(f_all(k), z_vals(k), X_vals(:,k));
        [f1,z1,x1]      = getKthEigFn(1);
        [fk,zk,xk]      = getKthEigFn(kth_eval);
        [f,z,x]         = deal(f1,z1,x1);
    end

    function [epsilon0,U0,V0] = computeDeltaFromA(f_in,z_in,x_in)
        [x,y,absyhx]    = eig_solver.eigentripleOfA(z_in,x_in);
        % Only update the state if the eigensolve succeeds
        z               = z_in;
        f               = f_in;
        Bhy             = B.applyHermitian(y);
        Cx              = C.apply(x);
        % epsilon = 0 means D has no effect
        [df_svs,U0,V0]  = computeDfSVSAndUVDzero();
        epsilon0        = -f/df_svs; % Newton step
    end

    function [f_out,z_out] = computeEigAUV(kth_eval)
        absyhx  = [];
        [z,x]   = eig_solver.eigenpairsOfAUV(BU_eI_T,ChV,kth_eval);
        f       = root_fn(z);
        f_out   = f;
        z_out   = z;
    end

    function [f,z] = getEigEpsilon(epsilon_new)
        updateEpsilon(epsilon_new);
        [f,z] = computeEigAUV(1);
    end

    function [f,z] = getEigUV(U_new,V_new)
        updateUV(U_new,V_new);
        [f,z] = computeEigAUV(1);
    end

    function [f,z] = getEigDelta(epsilon_new,U_new,V_new,kth_eval)
        if nargin < 4 || kth_eval < 1
            kth_eval = 1;
        end
        updateDelta(epsilon_new,U_new,V_new); 
        [f,z] = computeEigAUV(kth_eval);
    end

    function updateMatricesEpsilonDependent()
        I_epsVhDU   = I - epsilon * VhDU;
        I_T         = I + epsilon * ( VhDU / I_epsVhDU );
        BU_eI_T     = BU*(epsilon*I_T);
    end

    function updateMatricesUVDependent(U_new,V_new)
        U           = U_new;
        V           = V_new;

        BU          = B.apply(U);
        ChV         = C.applyHermitian(V);
        DU          = D.apply(U);
        VhDU        = V'*DU;
    end

    function updateUV(U_new,V_new)
        updateMatricesUVDependent(U_new,V_new);
        updateMatricesEpsilonDependent();
    end

    function updateEpsilon(epsilon_new)
        if epsilon_new ~= epsilon 
            if epsilon_new < 0 || epsilon_new >= epsilon_ub
                error('%s: epsilon_new is out of bounds.',mfilename());
            end
            epsilon = epsilon_new;
            updateMatricesEpsilonDependent();
        end   
    end

    function updateDelta(epsilon_new,U_new,V_new)
        if epsilon_new < 0 || epsilon_new >= epsilon_ub
            error('%s: epsilon_new is out of bounds.',mfilename());
        end
        updateMatricesUVDependent(U_new,V_new);
        epsilon = epsilon_new;
        updateMatricesEpsilonDependent();
    end
end

% This is a helper function for the two computePsi functions
function Uhat = computeUhat(U1,V1,U,V)
    UhU     = U'*U;
    UhU1    = U'*U1;
 
    VhV     = V'*V;
    VhV1    = V'*V1;
    
    eta     = trace( UhU1*VhV + UhU*VhV1 );
    Uhat    = U1 - eta*U;
end
