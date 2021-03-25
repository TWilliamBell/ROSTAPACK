function [d0,d_fin,d_co,d_u,inf_co,inf_u] = getEigenvalues(A,B,C,D,E,opts) 
%   getEigenvalues:
%       Computes the eigenvalues of (A,E) and returns a globablly rightmost
%       or outermost eigenvalue, subject to the parameters in opts.  The
%       entire (finite) spectrum is also returned, sorted by decreasing
%       real part or decreasing modulus, and optionally, as sorted subsets
%       of those that are both controllable and observable and those that
%       either uncontrollable or unobservable.
% 
%   INPUT:
%       A,B,C,D,E                   [matrices]
%           System matrices.  A must be provided while B,C,D,E can be given
%           explicitly or as [] for their shortcuts.
% 
%       opts                        [struct]
%       .discrete_time              [logical]
%           False:  sorts eigenvalues by decreasing real part
%           True:   sorts eigenvalues by decreasing modulus
%
%       .ignore_infinite            [logical]
%           When this is true, all infinite eigenvalues are excluded from
%           consideration as a globally rightmost/outermost eigenvalue. 
% 
%       .ignore_unperturbable       [logical]
%           When this is true, all unperturbable eigenvalues (either
%           uncontrollable or unobservable) are excluded from consideration
%           as a globally rightmost/outermost eigenvalue.  Note that 
%           additional computations are necessary to determine this.  
%           
%       .find_unperturbable         [optional, logical]
%           If enabled, this allows the unperturbable eigenvalues to still
%           be determined, even if opts.ignore_unperturbable is set to
%           false.  This may be desirable if one still wishes to know which
%           eigenvalues are perturbable and which are not, even though one
%           wants the rightmost/outermost eigenvalue regardless of its
%           perturbability or lackthereof.
%   
%   OUTPUT:
%       d0                          [complex scalar]
%           A globablly rightmost/outermost eigenvalue, as determined from
%           opts.ignore_unperturbable and opts.ignore_infinite.  If this
%           eigenvalue if infinite, it is returned as just inf, regardless
%           of whether it was actually +-inf, +-i1*inf, or +-inf +- 1i*inf.
%           
%       d_fin                       [vector of complex scalars]
%           All finite eigenvalues, sorted
%
%       d_co                        [vector of complex scalars]
%           All controllable and observable finite eigenvalues, sorted
%           (Equal to d_all when opts.ignore_unperturtable is false.)
%
%       d_u                         [vector of complex scalars]
%           All uncontrollable/unobservable finite eigenvalues, sorted
%           (If opts.ignore_unperturtable and opts.find_unperturbable are
%           both false, d_u is returned as [].)
% 
%       inf_co                      [logical]
%           True if any infinite eigenvalue of (A,E) is both controllable
%           and observable.            
%           (If opts.ignore_unperturtable and opts.find_unperturbable are
%           false, inf_co is true if (A,E) has any infinite eigenvalues,
%           regardless of whether they controllable and/or observable.)
%   
%       inf_u                       [logical]
%           True if any infinite eigenvalue of (A,E) is either 
%           uncontrollable and/or unobservable.
%           (If opts.ignore_unperturtable and opts.find_unperturbable are
%           false, inf_u is always false.)
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   getEigenvalues.m introduced in ROSTAPACK Version 2.0
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
    
    ignore_inf      = opts.ignore_infinite;
    ignore_unpert   = opts.ignore_unperturbable;
    
    if opts.discrete_time
        xFn     = @abs;
        yFn     = @angle;
    else
        xFn     = @real;
        yFn     = @imag;
    end
    
    [n,m,p]     = systemDimensions(A,B,C,D);
    if ~isempty(B) && isEye(B)
        B       = [];
    end
    if ~isempty(C) && isEye(C)
        C       = [];
    end
    
    if m == n && p == n && isempty(B) && isempty(C)
        % If B=C=I and m=p=n, all eigenvalues are controllable & observable
        find_co = false;
    elseif  ignore_unpert
        % The unperturbable eigenvalues can only be ignored if we do the
        % additional calculations to separate them out.
        find_co = true;
    elseif isfield(opts,'find_unperturbable')
        find_co = opts.find_unperturbable;
    end
                 
    A           = toFull(A);
    if isempty(E) || isEye(E)
        args    = {A};
    else
        E       = toFull(E);
        args    = {A,E};
    end
       
    if find_co
        [d,co]  = getEigenvaluesCO(args{:});  
    else
        [d,co]  = getEigenvaluesAll(args{:});
    end 
    fin         = ~isinf(d);
    inf_co      = any(~fin & co);
    inf_u       = any(~fin & ~co);
    % Get rid of infinite eigenvalues before sorting
    d_fin       = d(fin);
    co          = co(fin);
    [d_fin,co]  = sortEigenvalues(d_fin,co);
    
    % Controllable and observable eigenvalues
    d_co        = d_fin(co);
    % Uncontrollable or unobservable eigenvalues
    d_u         = d_fin(~co);
    
    % Rightmost/outermost eigenvalue
    if ~ignore_inf && (inf_co || (inf_u && ~ignore_unpert))
        d0      = inf;
    elseif ignore_unpert && ~isempty(d_co)
        d0      = d_co(1);
    elseif ~isempty(d_fin)
        d0      = d_fin(1);
    else
        d0      = [];
    end

    % Private functions 
        
    function [d,co] = getEigenvaluesAll(varargin)
        % Just assume all eigenvalues controllable and observable
        d       = eig(args{:});
        co      = true(length(d),1);
    end

    function [d,co] = getEigenvaluesCO(varargin)
        % Calculating which are controllable and observable also requires
        % the left eigenvectors.
        [X,d,Y] = eig(args{:},'vector');
        
        % It is critical to use vecnorm(M,2,1) since either B'*Y or C*X
        % could be a vector.  In this case, vecnorm instead returns the
        % norm of the (row or column) vector but what we need are the norms
        % of the columns of B'*Y and C*X.
        if isempty(B)
            BhY = vecnorm(Y(1:m,:),2,1);
        else
            BhY = vecnorm(B'*Y,2,1);
        end
        if isempty(C)
            CX  = vecnorm(X(1:p,:),2,1);
        else
            CX  = vecnorm(C*X,2,1);
        end
        % Indices of controllable and observable eigenvalues, as a column
        co      = (BhY > 0 & CX > 0).';
    end

    function [d,y] = sortEigenvalues(d,y) 
        [d,y]   = sortByX(yFn,d,y);
        [d,y]   = sortByX(xFn,d,y);
    end

    function [x,y] = sortByX(orderFn,x,y)
        [~,idx] = sort(orderFn(x),'descend');
        x       = x(idx);
        y       = y(idx);
    end

end


