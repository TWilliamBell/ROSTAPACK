function interpolation_model_fn = lineSearchInterpolation(type)
%   lineSearchInterpolation:
%       Returns a function for generating linear, quadratic, or cubic
%       interpolations for use with line searches.  
% 
%   INPUT:
%       type            [required value in {1,2,3}]
%           1 - request linear interpolation 
%           2 - request quadratic interpolation 
%           3 - request cubic interpolation
%       
%   OUTPUT:
%       interpolation_model_fn
%           A function to a compute an interpolatory model for a line
%           search function f at t = 0:
%       
%               computeNewT = interpolationModelFn(f0,f1,g0)
%   
%           where:
%          
%               f0 = f(0) is its value at t = 0
%               f1 = f(1) is its value at t = 1
%               g0 = f'(0), the derivative of f with respect t at t = 0.
%
%           Note that g0 is ignored when using linear interpolation.
%       
%       Then new interpolation-derived values of t can be computed by
%       successively calling:
%
%           [t_new,m_new,bisected] = quadraticModel(t_k,ft_k,gt_k)
%       
%       where:
%   
%           t_k         = the current value of t (usually t_1 is 1)
%           ft_k        = f(t_k), the value of line search function at t_k
%           gt_k        = f'(t_k), the derivative of f at t_k
% 
%           t_new       = the interpolated value to try for the next line
%                         search evaluation
%           m_new       = value of the interpolation model at t_new
%           bisected    = true if a simple bisection was done.  This is
%                         always for linear interpolation.  However, higher
%                         order interpolations may sometimes fallback to
%                         bisection, due to floating point inaccuracy
%                         issues that may arise.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   lineSearchInterpolation.m introduced in ROSTAPACK Version 1.0
%
% =========================================================================
% |  lineSearchInterpolation.m                                            |
% |  Copyright (C) 2016 Tim Mitchell                                      |
% |                                                                       |
% |  This file is originally from URTM.                                   |
% |                                                                       |
% |  URTM is free software: you can redistribute it and/or modify         |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  URTM is distributed in the hope that it will be useful,              |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================
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

    switch type
        case 1
            interpolation_model_fn = @getBisectionModel;
        case 2
            interpolation_model_fn = @getQuadraticModel;
        case 3
            interpolation_model_fn = @getCubicModel;
        otherwise     
            error('Unknown line search interpolation model type.');
    end
end

function fn = getBisectionModel(f0,f1,varargin) 
    fn = @bisect;
   
    function [t,ft,bisected] = bisect(t_k,varargin)
        t           = 0.5*t_k;
        ft          = t*f1 + (1-t)*f0;
        bisected    = true;
    end
end

function fn = getQuadraticModel(f0,f1,g0)
    bisect_fn       = getBisectionModel(f0,f1);
    fn              = @quadraticModel;
    
    function [t,qt,bisected] = quadraticModel(t1,ft1)
        bisected    = false;
        [t,qt]      = getQuadraticMaximizer(f0,g0,t1,ft1);
        
        if isnan(t) || t <= 0 || t >= t1
            [t,qt,bisected] = bisect_fn(t1);
        end
    end
end

function fn = getCubicModel(f0,f1,g0)

    t2          = [];
    ft2         = [];
    bisect_fn   = getBisectionModel(f0,f1);
    fn          = @cubicModel;
    
    function [t,ct,bisected] = cubicModel(t1,ft1)  
        bisected = false;
        if isempty(t2) % initial model is quadratic
            [t,ct] = getQuadraticMaximizer(f0,g0,t1,ft1);
        else            
            [t,ct] = getCubicMaximizer(f0,g0,t1,ft1,t2,ft2);     
        end
        if isnan(t) || t <= 0 || t >= t1
            [t,ct,bisected] = bisect_fn(t1);
        end
        t2  = t1;
        ft2 = ft1;
    end
end

function [t,qt] = getQuadraticMaximizer(f0,g0,t1,ft1)
    % This will return the maximizer in (0,t1) of a quadratic model
    % 
    %       qm(t) = c_2*t^2 + c_1*t + c0
    % 
    % interpolating (0,f0) and (t1,ft1), with 0 < t1, and such that 
    % 
    %       qm(0)' == g0 > 0.  
    % 
    % Note that: 
    %       c0  = f0
    %       c1  = g0
   
    % form quadratic model with respect to t
    c2      = ((ft1 - f0)/(t1^2)) - (g0/t1);
    qm      = @(t) c2*t.^2 + g0*t + f0;
    
    % maximizer and its corresponding quadratic model maximal value
    t       = -g0/(2*c2);
    qt      = qm(t);
end

function [t,ct] = getCubicMaximizer(f0,g0,t1,ft1,t2,ft2)
    % This will return the maximizer in (0,t1) of a cubic model
    % 
    %       cm(t) = c_3*t^3 + c_2*t^2 + c_1*t + c0
    % 
    % interpolating (0,f0), (t1,ft1), and (t2,ft2), with 0 < t1 < t2, and 
    % such that 
    % 
    %       cm(0)' == g0 > 0.  
    % 
    % Note that: 
    %       c0  = f0
    %       c1  = g0
   
    % compute coefficients c3 and c2 of the cubic model
    frac    = 1/(t1 - t2);
    t0sq    = t2^2;
    t1sq    = t1^2;
    t_mat   = frac * [ 1/t1sq -1/t0sq; -t2/t1sq t1/t0sq];
    c       = t_mat*[ft1 - f0 - g0*t1; ft2 - f0 - g0*t2];    
    c3      = c(1);
    c2      = c(2);
          
    % cubic model with respect to t
    cm      = @(t) c3*t.^3 + c2*t.^2 + g0*t + f0;
   
    % this should be purely real but numerically... well, let's be careful
    t       = (-c2 + real(sqrt(c2^2 - 3*c3*g0))) / (3*c3);
    ct      = cm(t);
end
