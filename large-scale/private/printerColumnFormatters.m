function s = printerColumnFormatters(opts)
%   printerColumnFormatters:
%       Sets up formatters for each column (or group of columns) needed by
%       the different printers implemented for various ROSTAPACK routines:
%       - printerContract
%       - printerExpand
%       - printerInit
%       - printerMain
%       Based on the parameters for the routines, the columns will be set
%       to appropriate widths, etc.
%
%   INPUT:
%       opts    A struct of the required/optional parameters:
%       .maxit          [positive integer, required]
%           The maximum possible number of iterations.
%       .phase_maxit    [positive integer, optional]
%           The maximum possible number of iterations in the subroutines
%           (or phases).
%       .step_type      [logical, optional]
%           A logical whether or not a step type column is requested.
%       .ls_maxit       [positive integer, optional]
%           The maximum possible number of line search iterations
%       .kth_eigenvalue [positive integer, optional]
%           The requested eigenvalue to use.
%      
%   OUTPUT:
%       A struct containing formatters for the following columns/sections:
%       .iter
%       .main
%       .phase          (if opts.phase_maxit is present)
%       .step_type      (if opts.step_type is present and true)
%       .line_search    (if opts.ls_maxit is present)
%       .init_pert      (if opts.kth_eigenvalue is present)
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   printerColumnFormatters.m introduced in ROSTAPACK Version 1.0
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

    s.iter              = getCount('Iter',opts.maxit);
    s.main              = getEpsilonAndF();
    
    if isfield(opts,'phase_maxit')
        s.phase         = getPhase(opts.phase_maxit);
    end
    
    if isfield(opts,'step_type') && opts.step_type
        s.step_type     = getCentered('Step\nType',4);
    end
    
    if isfield(opts,'ls_maxit') 
        s.line_search   = getLineSearch(opts.ls_maxit);
    end
    
    if isfield(opts,'kth_eigenvalue')
        s.init_pert     =  getInitMatrices(opts.kth_eigenvalue);
    end
end

function iter = getCount(label,maxit)
    width = max(nDigitsInWholePart(maxit),length(label));
    iter = struct(                                                      ...
        'label',            label,                                      ...
        'width',            width,                                      ...
        'format_fn',        @(x) sprintf('%*s',width,num2str(x)),       ...
        'blank_str',        blanks(width),                              ...
        'na_str',           centerString('-',width)                     );
end

function type = getCentered(label,width)
    type = struct(                                                      ...
        'label',                label,                                  ...
        'width',                width,                                  ...
        'format_fn',            @(x) centerString(x,width),             ...
        'na_str',               centerString('-',width)                 );
end

function s = getEpsilonAndF()

    persistent format_struct;
    
    if ~isempty(format_struct)
        s = format_struct;
        return
    end
    
    num_width   = 23;   % Always use 23 chars to display epsilon and f 
    diff_width  = 9;    % Always use 9 chars to display the differences

    get_num_fn  = double2FixedWidthStr(num_width);
    get_diff_fn = double2FixedWidthStr(diff_width);
    
    % allow for () around the numbers for rejected iterates
    % but epsilon is always nonnegative so we don't need a space for +/-
    eps_width = num_width + 1;  
    epsilon = struct(                                                   ...
        'label',                'Value',                                ...
        'width',                eps_width,                              ...
        'format_fn',            @(x) getEpsilonAccepted(x,get_num_fn),  ...
        'format_rejected_fn',   @(x) getEpsilonRejected(x,get_num_fn),  ...
        'blank_str',            blanks(eps_width)                       );
        
    % allow for () around the numbers for rejected iterates
    f = struct(                                                         ...
        'label',                'Value',                                ...
        'width',                num_width + 2,                          ...
        'format_fn',            @(x) getFAccepted(x,get_num_fn),        ...
        'format_rejected_fn',   @(x) getFRejected(x,get_num_fn)         );
    
    diff = struct(                                                      ...
        'label',                'Diff',                                 ...
        'width',                diff_width,                             ...
        'format_fn',            get_diff_fn,                            ...
        'blank_str',            blanks(diff_width),                     ...
        'na_str',               centerString('-',diff_width)            );
    
    format_struct   = struct('epsilon',epsilon,'f',f,'diff',diff);
    s               = format_struct;
end

function s = getPhase(maxit)
    s = struct( 'type',     getCentered('Type',4),                      ...
                'subiters', getCount('Iters',maxit),                    ...
                'term',     getCentered('TC',2)                         );
end

function s = getLineSearch(maxit)

    ls_itr_width    = max(nDigitsInWholePart(maxit),5);
    LS_T_WIDTH      = 9;
    get_ls_t        = double2FixedWidthStr(LS_T_WIDTH);
    
    iters = struct(                                                     ...
        'label',                '#',                                    ...
        'width',                ls_itr_width,                           ...
        'format_fn',            @(x) sprintf('%*d',ls_itr_width,x),     ...
        'na_str',               centerString('-',ls_itr_width)          );
    
    t_value = struct(                                                   ...
        'label',                't',                                    ...
        'width',                LS_T_WIDTH,                             ...
        'format_fn',            get_ls_t,                               ...
        'na_str',               centerString('-',LS_T_WIDTH)            );
    
    s = struct('iters',iters,'t',t_value);
end

function s = getInitMatrices(kth_eigenvalue)  
    s = struct(                                                         ...
        'type',                 getCentered('Type',4),                  ...
        'kth_eval',             getCount('Eig. #',kth_eigenvalue)       );
end

function s = getEpsilonAccepted(x,format_fn)
    s = sprintf('%s ',format_fn(x));
end

function s = getEpsilonRejected(x,format_fn)
    epsil_str   = format_fn(x);
    s           = sprintf('(%s)',epsil_str(2:end));
end

function s = getFAccepted(x,format_fn)
    s = sprintf(' %s ',format_fn(x));
end

function s = getFRejected(x,format_fn)
    s = sprintf('(%s)',format_fn(x));
end
