function [print_level,printer] = printerConfig(opts,type)
%   printerConfig:
%       This returns a printerContract, printerExpand, or printerInit
%       object, along with the desired print level, for use with
%       epsilonContract, epsilonExpand, uvExpand, or
%       initializePerturbation.
%      
%   INPUT:
%       opts                [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine, which determine the exact configuration of the printer
%           of the specified type:
%
%           .maxit      
%               Required when requesting either a printerContract or
%               printerExpand or printer.  This can also be optionally
%               supplied when requesting a printerInit printer, in order to
%               match the width of its iters column to another printer.
% 
%           .line_search_opts.maxit 
%               Required when requesting a printerExpand printer.
% 
%           .initial_perturbation.kth_eigenvalue.
%               Required when requesting a printerInit printer.
%
%           .print_level              
%               The desired print level.
% 
%           .print_ascii            
%               Optional.  Set to true to fallback to standard ASCII
%               character set for printing table borders.
%
%           .printer
%               If a printer[Type] object is already configured, it can
%               be provided here, to prevent instantiating a new one.  
%               Warning: it is assumed that the preconfigured printer is
%               compatible with the parameters stored in opts. 
%
%       type                [required: string]
%           One of {'contract','expand','init'} to choose which printer is
%           desired.  Note it is not possible to autodetect the printer
%           type from which fields are present in opts, since .maxit can be
%           optionally provided for printerInit.
%
%   OUTPUT:
%       print_level
%           The desired print level
%   
%       printer
%           A printerContract, printerExpand, or printerInit printer.
%
%   See also printerContract, printerExpand, and printerInit.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   printerConfig.m introduced in ROSTAPACK Version 1.0
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
        
    print_level = opts.print_level;
    printer     = [];
    
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
    ascii       = isfield(opts,'print_ascii') && opts.print_ascii;
    
    switch lower(type)
        case 'contract'
            getPrinterFn    = @printerContract;
            col_specs       = getColumnSpecsContract(opts);
        case 'expand'
            getPrinterFn    = @printerExpand;
            col_specs       = getColumnSpecsExpand(opts);
        case 'init'
            getPrinterFn    = @printerInit;
            col_specs       = getColumnSpecsInit(opts);
        otherwise
            error('getPrinter: type ''%s'' is unrecognized',type); 
    end
                                         
    cols        = printerColumnFormatters(col_specs);
    printer     = getPrinterFn(ascii,cols);
end

function c = getColumnSpecsContract(opts)
    c = struct('maxit',opts.maxit,'step_type',true);
end

function c = getColumnSpecsExpand(opts)
    maxit       = opts.maxit;
    ls_maxit    = opts.line_search_opts.maxit;
    c = struct('maxit',maxit,'step_type',true,'ls_maxit',ls_maxit);
end

function c = getColumnSpecsInit(opts)
    if isfield(opts,'maxit')
        maxit       = opts.maxit;
    else
        maxit       = 1000;
    end
    kth_eigenvalue  = opts.initial_perturbation.kth_eigenvalue;
    c = struct('maxit',maxit,'kth_eigenvalue',kth_eigenvalue);
end
