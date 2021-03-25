function [print_level,printers] = printersConfig(opts)
%   printersConfig:
%       This returns a struct of configured printer objects, along with 
%       the desired print level, for use with getStabRadBound and its two
%       main subroutines findUpperBound and hybridExpansionContraction.  It
%       is better to configure all necessary printers at once for two
%       reasons:
%
%       1) All the subroutines can reuse the printers, which is more
%          efficient than having the subroutines instantiate and configure
%          new printers every time they are called.
%
%       2) The printers can be configured so that all the column borders
%           across printer types are nicely aligned, which makes the
%           console output much easier to read.
%      
%   INPUT:
%       opts                [required: struct of parameters]
%           A required struct of settable parameters necessary to run this
%           routine, which determine which printers will be created and
%           their exact configurations:
%
%           .maxit 
%               Required when opts.print_level > 0.  This is the maximum
%               number of iterations allowed for the
%               hybridExpansionContraction routine, which uses a
%               printerMain printer.
%
%           .initial_perturbation.kth_eigenvalue
%               If this is present, a printerInit printer will be created.
%                
%           .upperbound_opts.maxit
%           .upperbound_opts.UV_steps_per_iter
%           .upperbound_opts.epsilon_steps_per_iter
%               Required when opts.print_level > 0 and a printerMain
%               printer needs to be configured for the findUpperBound
%               routine.
%
%           .upperbound_opts.epsilon_line_search_opts.maxit
%               Required when opts.print_level > 1 and printers are needed
%               for findUpperBound's subroutines.    
%
%           .expansion_opts.maxit   
%               Required when opts.print_level > 0.  The printerExpand
%               printer is needed by subroutines of both the findUpperBound
%               and hybridExpansionContraction routines.
%
%           .expansion_opts.line_search_opts.maxit;
%               Required when opts.print_level > 1 and printers are needed
%               for hybridExpansionContraction's subroutines.   
%
%           .contraction_opts.maxit
%               Required when opts.print_level > 0 and a printerMain
%               printer needs to be configured for the
%               hybridExpansionContraction routine.
%
%           .print_level              
%               The desired print level.
% 
%           .print_ascii            
%               Optional.  Set to true to fallback to standard ASCII
%               character set for printing table borders.
%
%           .printers
%               If all the printer[Type] objects are already configured,
%               they can be provided in this struct, to prevent
%               instantiating new ones.  The fields are:
%                   .main
%                   .init
%                   .expand
%                   .contract
%           Warning: it is assumed that the preconfigured printers are
%           compatible with the parameters stored in opts and that all
%           printers required are present in this struct.
%
%
%   OUTPUT:
%       print_level
%           The desired print level
%   
%       printers 
%           A struct of the configured printers:
%               .main       printerMain object
%               .init       printerInit object
%               .expand     printerExpand object
%               .contract   printerContract
%
%           Note that .init, .expand, and .contract may be [] if they are
%           not needed, as determined by opts.print_level and which
%           algorithm parameter fields are present in opts.
%
%   See also printerContract, printerExpand, printerInit, and printerMain.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   printersConfig.m introduced in ROSTAPACK Version 1.0
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

    print_level     = opts.print_level;
    printers        = [];
    
    % nothing to do if printing is disabled
    if print_level < 1
        return
    end
    
    % Check if the printers are already configured and if so, return 
    if isfield(opts,'printers')  
        printers = opts.printers;
        return
    end
    
    % Otherwise, we need to configure all the printers that are needed,
    % based on what options are present
    
    ascii           = isfield(opts,'print_ascii') && opts.print_ascii;
    
    [ub_maxit,ub_sub_maxit,epsilon_ls_maxit] = getUpperBoundMaxits(opts);
    contract_maxit  = getContractionMaxit(opts);
    expand_maxit    = opts.expansion_opts.maxit;
    
    maxit           = max(opts.maxit,ub_maxit);
    phase_maxit     = max([ub_sub_maxit contract_maxit expand_maxit]);      
    
    if isfield(opts,'initial_perturbation')
        kth_evals   = opts.initial_perturbation.kth_eigenvalue;
        init_opts   = {'kth_eigenvalue', kth_evals};  
    else
        init_opts   = {};
    end
       
    if opts.print_level > 1
        maxit       = max(maxit,phase_maxit);
        UV_ls_maxit = opts.expansion_opts.line_search_opts.maxit;
        ls_maxit    = max(epsilon_ls_maxit,UV_ls_maxit);
        sub_opts    = {'step_type', true, 'ls_maxit', ls_maxit};
    else
        sub_opts    = {};
    end
    
    col_opts        = getColumnSpecs(maxit,phase_maxit,init_opts,sub_opts);
    cols            = printerColumnFormatters(col_opts);  
    
    printers        = struct(   'main',     printerMain(ascii,cols),    ...
                                'init',     [],                         ...
                                'expand',   [],                         ...
                                'contract', []                          );
                            
    if ~isempty(init_opts)
        printers.init           = printerInit(ascii,cols);
    end
    
    if ~isempty(sub_opts)
        printers.expand         = printerExpand(ascii,cols);
        if contract_maxit > 0
            printers.contract   = printerContract(ascii,cols);
        end
    end
end

function [maxit,sub_maxit,epsilon_maxit] = getUpperBoundMaxits(opts)
    if isfield(opts,'upperbound_opts')
        ub_opts         = opts.upperbound_opts;
        maxit           = ub_opts.maxit;
        UV_maxit        = ub_opts.UV_steps_per_iter;
        epsilon_maxit   = ub_opts.epsilon_steps_per_iter;
        sub_maxit       = max(UV_maxit,epsilon_maxit);
        epsilon_maxit   = ub_opts.epsilon_line_search_opts.maxit;
    else
       	[maxit,sub_maxit,epsilon_maxit] = deal(0,0,0);
    end
end

function maxit = getContractionMaxit(opts)
    if isfield(opts,'contraction_opts')
        maxit = opts.contraction_opts.maxit;
    else
       	maxit = 0;
    end
end

function c = getColumnSpecs(maxit,phase_maxit,init_opts,sub_opts)
    c = struct( 'maxit',        maxit,          ...
                'phase_maxit',  phase_maxit,    ...
            	init_opts{:},                   ...
                sub_opts{:}                     );
end
