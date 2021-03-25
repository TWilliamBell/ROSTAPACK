function printer = printerInit(ascii,cols)
%   printerInit:
%       Object for handling initializePerturbation's printing.
%       
%   INPUT:
%       ascii  
%           A logical indicating whether or not to restrict to the ASCII
%           character set.  If it is false, extended chars will be used to
%           produce nicer looking tables.
% 
%       cols
%           A struct of formatters created by printerColumnFormatters.  The
%           required formatters are:
%           .iter
%           .main
%           .init_pert
%           See printerColumnFormatters for how each of these are defined.
%   
%   OUTPUT:
%       An "object", a struct containing the following functions for
%       printing tasks:
%       .msg(a_msg_str_or_cell_array)               
%           Prints a message inside the table; the full width is available
%           and cell arrays can be used to print multiple lines.
%       .close()     
%           Closes the table.
%       .matrixA(f,kth_eigenvalue)            
%           Prints the value of the root function f for the kth eigenvalue
%           of A, which is either the spectral abscissa or spectral radius
%           minus 1.
%       .initPert(epsilon,f,f_diff,kth_eigenvalue)
%           Print the info for the first perturbation: epsilon*U*V'. 
%           Requires arguments:
%               epsilon         The initial value of epsilon 
%               f               The value of root function for epsilon*U*V'
%               f_diff          The difference between f and the value of 
%                               the user-selected eigenvalue of A 
%               kth_eigenvalue  Which eigenvalue was requested to be used 
%                               for this perturbation. 
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   printerInit.m introduced in ROSTAPACK Version 1.0
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
               
    iter_c          = cols.iter;
    epsilon_c       = cols.main.epsilon;
    f_c             = cols.main.f;
    diff_c          = cols.main.diff;
    type_c          = cols.init_pert.type;
    kth_eval_c      = cols.init_pert.kth_eval;
              
    labels          = { iter_c.label                                    ...
                        epsilon_c.label diff_c.label                    ...
                        f_c.label       diff_c.label                    ...
                        type_c.label    kth_eval_c.label                };                              
            
    widths          = [ iter_c.width                                    ...
                        epsilon_c.width diff_c.width                    ...
                        f_c.width       diff_c.width                    ...
                        type_c.width    kth_eval_c.width                ]; 
            
    span_labels     = { {'Epsilon',2,3} ...
                        {'Root Fn',4,5} ...
                        {'Matrix',6,7}  };
                    
            
    table_printer   = tablePrinter(ascii,[],labels,widths,1,span_labels);        
        
    printer         = struct('msg',             @table_printer.msg,     ...
                             'close',           @table_printer.close,   ...
                             'matrixA',         @matrixA,               ...
                             'initPert',        @initPert               );
                    
    print_header    = true;
    
    function matrixA(f,kth,varargin)
        if print_header
            table_printer.header();
            print_header = false;
        end
        table_printer.row(  iter_c.blank_str,               ...
                            epsilon_c.blank_str,            ...
                            diff_c.blank_str,               ...
                            f_c.format_fn(f),               ...
                            diff_c.blank_str,               ...
                            type_c.format_fn('A'),          ...
                            kth_eval_c.format_fn(kth)       );
    end

    function initPert(epsilon,f,f_diff,kth,varargin)
        if print_header
            table_printer.header();
            print_header = false;
        end
        if isempty(f_diff) || isnan(f_diff)
            f_diff_str = diff_c.blank_str;
        else
            f_diff_str = diff_c.format_fn(f_diff);
        end
        table_printer.row(  iter_c.blank_str,               ...
                            epsilon_c.format_fn(epsilon),   ...
                            diff_c.blank_str,               ...
                            f_c.format_fn(f),               ...
                            f_diff_str,                     ...
                            type_c.format_fn('M0'),         ...
                            kth_eval_c.format_fn(kth)       );
    end
end
