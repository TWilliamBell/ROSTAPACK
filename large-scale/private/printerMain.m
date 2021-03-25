function printer = printerMain(ascii,cols)
%   printerMain:
%       Object for handling printing out the info for each iteration of
%       either findUpperBound or hybridExpansionContraction.
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
%           .phase
%           See printerColumnFormatters for how each of these are defined.
%   
%   OUTPUT:
%       An "object", a struct containing the following functions for
%       printing tasks:
%       .msg(a_msg_str_or_cell_array)               
%           Prints a message inside the table; the full width is available
%           and cell arrays can be used to print multiple lines.
%       .header()
%           Prints all the labels for each column of data.
%       .close()     
%           Closes the table.
%       .init(epsilon,f)            
%           Prints the initial value of epsilon and f, the value of the
%           root function, which is either the spectral abscissa or
%           spectral radius minus 1.
%       .stepF(iter,accepted,f,f_diff,type,siters,tc)
%           Print the info for the current phase where epsilon is held
%           constant but perturbation U*V' is updated by the uvExpand
%           subroutine.  Requires arguments:
%               iter            The current iteration number
%               accepted        Logical indicating whether or not progress 
%                               was made by uvExpand
%               f               The possibly new value of root function 
%               f_diff          The difference between f and its old value
%               type            At most 3 char string to indicate the type
%                               of the phase (e.g. 'EXP' for expansion)
%               siters          Number of iterations incurred by uvExpand 
%               tc              Integer, uvExpand's halt status. 
%       .stepEF(iter,accepted,epsilon,epsilon_diff,f,f_diff,type,siters,tc) 
%           Print the info for the current phase where epsilon is updated
%           by either epsilonContract or epsilonExpand, with the
%           perturbation U*V' fixed.  Requires arguments:
%               iter            The current iteration number
%               accepted        Logical indicating whether or not progress 
%                               was made by epsilonContract or
%                               epsilonExpand
%               epsilon         The possibly new value of epsilon 
%               epsilon_diff    The difference between epsilon and its old
%                               value
%               f               The possibly new value of root function 
%               f_diff          The difference between f and its old value
%               type            At most 3 char string to indicate the type
%                               of the phase (e.g. 'EXP' for expansion or
%                               'CON' for contraction)
%               siters          Number of iterations incurred by 
%                               epsilonContract or epsilonExpand
%               tc              Integer, the halt status of epsilonContract
%                               or epsilonExpand.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   printerMain.m introduced in ROSTAPACK Version 1.0
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
    type_c          = cols.phase.type;
    siters_c        = cols.phase.subiters;
    term_c          = cols.phase.term;
              
    labels          = { iter_c.label                                    ...
                        epsilon_c.label diff_c.label                    ...
                        f_c.label       diff_c.label                    ...
                        type_c.label    siters_c.label  term_c.label    };                              
            
    widths          = [ iter_c.width                                    ...
                        epsilon_c.width diff_c.width                    ...
                        f_c.width       diff_c.width                    ...
                        type_c.width    siters_c.width  term_c.width    ]; 
            
    span_labels     = { {'Epsilon',2,3} ...
                        {'Root Fn',4,5} ...
                        {'Phase',6,8}   };
                    
    format_tc_fn    = @(x) term_c.format_fn(num2str(x));
            
    table_printer   = tablePrinter(ascii,[],labels,widths,1,span_labels);        
         
    printer         = struct('msg',             @table_printer.msg,     ...
                             'header',          @table_printer.header,  ...
                             'close',           @table_printer.close,   ...
                             'init',            @init,                  ...
                             'stepF',           @stepF,                 ...
                             'stepEF',          @stepEF                 );
  
    function init(epsilon,f,varargin)
        table_printer.header();
        table_printer.row(  iter_c.format_fn(0),            ...
                            epsilon_c.format_fn(epsilon),   ...
                            diff_c.na_str,                  ...
                            f_c.format_fn(f),               ...
                            diff_c.na_str,                  ...
                            type_c.na_str,                  ...
                            siters_c.na_str,                ...
                            term_c.na_str                   );
    end

    function stepF(iter,accepted,f,f_diff,type,siters,tc,varargin)
        if accepted
            format_f_fn = f_c.format_fn;
        else
            format_f_fn = f_c.format_rejected_fn;
        end
        table_printer.row(  iter_c.format_fn(iter),         ...
                            epsilon_c.blank_str,            ...
                            diff_c.blank_str,               ...
                            format_f_fn(f),                 ...
                            diff_c.format_fn(f_diff),       ...
                            type_c.format_fn(type),         ...
                            siters_c.format_fn(siters),     ...
                            format_tc_fn(tc)                );
    end

    function stepEF(iter,               accepted,           ...     
                    epsilon,            epsilon_diff,       ...
                    f,                  f_diff,             ...
                    type,               siters,         tc  )
        if accepted
            format_epsilon_fn   = epsilon_c.format_fn;
            format_f_fn         = f_c.format_fn;
        else
            format_epsilon_fn   = epsilon_c.format_rejected_fn;
            format_f_fn         = f_c.format_rejected_fn;
        end
        table_printer.row(  iter_c.format_fn(iter),         ...
                            format_epsilon_fn(epsilon),     ...
                            diff_c.format_fn(epsilon_diff), ...
                            format_f_fn(f),                 ...
                            diff_c.format_fn(f_diff),       ...
                            type_c.format_fn(type),         ...
                            siters_c.format_fn(siters),     ...
                            format_tc_fn(tc)                );
    end
end
