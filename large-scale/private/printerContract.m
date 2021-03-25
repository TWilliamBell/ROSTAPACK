function printer = printerContract(ascii,cols)
%   printerContract:
%       Object for handling printing out the info for each iteration of
%       epsilonContract.
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
%           .step_type
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
%       .init(epsilon,f)            
%           Prints the initial value of epsilon and f, the value of the
%           root function, which is either the spectral abscissa or
%           spectral radius minus 1.
%       .stepEF(iter,epsilon,epsil_diff,f,f_diff,type)
%           Print the info for the current step where epsilon is updated. 
%           Requires arguments:
%               iter            The current iteration number
%               epsilon         The new value of epsilon 
%               epsilon_diff    The difference between epsilon and its old
%                               value
%               f               The new value of root function 
%               f_diff          The difference between f and its old value
%               step_type       At most 3 char string to indicate the type
%                               of step taken
%       .stepEFReject(iter,epsil,epsil_diff,f,f_diff,type)
%           Same as .stepEF but to be used when the step is not accepted, 
%           so that the values can be printed inside parentheses, to 
%           indicate that this candidate step was not accepted.  
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   printerContract.m introduced in ROSTAPACK Version 1.0
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
    type_c          = cols.step_type;
            
    labels          = { iter_c.label                                ...
                        epsilon_c.label     diff_c.label            ...
                        f_c.label           diff_c.label            ...
                        type_c.label                                };  
            
    widths          = [ iter_c.width                                ...
                        epsilon_c.width diff_c.width                ...
                        f_c.width       diff_c.width                ...
                        type_c.width                                ];
            
    span_labels     = { {'Epsilon',2,3} {'Root Fn',4,5} };
            
    table_printer   = tablePrinter(ascii,[],labels,widths,1,span_labels);        
         
    printer         = struct('msg',             @table_printer.msg,     ...
                             'close',           @table_printer.close,   ...
                             'init',            @init,                  ...
                             'stepEF',          @stepEFAccepted,        ...
                             'stepEFRejected',  @stepEFRejected         );
  
    function init(epsilon,f,varargin)
        table_printer.header();
        table_printer.row(  iter_c.format_fn(0),            ...
                            epsilon_c.format_fn(epsilon),   ...
                            diff_c.na_str,                  ...
                            f_c.format_fn(f),               ...
                            diff_c.na_str,                  ...
                            type_c.na_str                   );
    end

    function stepEF(  format_epsilon_fn,  format_f_fn,        ...
                    iter,                                   ...
                    epsilon,            epsilon_diff,       ...
                    f,                  f_diff,             ...
                    type,               varargin            )
        if mod(iter,20) == 0
            table_printer.header();
        end
        table_printer.row(  iter_c.format_fn(iter),         ...
                            format_epsilon_fn(epsilon),     ...
                            diff_c.format_fn(epsilon_diff), ...
                            format_f_fn(f),                 ...
                            diff_c.format_fn(f_diff),       ...
                            type_c.format_fn(type)          );
    end

    function stepEFAccepted(varargin)
        stepEF( @epsilon_c.format_fn,           ...
                @f_c.format_fn,                 ...
                varargin{:}                     );
    end

    function stepEFRejected(varargin)
        stepEF( @epsilon_c.format_rejected_fn,  ...
                @f_c.format_rejected_fn,        ...
                varargin{:}                     );
    end
end
