function printer = printerExpand(ascii,cols)
%   printerExpand:
%       Object for handling printing out the info for each iteration of
%       either uvExpand or epsilonExpand.
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
%           .line_search
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
%       .stepF(iter,f,f_diff,step_type,ls_iters,t)
%           Print the info for the current step where epsilon is held
%           constant but perturbation U*V' is updated.  Requires arguments:
%               iter            The current iteration number
%               f               The new value of root function 
%               f_diff          The difference between f and its old value
%               step_type       At most 3 char string to indicate the type
%                               of step taken
%               ls_iters        Number of evaluations incurred by the line 
%                               search        
%               t               The value of t accepted by the line search
%       .stepFRejected(iter,f,f_diff,step_type,ls_iters,t)
%           Same as .stepF but to be used when f_diff <= 0, so that the
%           values can be printed inside parentheses, to indicate that this
%           candidate step was not accepted.
%       .stepEF(iter,epsilon,epsil_diff,f,f_diff,type,ls_iters,t)
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
%               ls_iters        Number of evaluations incurred by the line 
%                               search        
%               t               The value of t accepted by the line search
%       .stepEFReject(iter,epsil,epsil_diff,f,f_diff,type,ls_iters,t)
%           Same as .stepEF but to be used when f_diff <= 0, so that the
%           values can be printed inside parentheses, to indicate that this
%           candidate step was not accepted.
%       .lineSearchF(ft,ft_diff,iter,t)
%           Print the info for a rejected line search evaluation where
%           epsilon is held constant but perturbation U*V' is updated.
%           Requires arguments:
%               ft              The value of root function at parameter t
%                               in [0,1]
%               ft_diff         f(t) - f(0) 
%               iter            Line search iteration number 
%               t               The current value of t attempted.
%       .lineSearchEF(epsilont,epsilont_diff,ft,ft_diff,iter,t)
%           Print the info for a rejected line search evaluation where 
%           epsilon.  Requires arguments:
%               epsilont        The value of epsilon for the line search
%                               parameter t in [0,1]
%               epsilont_diff   epsilon(t) - epsilon(0) 
%               ft              The value of root function at parameter t
%                               in [0,1]
%               ft_diff         f(t) - f(0) 
%               iter            Line search iteration number 
%               t               The current value of t attempted.
%       .interp(ft,ft_diff,t)
%           Print the info for an accepted interpolation where epsilon is
%           held constant.  Requires arguments:
%               ft              The value of root function at interpolation
%                               parameter t in [0,1]
%               ft_diff         f(t) - f(0) 
%               t               The interpolation value t.
%       .interpRejected(ft,ft_diff,t)
%           Same as .interp, but to be used when ft_diff <= 0, so that the
%           values can be printed inside parentheses, to indicate that this
%           interpolation was not accepted.
%       .extrap(f,f_diff)
%           Print the info for an accepted extrapolation where epsilon is
%           held constant.  Requires arguments:
%               f               The value of root function at the 
%                               extrapolation
%               f_diff          f minus the old value of f.
%       .extrapRejected(f,f_diff)
%           Same as .extrap, but to be used when f_diff <= 0, so that the
%           values can be printed inside parentheses, to indicate that this
%           extrapolation was not accepted.
%
%
%   For comments/bug reports, please visit the ROSTAPACK GitLab webpage:
%   https://gitlab.com/timmitchell/ROSTAPACK
%
%   printerExpand.m introduced in ROSTAPACK Version 1.0
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
    ls_iters_c      = cols.line_search.iters;
    ls_t_c          = cols.line_search.t;
            
    labels          = { iter_c.label                                ...
                        epsilon_c.label     diff_c.label            ...
                        f_c.label           diff_c.label            ...
                        type_c.label                                ...
                        ls_iters_c.label    ls_t_c.label            };  
            
    widths          = [ iter_c.width                                ...
                        epsilon_c.width diff_c.width                ...
                        f_c.width       diff_c.width                ...
                        type_c.width                                ...
                        ls_iters_c.width    ls_t_c.width            ];
            
    span_labels     = { {'Epsilon',2,3}         ...
                        {'Root Fn',4,5}         ...
                        {'Line Search',7,8}     }; 
            
    table_printer   = tablePrinter(ascii,[],labels,widths,1,span_labels);        
         
    interp_str      = type_c.format_fn('I');
    extrap_str      = type_c.format_fn('E');
    ls_str          = type_c.format_fn('LS');
    
    printer         = struct('msg',             @table_printer.msg,     ...
                             'close',           @table_printer.close,   ...
                             'init',            @init,                  ...
                             'stepF',           @stepFAccepted,         ...
                             'stepFRejected',   @stepFRejected,         ...
                             'stepEF',          @stepEFAccepted,        ...
                             'stepEFRejected',  @stepEFRejected,        ...
                             'lineSearchF',     @lineSearchF,           ...
                             'lineSearchEF',    @lineSearchEF,          ...
                             'interp',          @interpAccepted,        ...
                             'interpRejected',  @interpRejected,        ...
                             'extrap',          @extrapAccepted,        ...
                             'extrapRejected',  @extrapRejected         );
  
    function init(epsilon,f,varargin)
        table_printer.header();
        table_printer.row(  iter_c.format_fn(0),            ...
                            epsilon_c.format_fn(epsilon),   ...
                            diff_c.na_str,                  ...
                            f_c.format_fn(f),               ...
                            diff_c.na_str,                  ...
                            type_c.na_str,                  ...
                            ls_iters_c.na_str,              ...
                            ls_t_c.na_str                   );
    end

    function stepF(format_f_fn,iter,f,f_diff,type,ls_iters,ls_t,varargin)
        if mod(iter,20) == 0
            table_printer.header();
        end
        table_printer.row(  iter_c.format_fn(iter),         ...
                            epsilon_c.blank_str,            ...
                            diff_c.blank_str,               ...
                            format_f_fn(f),                 ...
                            diff_c.format_fn(f_diff),       ...
                            type_c.format_fn(type),         ...
                            ls_iters_c.format_fn(ls_iters), ...
                            ls_t_c.format_fn(ls_t)          );
    end

    function stepFAccepted(varargin)
        stepF(@f_c.format_fn,varargin{:});
    end

    function stepFRejected(varargin)
        stepF(@f_c.format_rejected_fn,varargin{:});
    end

    function stepEF(    format_epsilon_fn,  format_f_fn,                ...
                        iter,                                           ...
                        epsilon,            epsilon_diff,               ...
                        f,                  f_diff,                     ...
                        type,                                           ...
                        ls_iters,           ls_t,           varargin    )
                    
        if mod(iter,20) == 0
            table_printer.header();
        end
        table_printer.row(  iter_c.format_fn(iter),         ...
                            format_epsilon_fn(epsilon),     ...
                            diff_c.format_fn(epsilon_diff), ...
                            format_f_fn(f),                 ...
                            diff_c.format_fn(f_diff),       ...
                            type_c.format_fn(type),         ...
                            ls_iters_c.format_fn(ls_iters), ...
                            ls_t_c.format_fn(ls_t)          );
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

    function lineSearchF(f,f_diff,ls_iter,ls_t,varargin)
        table_printer.row(  iter_c.blank_str,               ...
                            epsilon_c.blank_str,            ...
                            diff_c.blank_str,               ...
                            f_c.format_rejected_fn(f),      ...
                            diff_c.format_fn(f_diff),       ...
                            ls_str,                         ...
                            ls_iters_c.format_fn(ls_iter),  ...
                            ls_t_c.format_fn(ls_t)          );
    end

    function lineSearchEF(epsilon,epsilon_diff,f,f_diff,ls_iter,ls_t)
        table_printer.row(  iter_c.blank_str,                       ...
                            epsilon_c.format_rejected_fn(epsilon),  ...
                            diff_c.format_fn(epsilon_diff),         ...
                            f_c.format_rejected_fn(f),              ...
                            diff_c.format_fn(f_diff),               ...
                            ls_str,                                 ...
                            ls_iters_c.format_fn(ls_iter),          ...
                            ls_t_c.format_fn(ls_t)                  );
    end

    function interp(format_f_fn,f,f_diff,ls_t,varargin)
        table_printer.row(  iter_c.blank_str,           ...
                            epsilon_c.blank_str,        ...
                            diff_c.blank_str,           ...
                            format_f_fn(f),             ...
                            diff_c.format_fn(f_diff),   ...
                            interp_str,                 ...
                            ls_iters_c.na_str,          ...
                            ls_t_c.format_fn(ls_t)      );
    end
    
    function interpAccepted(varargin)
        interp(@f_c.format_fn,varargin{:});
    end

    function interpRejected(varargin)
        interp(@f_c.format_rejected_fn,varargin{:});
    end

    function extrap(format_f_fn,f,f_diff,varargin)
        table_printer.row(  iter_c.blank_str,           ...
                            epsilon_c.blank_str,        ...
                            diff_c.blank_str,           ...
                            format_f_fn(f),             ...
                            diff_c.format_fn(f_diff),   ...
                            extrap_str,                 ...
                            ls_iters_c.na_str,          ...
                            ls_t_c.na_str               );
    end

    function extrapAccepted(varargin)
        extrap(@f_c.format_fn,varargin{:});
    end

    function extrapRejected(varargin)
        extrap(@f_c.format_rejected_fn,varargin{:});
    end
end
