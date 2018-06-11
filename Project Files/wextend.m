## Copyright (C) 2013   Lukas F. Reichlin
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{y} =} wextend (@var{type}, @var{mode}, @var{x}, @var{len})
## @deftypefnx {Function File} {@var{y} =} wextend (@var{type}, @var{mode}, @var{x}, @var{len}, @var{loc})
## Extend a signal.
## 
## @strong{Inputs}
## @table @var
## @item type
## Type.
## @item mode
## Mode.
## @item x
## Signal as a vector (1-D) or matrix (2-D).
## @item len
## Length of the extension.
## @item loc
## Location of the extension.
## @end table
##
## @strong{Outputs}
## @table @var
## @item y
## Extended signal.
## @end table
## @end deftypefn

## Author: Lukas Reichlin <lukas.reichlin@gmail.com>
## Created: June 2013
## Version: 0.1

function y = wextend (type, mode, x, len, location)

  if (nargin < 4 || nargin > 5)
    print_usage ();
  endif

endfunction
