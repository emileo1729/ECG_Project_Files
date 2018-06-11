function y = wconv1 (x, f, shape = "full")

  if (nargin < 2 || nargin > 3)
    print_usage ();
  endif

  y = conv2 (x(:).', f(:).', shape);

  if (rows (x) > 1)
    y = y.';
  endif

endfunction
