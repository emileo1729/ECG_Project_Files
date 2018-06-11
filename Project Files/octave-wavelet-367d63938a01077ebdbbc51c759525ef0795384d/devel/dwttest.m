t = 0:10:360;
x = sind (t);

[u, v] = dwt (x, 'db1')

[u, v] = dwt (x, 'db1', 'mode', 'n')
