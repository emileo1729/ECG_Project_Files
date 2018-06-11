%{
db = {};

for k = 1:10

  str = ["db", num2str(k)];
  load (str);
  db = [db; eval(str)];

endfor

db
%}

%{
coif = {};

for k = 1:5

  str = ["coif", num2str(k)];
  load (str);
  coif = [coif; eval(str)];

endfor

coif

save coif
%}

%{
sym = {[]};

for k = 2:8

  str = ["sym", num2str(k)];
  load (str);
  sym = [sym; eval(str)];

endfor

sym

save sym
%}

% [LO_D, HI_D, LO_R, HI_R] = wfilters ('wname')

[LO_D, HI_D, LO_R, HI_R] = arrayfun (@(x) wfilters (['db', num2str(x)]), 1:45, 'uniformoutput', false)

save db.mat LO_D HI_D LO_R HI_R
