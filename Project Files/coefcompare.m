%[LO_D, HI_D, LO_R, HI_R] = arrayfun (@(x) wfilters (['db', num2str(x)]), 1:45, 'uniformoutput', false)

load db.mat

[Hp, Gp, H, G] = arrayfun (@(x) wfilters (['db', num2str(x)]), 1:45, 'uniformoutput', false);

a = isequal (H, LO_R)
b = isequal (G, HI_R)

c = isequal (Hp, LO_D)
d = isequal (Gp, HI_D)

dH = cellfun (@(x, y) abs (x-y), H, LO_R, "uniformoutput", false)
dG = cellfun (@(x, y) abs (x-y), G, HI_R, "uniformoutput", false)

dHp = cellfun (@(x, y) abs (x-y), Hp, LO_D, "uniformoutput", false)
dGp = cellfun (@(x, y) abs (x-y), Gp, HI_D, "uniformoutput", false)


dHmax = cellfun (@(x, y) max (abs (x-y)), H, LO_R)
dGmax = cellfun (@(x, y) max (abs (x-y)), G, HI_R)

dHpmax = cellfun (@(x, y) max (abs (x-y)), Hp, LO_D)
dGpmax = cellfun (@(x, y) max (abs (x-y)), Gp, HI_D)
