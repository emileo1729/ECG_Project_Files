load db.mat

H = LO_R;

G = cellfun (@fliplr, H, "uniformoutput", false);

for k = 1 : numel (G)
  G{k}(2:2:end) = -G{k}(2:2:end);
endfor


Hp = cellfun (@fliplr, H, "uniformoutput", false);

Gp = cellfun (@fliplr, G, "uniformoutput", false);


a = isequal (H, LO_R)
b = isequal (G, HI_R)

c = isequal (Hp, LO_D)
d = isequal (Gp, HI_D)

n = cellfun (@norm, H)