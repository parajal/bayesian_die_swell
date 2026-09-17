% Reads the points of a TFEM pcurve_t and fills the space in between.

for i=1:100
  [X,Y] = textread(['pcurvecoor' sprintf('%03d',i) '.txt']);
  fill(X,Y,'black')
  axis([0 1 0 1])
  pause(0.2);
end
