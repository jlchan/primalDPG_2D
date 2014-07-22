
function PlotMeshQuad2D()

  Globals2D;

  quads = EToV(:,[1 2 3 4 1])';
  ha = plot(VX(quads), VY(quads), 'LineWidth', 1.2);
  set(ha, 'color', 'black');
