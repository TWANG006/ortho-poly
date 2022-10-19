function show_surface(X, Y, Z, unit, unit_str, title_str)

X = X * 1e3;
Y = Y * 1e3;
Z = Z * unit;

ratio = (max(X(:)) - min(X(:))) / (max(Y(:)) - min(Y(:)));

% display the map
surf(X, Y, Z, 'EdgeColor', 'none'); alpha(gca, 0.5); hold on;
p = pcolor(X, Y, Z); p.EdgeColor = 'none';
axis(gca, 'tight');
c = colorbar;
c.Label.String = unit_str;
title({title_str,...
      ['RMS = ' num2str(round(nanstd(Z(:),1), 2)) ' ' unit_str]},...
    'FontWeight', 'normal',...
    'FontSize', 12 ...
);
xlabel('x [mm]', 'FontSize', 12);
ylabel('y [mm]', 'FontSize', 12);
colormap(gca, 'jet');
shading(gca, 'flat');
pbaspect([ratio, 1, 1]);

end