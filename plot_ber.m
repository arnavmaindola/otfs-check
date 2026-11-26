function h = plot_ber(x, B, nbits_total, ylo, styles)
% PLOT_BER  Semilog BER curves that keep zero-error points visible.
%   x           : abscissa (nx x 1)
%   B           : BER matrix (nx x nc)
%   nbits_total : bits simulated per point (floor = 1/nbits_total)
%   ylo         : lower y-limit; points below it are drawn at ylo
%   styles      : cell array of {color, 'linespec'} per column
% Zero-error points are drawn at the floor with a hollow marker.
% Returns the line handles (one per column) for the legend.
    floor_ber = 1 / nbits_total;
    nc = size(B, 2); h = gobjects(nc, 1);
    hold on;
    for k = 1:nc
        c = styles{k}{1}; ls = styles{k}{2};
        y = max(B(:, k), floor_ber);
        y = max(y, ylo);
        h(k) = semilogy(x, y, ls, 'Color', c, 'LineWidth', 1.5, 'MarkerFaceColor', c);
        z = B(:, k) == 0;
        if any(z)
            semilogy(x(z), y(z), ls(end), 'Color', c, 'LineStyle', 'none', ...
                'MarkerFaceColor', 'w', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    set(gca, 'YScale', 'log'); grid on; ylim([ylo 1]);
end
