%Plotting function
function plotReq(xAx, yAx, xlab, ylab, titl, lim, name)
    plot(xAx, yAx, 'DisplayName', name, 'LineWidth', 1.5);
    xlabel(xlab);
    ylabel(ylab);
    title(titl);
    ylim(lim);
end