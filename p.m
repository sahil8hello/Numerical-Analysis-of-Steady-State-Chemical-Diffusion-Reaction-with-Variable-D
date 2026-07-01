function p(x, a, b, c, d)

hold on; grid on;
title('Numerical solution Cubic Spline');

hSpline = gobjects(1,1);  % reserve graphics handle

for i = 1:length(x)-1

    Si = @(X) a(i) ...
        + b(i).*(X - x(i)) ...
        + c(i).*(X - x(i)).^2 ...
        + d(i).*(X - x(i)).^3;

    % Plot spline
    h = fplot(Si, [x(i), x(i+1)], 'y');

    if i == 1
        hSpline = h; % store 1st plotted yellow curve as representative
    end
end

legend(hSpline, 'Numerical Solution (Yellow)');
xlabel('L (cm)');
ylabel('Concentration');

end
