free1 = load('freeBoundary2.dat');
free2 = load('freeBoundary2.dat');

plot(free1(:,2), free1(:,1), 'k-.');
hold on;
plot(free2(:,2), free2(:,1), 'k-');
hold off;
xlabel('spot', 'fontsize', 12);
ylabel('time to maturity', 'fontsize', 12);
title('Exercise boundary VG vs. GBM');
axis tight;
legend('GBM', 'VG', 'Location', 'Best');
eval(['print -deps2 ', 'freeboundaryVG_GBM']);