%% int_0^1 int_0^1 x^3 (1+y^2) dx dy = 1/3

n = 2000;

%X = lhsdesign(n,p)
%X = lhsdesign(...,'smooth','off')
%X = lhsdesign(...,'criterion',criterion)
%X = lhsdesign(...,'iterations',k)

p = haltonset(2, 'Skip', 1e3, 'Leap', 1e2);
p = scramble(p, 'RR2');
X0 = net(p, n);
figure(1);
scatter(X0(:,1), X0(:,2), 5, 'k');
axis square;
% int_0^1 int_0^1 x^3 (1+y^2) dx dy = 1/3
f1 = (X0(:,1).^3).*(1+X0(:,2).^2);
I1 = mean(f1);
% title(['low discrepancy sequence (Halton set) ', num2str(I1)]);
legend('low discrepancy sequence (Halton set)', 'Location','Best');
print -deps HaltonSet;
% saveas(gcf,'quasiMC1.png')

U0 = rand(n,2);
figure(2);
scatter(U0(:,1), U0(:,2), 5, 'k');
axis square;
f2 = (U0(:,1).^3).*(1+U0(:,2).^2);
I2 = mean(f2);
% title(['uniform random variables ', num2str(I2)]);
legend('uniform random variables', 'Location','Best');
print -deps uniformRand;
% saveas(gcf,'quasiMC2.png');

disp([1/3 I1 I2]);



