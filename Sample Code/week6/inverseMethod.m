x = -10:0.01:10;
y1 = normpdf(x, 0, 1);
y2 = normcdf(x, 0, 1);

y1h = (normpdf(x, 0, 1) + normpdf(x, -2, 2) + normpdf(x, 2, 0.75))/3;
y2h = (normcdf(x, 0, 1) + normcdf(x, -2, 2) + normcdf(x, 2, 0.75))/3;

u1 = 0.82; %rand;

i1 = find(y2 <=u1);
x1 = x(i1(end));
i2 = find(y2h<=u1);
x2 = x(i2(end));

subplot(2,1,1);
plot(x,y1);
hold on;
plot(x,y1h);
hold off;
grid on;
ylabel('pdf');
xlabel('x');

subplot(2,1,2);
plot(x,y2);
hold on;
plot(x,y2h);
ylabel('cdf');
xlabel('x');


xX1 = 0:0.01:x1;
xX2 = 0:0.01:x2;

plot(xX1, u1*ones(size(xX1)), 'k:');
plot(xX2, u1*ones(size(xX2)), 'k:');

uU1 = 0:0.01:u1;

plot(x1*ones(size(uU1)),uU1, 'k:');
plot(x2*ones(size(uU1)),uU1, 'k:');


plot(x1,0, 'ko');
plot(x2,0, 'ko');

%plot(x1,u1, 'ko');
%plot(x2,u1, 'ko');

text(x1,-0.2, 'x_1');
text(x2,-0.2, 'x_2');
text(-0.2,u1, 'u_1');
%annotation('arrow', [.3 .5], [.6 .5]); 

hold off;
grid on;
saveas(gcf,'inveseMethod.png')
