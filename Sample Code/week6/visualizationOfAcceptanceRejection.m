clear;
close all;

xX = 0:0.01:5;
f1 = exp(-xX.*xX/2)/sqrt(2*pi);
f = sqrt(2/pi)*exp(-xX.*xX/2);

g2 = exp(-xX);

factor = 3;
c = factor*sqrt(2/pi)*exp(0.5);
g1 = c*g2;

figure(1);
h1 = plot(xX,f, 'b', 'linewidth', 2);
hold on;
h2 = plot(xX,g1, 'r', 'linewidth', 2);
plot(xX,g2, 'c', 'linewidth', 2);

nTrials = 10000000;
xArray = zeros(nTrials,1);

nSamples = 50000;
count = 0;
i = 0;
flagPlot = 'on';

indicators = zeros(nTrials,1);

while 1
    
    i = i+1;
    
    u = rand;
    x = exprnd(1);
    
    fX = sqrt(2/pi)*exp(-x*x/2);
    gX = exp(-x);
    ratio = fX/(c*gX);
    
    if strcmp(flagPlot,'on')
        h3 = plot([x, x], [0, fX], 'g-');
        h4 = plot([x, x], [fX, c*gX], 'r-');
    end
    
    if u <= ratio
        % accept
        count = count+1;
        indicators(i) = 1;
        if strcmp(flagPlot,'on')
            h5 = plot(x, u*c*gX, 'go'); 
        end
        % which side?
        u2 = rand;
        if u2 >0.5
            x = -x;
        end
    else  
        % reject
        indicators(i) = 0;
        if strcmp(flagPlot,'on')
            h6 = plot(x,  u*c*gX, 'ro');
        end
    end
    xArray(i) = x;
    %pause;
    
    if count >=5
        flagPlot = 'off';
    end
    
    if count >= nSamples
        break;
    end
end

legend([h1 h2 h3], {'f(x)', 'cg(x)', 'g(x)'});
saveas(gcf,'a-c3.png')


indicators = indicators(1:i);
xArray = xArray(1:i);
z = -5:0.01:5;
fz = normpdf(z,0,1);

figure(2);

subplot(1,2,1);
histogram(xArray,'Normalization','pdf');
hold on;
plot(z, fz, 'r','linewidth', 2);
hold off;

subplot(1,2,2);
histogram(xArray(indicators==1),'Normalization','pdf');
hold on;
plot(z, fz, 'r','linewidth', 2);
hold off;
saveas(gcf,'a-c4.png')

