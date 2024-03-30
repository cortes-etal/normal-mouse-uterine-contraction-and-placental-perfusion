function [bestFit,c,lags,r] = findhistfit(regions,animal)


figure;

subplot(1,2,1);

histogram(regions(regions > 0));
title('Whole Placenta Perfusion Histogrm');
ylabel('Frequency');
xlabel('Perfusion mL/min/100mL');

hold on;
xline(mean(regions(regions > 0)),'r','LineWidth',2);
hold off;

subplot(1,2,2);

h2= histogram(log10(regions(regions >0)));
tdata = h2.Values;

title('Whole Placenta Perfusion Histogram - Log Transform');
ylabel('Frequency');
xlabel('Log10(Perfusion mL/min/100mL)');

hold on;
xline(log10(mean(regions(regions > 0))),'r','LineWidth',2);
hold off;

sgtitle(animal);

data =tdata;%( (tdata - min(tdata)) / (range(tdata))) .* (2-0) + 0;
data = histcounts(log10(regions(regions > 0)), 'Normalization','pdf');
figure;
bar(data);
%
% pause;

% bestfit = 0;
% while bestfit <= .90
  for ii = 1:50  
    [f,xi] = ksdensity(log10(regions(regions > 0)));  % x is your data
    % Here I generate a function from two Gaussians and output
    % the rms of the estimation error from the values obtained from ksdensity
    fun = @(xx,t,y)peak2rms(y-(xx(5)*1./sqrt(xx(1)^2*2*pi).*exp(-(t-xx(2)).^2/(2*xx(1)^2))+...
        xx(6)*1./sqrt(xx(3)^2*2*pi).*exp(-(t-xx(4)).^2/(2*xx(3)^2))   )  );
    
    % Get the parameters with the minimum error. To improve convergence,choose reasonable initial values
 params = max(log10(regions(regions > 0))) .* rand(1,6);
    options = optimset('MaxFunEvals',50000,'Display','off','MaxIter',15000,'Tolfun',1E-20,'TolX',1E-20);
    [x,fval] = fminsearch(@(r)fun(r,xi,f),params,options);%[1.3 0.3 0.25 1 1.9 0.7]);
    
    % Make sure sigmas are positive
    x([1,3]) = abs(x([1,3]));
    
    disp(x)
    % Generate the Parametric functions
    pd1 = makedist('Normal','mu',x(2),'sigma',x(1));
    pd2 = makedist('Normal','mu',x(4),'sigma',x(3));
    % Get the probability values
    y1 = pdf(pd1,xi)*x(5); % x(5) is the participation factor from pdf1
    y2 = pdf(pd2,xi)*x(6); % x(6) is the participation factor from pdf2
    

   
    figure(111);
    plot(xi,f);
    hold on;
    plot(xi,y1);
    plot(xi,y2);
    plot(xi,y2+y1);
    legend({'ksdensity',['\mu : ',num2str(x(2)),'. \sigma :',num2str(x(1))],...
        ['\mu : ',num2str(x(4)),'. \sigma :',num2str(x(3))],'pdf1+pdf2'})
    title(['iter: ', num2str(ii)]);
    pause(0.25);
    hold off;
  end
    
    [c,lags] =xcorr(xi, y1+y2,'normalized');
    r = c(lags ==0);
    
%     pause;
    
    bestFit = r;
    
% end

end