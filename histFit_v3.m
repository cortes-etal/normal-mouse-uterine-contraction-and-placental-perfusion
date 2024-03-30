function [pdfDat, skew] = histFit_v3(regions,animal,fname,flag)

dat = log10(regions(regions>0)); % log transform of data
options = optimset('Display','off','MaxIter',1000,'TolFun',1E-8);

if flag == 0
    comp = 2;
    
    mdl= fitgmdist(dat,comp,'Options',options); % fit the mixed gaussian
    
    figure(112); % inti fig
    paramEsts= mdl;
    %mu of mixture of gaussians of 2 components % gaussian parameters
    MU=[paramEsts.mu(1);paramEsts.mu(2)];
    SIGMA=cat(3,[paramEsts.Sigma(1)],[paramEsts.Sigma(2)]);
    PPp=[paramEsts.PComponents(1),paramEsts.PComponents(2)];
    
    objA = gmdistribution(MU,SIGMA,PPp); % creating fit objects
    objA1 = gmdistribution(MU(1),SIGMA(:,:,1),PPp(1));
    objA2 = gmdistribution(MU(2),SIGMA(:,:,2),PPp(2));
    
    xi = linspace(min(dat(:)), max(dat(:)),100); % index var
    
    y1 = pdf(objA1,xi.').*PPp(1); % first pdf (low perf?)
    y2 = pdf(objA2,xi.').*PPp(2); % (high perf)?
    ycomb = pdf(objA,xi.'); % combined pdf (looks like histogram)
    
    skew = skewness(regions(regions>0));
    
    %     begin plotting
    plot(xi,ycomb,'r-','linewidth',2)
    hold on;
    plot(xi,y1,'g-','Linewidth',2);
    plot(xi,y2,'b-','Linewidth',2);
 %   histogram(log10(regions(regions > 0)),'Normalization','pdf','NumBins',numel(xi),'FaceAlpha',0);
    lgd = legend('Combined PDF',sprintf('Mu = %0.2f, Sigma = %0.3f',MU(1), SIGMA(:,:,1)),sprintf('Mu = %0.2f, Sigma = %0.3f',MU(2), SIGMA(:,:,2)),'histogram');
    lgd.FontSize = 16;
    title(animal);
    xlim([0 2.5]);
    ylim([0 6]);
    xlabel('Log_{10}( Perfusion [mL / min / 100 mL]');
    ylabel('PDF');
    hold off;
    saveas(gcf,fname,'png'); % saving historam figure
    
%     pause;
    aoc1=trapz(xi,pdf(objA1,xi.').*PPp(1)); % getting area under curves
    aoc2=trapz(xi,pdf(objA2,xi.').*PPp(2));
    
    hh1 = (min(y1) + max(y1)) / 2;
    hh2 = (min(y2) + max(y2)) / 2;
    
    xin1 = find(y1 >= hh1,1,'first');
    xin2 = find(y1 >= hh1,1,'last');
    
    x1 = xi(xin1);
    x2 = xi(xin2);
    fwhm = x2-x1;
    
    xin12 = find(y2 >= hh2,1,'first');
    xin22 = find(y2 >= hh2,1,'last');
    
    x12 = xi(xin12);
    x22 = xi(xin22);
    fwhm2 = x22-x12;
    
    
    
    %create output var
    pdfDat(1).y=y1;
    pdfDat(2).y=y2;
    
    pdfDat(1).peak = max(y1);
    pdfDat(2).peak = max(y2);
    
    pdfDat(1).median = median(y1);
    pdfDat(2).median = median(y2);
    
    pdfDat(1).halfwidth = fwhm;
    pdfDat(2).halfwidth = fwhm2;
    
    pdfDat(1).xi=xi;
    pdfDat(2).xi=xi;
    
    pdfDat(1).AOC = aoc1;
    pdfDat(2).AOC = aoc2;
    
    pdfDat(1).objSingle = objA1;
    pdfDat(2).objSingle = objA2;
    pdfDat(1).objComb = objA;
    pdfDat(2).objComb = objA;
    
else
    comp = 1;
    
    mdl= fitgmdist(dat,comp,'Options',options); % fit the mixed gaussian
    
    figure(112); % inti fig
    paramEsts= mdl;
    %mu of mixture of gaussians of 2 components % gaussian parameters
    MU=[paramEsts.mu(1)];
    SIGMA=[paramEsts.Sigma(1)];
    PPp=[paramEsts.PComponents(1)];
    
    objA = gmdistribution(MU,SIGMA,PPp); % creating fit objects
    xi = linspace(min(dat(:)), max(dat(:)),100); % index var
    
    y1 = pdf(objA,xi.').*PPp(1); % first pdf (low perf?)
    
    skew = skewness(regions(regions > 0));
    %     y2 = pdf(objA2,xi.').*PPp(2); % (high perf)?
    %     ycomb = pdf(objA,xi.'); % combined pdf (looks like histogram)
    SIGMA
    %begin plotting
    % plot(xi,ycomb,'r-','linewidth',2)
    % hold on;
    plot(xi,y1,'g-','Linewidth',2);
    % plot(xi,y2,'b-','Linewidth',2);
    hold on;
    histogram(log10(regions(regions > 0)),'Normalization','pdf','NumBins',numel(xi),'FaceAlpha',0);
    lgd=legend(sprintf('Mu = %0.2f, Sigma = %0.3f',double(MU),SIGMA(end)),'Histogram');
    lgd.FontSize = 16;
    title(animal);
    xlim([0 2.5]);
    ylim([0 6]);
    xlabel('Log_{10}( Perfusion [mL / min / 100 mL]');
    ylabel('PDF');
    hold off;
    saveas(gcf,fname,'png'); % saving historam figure
%     pause;
    aoc1=trapz(xi,pdf(objA,xi.').*PPp(1)); % getting area under curves
    
    
    
     hh1 = (min(y1) + max(y1)) / 2;
%     hh2 = (min(2) + max(y2)) / 2;
    
    xin1 = find(y1 >= hh1,1,'first');
    xin2 = find(y1 >= hh1,1,'last');
    
    x1 = xi(xin1);
    x2 = xi(xin2);
    fwhm = x2-x1;
  
    
    %create output var
    pdfDat(1).y=y1;
   
    pdfDat(1).xi=xi;
   
    pdfDat(1).AOC = aoc1;

    pdfDat(1).halfwidth = fwhm;

    pdfDat(1).peak = max(y1);
    pdfDat(1).median = median(y1);
    pdfDat(1).objSingle = objA;
    
    
    
end
end

