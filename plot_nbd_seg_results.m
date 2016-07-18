function plot_nbd_seg_results(displaycharacter, file)

load(file)

yaxistexts = {'Distance / km', 'Neighbourhood Percentage Non-Western', 'Neighbourhood Mean Disposable Income / 10k SEK/yr', 'Neighbourhood Percentage With Children'};

%displaycharacter = 2; %change this to display utility as a function of the neighbourhood characteristics above in yaxistexts

yaxistext = yaxistexts{displaycharacter};

colorscale = [min(f{displaycharacter}(:)), max(f{displaycharacter}(:))];


figure
colormap(hot(128))
subplot(2,2,1)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(1, :, 1,1,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(1, :, 1,1, :))' - 1.96*squeeze(se{displaycharacter}(1, :, 1,1, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(1, :, 1,1, :))' + 1.96*squeeze(se{displaycharacter}(1, :, 1,1, :))' < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');

axis xy
%ylabel(yaxistext)
%xlabel('Disposable income')
title('(A) Swedish Young No Child')
if displaycharacter ==1
    set(gca, 'YLim', [0, 5e5], 'YTick', linspace(0, 5e5, 6), 'YTickLabel', linspace(0, 500, 6))
    
end

colorbar

subplot(2,2,2)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(2, :, 1,1,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(2, :, 1,1, :))' - 1.96*squeeze(se{displaycharacter}(2, :, 1,1, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(2, :, 1,1, :))' + 1.96*squeeze(se{displaycharacter}(2, :, 1,1, :))' < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');

axis xy
%ylabel(yaxistext)
%xlabel('Disposable income')
title('(B) Non-Swedish Young, No Child')
colorbar
if displaycharacter ==1
    set(gca, 'YLim', [0, 5e5], 'YTick', linspace(0, 5e5, 6), 'YTickLabel', linspace(0, 500, 6))
    
end

subplot(2,2,3)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(1, :, 2,2,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(1, :, 2,2, :))' - 1.96*squeeze(se{displaycharacter}(1, :, 2,2, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(1, :, 2,2, :))' + 1.96*squeeze(se{displaycharacter}(1, :, 2,2, :))' < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');
axis xy
%ylabel(yaxistext)
%xlabel('Disposable income')
title('(C) Swedish Middle-Aged w/ Child')
colorbar
if displaycharacter ==1
    set(gca, 'YLim', [0, 5e5], 'YTick', linspace(0, 5e5, 6), 'YTickLabel', linspace(0, 500, 6))
    
end

subplot(2,2,4)
imagesc(X{2}, Y{displaycharacter}, squeeze(f{displaycharacter}(2, :, 2,2,:))')
caxis(colorscale);
hold on
Zhigh = double(squeeze(f{displaycharacter}(2, :, 2,2, :))' - 1.96*squeeze(se{displaycharacter}(2, :, 2,2, :))' > 0);
Zlow = double(squeeze(f{displaycharacter}(2, :, 2,2, :))' + 1.96*squeeze(se{displaycharacter}(2, :, 2,2, :))' < 0);
contour(X{2}, Y{displaycharacter}, Zhigh, [1 1],'k');
contour(X{2}, Y{displaycharacter}, Zlow, [1 1],'w');
axis xy
%ylabel(yaxistext)
%xlabel('Disposable income')
title('(D) Non-Swedish Middle-Aged w/ Child')
colorbar
if displaycharacter ==1
    set(gca, 'YLim', [0, 5e5], 'YTick', linspace(0, 5e5, 6), 'YTickLabel', linspace(0, 500, 6))
    
end

suplabel('Disposable Income / 10k SEK/yr', 'x');
suplabel(yaxistext, 'y');
