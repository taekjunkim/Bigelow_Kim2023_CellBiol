%% Code to make_FigureS3

% A. Correlation between RF/36 and other dX conditions: boxplot, histogram
% B. significant DI matrix
% C. Histogram of mean aggregated correlation
% D. raster plot: DI (object) vs. DI (surface)
% - RF/3, RF/6, RF/12, RF/36


%% A. Correlation between RF/36 and other dX conditions: boxplot, histogram
load ../dataFiles/cellData_NPX_dX.mat

conditions = [{'RF/3'},{'RF/6'},{'RF/12'},{'RF/36'},{'Local'}]; 
corrMtx = []; 
pvalMtx = []; 

for c=1:length(cellData_NPX_dX)
    % RF/36 condition is the reference direction tuning
    ref_tuning = cellData_NPX_dX(c).dir_tuning(:,4); 
    for n=1:5
        if n==4
            corrMtx(c,n) = 1; 
            pvalMtx(c,n) = 0; 
            continue; 
        end
        [r,p] = corr(cellData_NPX_dX(c).dir_tuning(:,n),ref_tuning); 
        corrMtx(c,n) = r; 
        pvalMtx(c,n) = p; 
    end
end

figure; 
set(gcf,'Position',[100 100 800 300]); 
sgtitle('Figure S3A');              

subplot(1,2,1); 
boxplot(corrMtx);   hold on; 
for n=1:5
    plot(n,median(corrMtx(:,n)),'ro');
end
xticklabels(conditions); 
set(gca,'box','off','TickDir','out','Ylim',[-1,1]); 

for n=1:4
    subplot(1,8,n+4);
    if n<4
        histogram(corrMtx(:,n),-1:0.05:1,'FaceColor',[0.5,0.5,0.5], ...
                  'EdgeColor',[1,1,1],'FaceAlpha',1);    hold on;
        signi = find((corrMtx(:,n)>0) & pvalMtx(:,n)<0.05); 
        histogram(corrMtx(signi,n),-1:0.05:1,'FaceColor',[0,0,0],...
                  'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;
    else
        histogram(corrMtx(:,n+1),-1:0.05:1,'FaceColor',[0.5,0.5,0.5],...
                  'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;
        signi = find((corrMtx(:,n+1)>0) & pvalMtx(:,n+1)<0.05); 
        histogram(corrMtx(signi,n+1),-1:0.05:1,'FaceColor',[0,0,0],...
                  'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;
    end
    title(['N = ',num2str(length(signi))]);
    if n>1
        xticklabels([]);
    end
    set(gca,'box','off','TickDir','out','Xlim',[-1, 1]); 
    set(gca,'View',[90 -90]); 
end

%% B. significant DI matrix
dft_di_base = []; 
ppp = []; 
for i=1:78
    dft_di_base(i,:) = cellData_NPX_dX(i).DI_base;     
    ppp(i,:) = cellData_NPX_dX(i).anova(2,:); 
end
ppp(find(ppp(:)<0.05)) = 1; 
ppp(find(ppp(:)<1)) = 0; 

figure; 
set(gcf,'Position',[100 100 400 150]);
sgtitle('Figure S3B');              

highDI = zeros(78,5); 
highDI(find(dft_di_base(:)>=0.5)) = 1; 

sigMtx = sortrows(highDI.*ppp,[-4,-3,-2,-1,-5]);

subplot(1,1,1); 
imagesc(1-sigMtx');   colormap("gray"); 
yticks([1,2,3,4,5]); 
yticklabels([{'RF/3'},{'RF/6'},{'RF/12'},{'RF/36'},{'Drift'}]);
xlabel('Units'); 
ylabel('dX condition'); 
title('DI ≥ 0.5, P < 0.05 (one-way ANOVA)'); 
set(gca,'TickDir','out');


%% Figure S3C: Histogram of mean aggregated correlation
num_units = length(cellData_NPX_dX); 
pairs = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4]; 

consistency = []; 
for n=1:num_units
    numerators = []; 
    denominators = []; 
    for p=1:length(pairs)
        tuning1 = cellData_NPX_dX(n).dir_tuning(:,pairs(p,1));  
        tuning2 = cellData_NPX_dX(n).dir_tuning(:,pairs(p,2)); 

        % mean subtracted
        tuning1 = tuning1 - mean(tuning1); 
        tuning2 = tuning2 - mean(tuning2);         

        cov_mtx = cov(tuning1,tuning2); 
        numerators = [numerators; cov_mtx(1,2)]; 
        denominators = [denominators; std(tuning1)*std(tuning2)]; 
    end
    DS_consistency(n).numerators = numerators; 
    DS_consistency(n).denominators = denominators; 

    DS_consistency(n).consistency = sum(numerators)/sum(denominators); 

    % 100 simulation to compute statistical significance
    for r=1:100
        numerators = []; 
        denominators = []; 
        for p=1:length(pairs)
            tuning1 = cellData_NPX_dX(n).dir_tuning(:,pairs(p,1));  
            tuning2 = cellData_NPX_dX(n).dir_tuning(:,pairs(p,2)); 
    
            % mean subtracted
            tuning1 = tuning1 - mean(tuning1); 
            tuning2 = tuning2 - mean(tuning2);         
    
            tuning1_shuffled = tuning1(randperm(length(tuning1)));
            tuning2_shuffled = tuning2(randperm(length(tuning2)));

            cov_mtx = cov(tuning1_shuffled,tuning2_shuffled); 
            numerators = [numerators; cov_mtx(1,2)]; 
            denominators = [denominators; std(tuning1)*std(tuning2)]; 
        end
        DS_consistency(n).random_consistency(r) = sum(numerators)/sum(denominators); 
    end
    bigger_than_random = length(find(DS_consistency(n).random_consistency>DS_consistency(n).consistency)); 
    if bigger_than_random < 5;
        DS_consistency(n).signi = 1; 
    else
        DS_consistency(n).signi = 0; 
    end

    consistency = [consistency; DS_consistency(n).consistency DS_consistency(n).signi]; 
end

figure; 
set(gcf,'Position',[100 100 250 150]);
sgtitle('Figure S3C');              

subplot(1,1,1); 
histogram(consistency(:,1),-1:0.05:1,'Facecolor',[0.8,0.8,0.8],...
          'Edgecolor',[1,1,1],'FaceAlpha',1,'EdgeAlpha',1);   hold on; 
signi = find(consistency(:,2)==1); 
histogram(consistency(signi,1),-1:0.05:1,'Facecolor',[0,0,0],...
          'EdgeColor',[1,1,1],'FaceAlpha',1,'EdgeAlpha',1);   
plot(mean(consistency(:,1)),13,'kv','MarkerFaceColor','k'); 
set(gca,'TickDir','out','box','off');
xlabel('Mean aggregated correlation');
ylabel('Number of units');


%% Figure S3D: raster plot: DI (object) vs. DI (surface)
% - RF/3, RF/6, RF/12, RF/36

dft_di = []; 
ppp = []; 
for i=1:78
    dft_di(i,:) = cellData_NPX_dX(i).DI_base; 
    ppp(i,:) = cellData_NPX_dX(i).anova(2,:); 
end

ppp(find(ppp(:)<0.05)) = 1; 
ppp(find(ppp(:)<1)) = 0; 


figure; 
set(gcf,'Position',[100 100 500 500]); 
sgtitle('Figure S3D');          

dX_cond = [{'RF/3'},{'RF/6'},{'RF/12'},{'RF/36'}]

for n=1:4   % RF/3, RF/6, RF/12, RF/36
    subplot(2,2,n);  % RF/3
    plot(dft_di(:,n), dft_di(:,5),'.','color',[0.5 0.5 0.5],'MarkerSize',15);   hold on;
    [r, p] = corr(dft_di(:,n),dft_di(:,5))
    plot([0,2],[0,2],'k')
    title(dX_cond(n),['r = ',num2str(round(r,2)),', p = ',num2str(round(p,2))]); 
    xticks([0,0.5,1,1.5,2]); 
    yticks([0,0.5,1,1.5,2]); 
    xlabel('Object motion DI'); 
    ylabel('Surface motion DI'); 
    set(gca,'box','off','TickDir','out')
end

clearvars -except cellData_NPX_dX; 


