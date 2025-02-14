%% Figure 6 & Figure S5

% This is the code used to report C and C45 in the paper. 
% Fit a decision tree based on one speed and 
% decoded direction for the other speeds.

%% load data
load ../dataFiles/cellData_NPX_speed.mat; 

% columns: 27 neurons
% rows 1-160   : dT = 100 ms(the slowest speed) 
% rows 161-320 : dT = 50  ms. 
% rows 321-480 : dT = 25  ms. 
% rows 481-640 : dT = 8.3 ms. 
% the data is organized such that each block of speed is organized from 0 to 315° 
%(e.g., for slowest speed rows 1-20 is 0°,  21-40 is 45°, etc.). 

% speed 
RFd = 2.73; 
step = RFd/6;
dT = [0.1, 0.05, 0.025, 0.0083]; 
speeds = step./dT;


% trainData
% this is the training data with dt = 25 ms -- its the most segregated data
trainData = cellData_NPX_speed(321:480, :);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Direction Index (DI) from trainData. 
% No baseline activity available

% direction tuning
for i = 1:8
    dir_tuning(i,:) = mean(trainData((i-1)*20+1:i*20,:), 1);
end

% find preferred, non-preferred direction
[M, prefI] = max(dir_tuning, [], 1);   % preferred direction index
nprefI = rem(prefI + 3, 8) + 1;        % non-preferred direction index

% compute DI
for i = 1:size(dir_tuning,2)
    di_vec(i) = (dir_tuning(prefI(i), i) - dir_tuning(nprefI(i), i)) ...
                /dir_tuning(prefI(i),i);
end


%% train and fit

% DI threshold
thr = [0.0, 0.245, 0.495 0.745];

% Lists for output 
perfD = [];     % report C for Dir
perf45D = [];   % report C45 for Dir
perfND = [];    % report C for nonDir
perf45ND = [];  % report C45 for nonDir

numDir = [];    % number of directional units

% colormap setting for figure
nreps = 20;
clrs = colormap(jet);
clrs = flipud(vertcat(clrs, flipud(clrs)));
nJ =  55;
close all; 

% Train and Fit
for t1 = 1:length(thr)
    
    % directional unit vs. non-directional unit
    dir = find(di_vec >= thr(t1));
    nDir = find(di_vec < thr(t1));
    
    %number of direction selective neurons for the different thresholds tested.
    numDir(t1) = length(dir); 

    % run PCA
    [cfDir, sDir, ~,~, expl, mu] = pca(trainData(:, dir));
    if(~isempty(nDir))
        [cfNdir, sNdir, ~,~,exp2, mu2] = pca(trainData(:, nDir));
    end
    idx = min(4, size(cfDir, 2)); % ideally we want 4 dims but for some thresholds this is smaller

    % Now lets project the other speeds onto the PCA based on direction selective cells
    testData1 = cellData_NPX_speed(1:160, dir);
    testData2 = cellData_NPX_speed(161:320, dir);
    testData3 = cellData_NPX_speed(481:640, dir);

    % First project data on to idx PCs
    testS1 = (testData1-mu)*cfDir(:,1:idx);
    testS2 = (testData2-mu)*cfDir(:,1:idx);
    testS3 = (testData3-mu)*cfDir(:,1:idx);
    

    % projection of Data on PC space
    if thr(t1)==0.245 % this is Figure S5B
        figure; 
        set(gcf,'position',[100,100,250,250]);  
        sgtitle('Figure S5B');         

        subplot(1,1,1);  hold on;
        for i = 1:8
            if i ~= 7
                scatter(sDir((i-1)*nreps+1:i*nreps,1), sDir((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
            else
                scatter(sDir((i-1)*nreps+1:i*nreps,1), sDir((i-1)*nreps+1:i*nreps,2), 20, [0.2 0.7 0.2], 'filled');
            end
        end
        axis([-50 100 -50 50]);
        legend({'0°','45°','90°','135°','180°','225°','270°','315°'},'position',[0.75,0.3,0.05,0.15]); 
        set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
        xlabel('PC1');   ylabel('PC2'); 
        pbaspect([1 1 1]); 
        title(sprintf("Speed: %.1f°/s Threshold: %.2f Directional",speeds(3), round(thr(t1)*100)/100)); 
    end


    if thr(t1)==0.495 % this is Figure 6
        figure; 
        set(gcf,'position',[100,100,800,500]); 
        sgtitle('Figure 6');                 

        subplot(2,3,1);   hold on;
        for i = 1:8
            if i ~= 7
                scatter(sDir((i-1)*nreps+1:i*nreps,1), sDir((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
            else
                scatter(sDir((i-1)*nreps+1:i*nreps,1), sDir((i-1)*nreps+1:i*nreps,2), 20, [0.2 0.7 0.2], 'filled');
            end
        end
        axis([-50 100 -50 50]);
        legend({'0°','45°','90°','135°','180°','225°','270°','315°'},'position',[0.31,0.65,0.05,0.15]); 
        set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
        xlabel('PC1');   ylabel('PC2'); 
        pbaspect([1 1 1]); 
        title(sprintf("Speed: %.1f°/s Threshold: %.2f Directional",speeds(3), round(thr(t1)*100)/100)); 

        if(~isempty(nDir))
            subplot(2,3,2);   hold on;
            axis([-50 100 -50 50]);
            for i = 1:8
                if i ~= 7
                    scatter(sNdir((i-1)*nreps+1:i*nreps,1), sNdir((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
                else
                    scatter(sNdir((i-1)*nreps+1:i*nreps,1), sNdir((i-1)*nreps+1:i*nreps,2), 20, [0.2 0.7 0.2], 'filled');
                end
            end
            set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
            title(sprintf("Non Directional")); 
            xlabel('PC1');
            pbaspect([1 1 1]); 
        end        

        subplot(2,3,4);   hold on; 
        axis([-50 100 -50 50]);
        c7 = [0.2 0.7 0.2];
        for i = 1:8
            if i ~= 7
                scatter(testS1((i-1)*nreps+1:i*nreps,1), testS1((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
            else
                scatter(testS1((i-1)*nreps+1:i*nreps,1), testS1((i-1)*nreps+1:i*nreps,2), 20,c7, 'filled');
            end
        end
        title(sprintf("Speed: %.1f°/s",speeds(1))); 
        pbaspect([1 1 1]); 
        set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
        xlabel('PC1');   ylabel('PC2'); 
            
        subplot(2,3,5);   hold on; 
        axis([-50 100 -50 50]);
        for i = 1:8
            if i ~= 7
                scatter(testS2((i-1)*nreps+1:i*nreps,1), testS2((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
            else
                scatter(testS2((i-1)*nreps+1:i*nreps,1), testS2((i-1)*nreps+1:i*nreps,2), 20, c7, 'filled');
            end
        end
        title(sprintf("Speed: %.1f°/s",speeds(2))); 
        pbaspect([1 1 1]); 
        set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
        xlabel('PC1');   ylabel('PC2');         
        
        subplot(2,3,6);   hold on; 
        axis([-50 100 -50 50]);
        for i = 1:8
            if i ~= 7
                scatter(testS3((i-1)*nreps+1:i*nreps,1), testS3((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
            else
                scatter(testS3((i-1)*nreps+1:i*nreps,1), testS3((i-1)*nreps+1:i*nreps,2), 20, c7, 'filled');
            end
        end
        title(sprintf("Speed: %.1f°/s",speeds(4))); 
        pbaspect([1 1 1]); 
        set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
        xlabel('PC1');   ylabel('PC2');         

    end

    % Make classification tree based on first N principal components
    numPCs = idx;
    mdl = fitctree(sDir(:,1:numPCs),[repmat(0,20,1); repmat(45,20,1); repmat(90,20,1); repmat(135,20,1); repmat(180,20,1); repmat(225,20,1); repmat(270,20,1); repmat(315,20,1)]);

    btruth = [repmat(0,20,1); repmat(45,20,1); repmat(90,20,1); repmat(135,20,1); repmat(180,20,1); repmat(225,20,1); repmat(270,20,1); repmat(315,20,1)]; 

    dir_pred1 = predict(mdl,testS1(:,1:numPCs));
    dir_pred2 = predict(mdl,testS2(:,1:numPCs)); 
    dir_pred3 = predict(mdl,testS3(:,1:numPCs)); 

    acc1 = sum(dir_pred1 == btruth)/160; 
    acc2 = sum(dir_pred2 == btruth)/160; 
    acc3 = sum(dir_pred3 == btruth)/160;

    perfD(t1,:) = [acc1 acc2 acc3]; %this is C performance as a function of threshold
    fprintf("S1: %.2f  S2: %.2f  S3: %.2f\n",acc1,acc2,acc3)

    diff1 = abs(dir_pred1 - btruth); 
    aacc1 = length(find(diff1 <= 45))/160; 

    diff2 = abs(dir_pred2 - btruth); 
    aacc2 = length(find(diff2 <= 45))/160; 

    diff3 = abs(dir_pred3 - btruth); 
    aacc3 = length(find(diff3 <= 45))/160; 
    
    perf45D(t1,:) = [aacc1 aacc2 aacc3]; % this is C45 performance as a function of threshold
    
    fprintf("S1: %.2f  S2: %.2f  S3: %.2f\n",aacc1,aacc2,aacc3)

    if thr(t1)==0.495 % this is for Figure 6
        subplot(2,3,4)
        t = text(0,-45,sprintf('Acc: %.2f, Acc < 45: %.2f',acc1,aacc1),'fontname','arial','fontsize',12);
    
        subplot(2,3,5)
        t = text(0,-45,sprintf('Acc: %.2f, Acc < 45: %.2f',acc2,aacc2),'fontname','arial','fontsize',12);
    
        subplot(2,3,6)
        t = text(0,-45,sprintf('Acc: %.2f, Acc < 45: %.2f',acc3,aacc3),'fontname','arial','fontsize',12);

        disp('Done generating Figure 6');         
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now for non-dir selective cells 
    %
    % Let's project the other speeds onto the PCA based on direction
    % selective cells


    if(~isempty(nDir))
        testData1 = cellData_NPX_speed(1:160, nDir);
        testData2 = cellData_NPX_speed(161:320, nDir);
        testData3 = cellData_NPX_speed(481:640, nDir);

        % First project data on to idx PCs
        idx = min(4, size(cfNdir, 2)); % ideally we want 4 dims but for some thresholds this is smaller
        testS1 = (testData1-mu2)*cfNdir(:,1:idx);
        testS2 = (testData2-mu2)*cfNdir(:,1:idx);
        testS3 = (testData3-mu2)*cfNdir(:,1:idx);

        numPCs = idx;
        mdl = fitctree(sNdir(:,1:numPCs),[repmat(0,20,1); repmat(45,20,1); repmat(90,20,1); repmat(135,20,1); repmat(180,20,1); repmat(225,20,1); repmat(270,20,1); repmat(315,20,1)]);

        btruth = [repmat(0,20,1); repmat(45,20,1); repmat(90,20,1); repmat(135,20,1); repmat(180,20,1); repmat(225,20,1); repmat(270,20,1); repmat(315,20,1)]; 
        
        dir_pred1 = predict(mdl,testS1(:,1:numPCs));
        dir_pred2 = predict(mdl,testS2(:,1:numPCs)); 
        dir_pred3 = predict(mdl,testS3(:,1:numPCs)); 

        acc1 = sum(dir_pred1 == btruth)/160; 
        acc2 = sum(dir_pred2 == btruth)/160; 
        acc3 = sum(dir_pred3 == btruth)/160;
        perfND(t1,:) = [acc1 acc2 acc3]; %this is C performance as a function of threshold
        fprintf("Acc1: %.2f Acc2: %.2f  Acc3: %.2f",acc1,acc2,acc3)


        diff1 = abs(dir_pred1 - btruth); 
        aacc1 = length(find(diff1 <= 45))/160; 

        diff2 = abs(dir_pred2 - btruth); 
        aacc2 = length(find(diff2 <= 45))/160; 

        diff3 = abs(dir_pred3 - btruth); 
        aacc3 = length(find(diff3 <= 45))/160; 

        fprintf("Acc1: %.2f Acc2: %.2f  Acc3: %.2f",aacc1,aacc2,aacc3); 

        perf45ND(t1,:) = [aacc1 aacc2 aacc3]; % this is C45 performance as a function of threshold
        
        if thr(t1)==0.495   % this is for Figure S5A
            figure; 
            set(gcf,'position',[100,100,800,250]);             
            sgtitle('Figure S5A');                             

            subplot(1,3,1);   hold on; 
            axis([-50 100 -50 50]);
            for i = 1:8
                if i ~= 7
                    scatter(testS1((i-1)*nreps+1:i*nreps,1), testS1((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
                else
                    scatter(testS1((i-1)*nreps+1:i*nreps,1), testS1((i-1)*nreps+1:i*nreps,2), 20, c7, 'filled');
                end
            end
            title(sprintf("Speed: %.1f°/s Threshold: %.2f Non Directional",speeds(1), round(thr(t1)*100)/100)); 
            t = text(0,-45,sprintf('Acc: %.2f  Acc within 45°: %.2f',acc1,aacc1)); 
            pbaspect([1 1 1]); 
            set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
            legend({'0°','45°','90°','135°','180°','225°','270°','315°'},'position',[0.31,0.4,0.05,0.15]); 
            xlim([-50 100]); ylim([-50 50]); 
            xlabel('PC1');   ylabel('PC2'); 
            
            subplot(1,3,2);   hold on;
            axis([-50 100 -50 50]);
            for i = 1:8
                if i ~= 7
                    scatter(testS2((i-1)*nreps+1:i*nreps,1), testS2((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
                else
                    scatter(testS2((i-1)*nreps+1:i*nreps,1), testS2((i-1)*nreps+1:i*nreps,2), 20, c7, 'filled');
                end
            end
            title(sprintf("Speed: %.1f°/s",speeds(2))); 
            t = text(0,-45,sprintf('Acc: %.2f  Acc within 45°: %.2f',acc2,aacc2)); 
            pbaspect([1 1 1]); 
            set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[])
            xlim([-50 100]); ylim([-50 50]);             

            subplot(1,3,3);   hold on;
            axis([-50 100 -50 50]);
            for i = 1:8
                if i ~= 7
                    scatter(testS3((i-1)*nreps+1:i*nreps,1), testS3((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
                else
                    scatter(testS3((i-1)*nreps+1:i*nreps,1), testS3((i-1)*nreps+1:i*nreps,2), 20, c7, 'filled');
                end
            end
            title(sprintf("Speed: %.1f°/s",speeds(4)));
            t = text(0,-45,sprintf('Acc: %.2f  Acc within 45°: %.2f',acc3,aacc3)); 
            pbaspect([1 1 1]); 
            set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[])
            xlim([-50 100]); ylim([-50 50]); 
        end
    end

end

%% Figure S5C
figure; 
set(gcf,'position',[100 100 280 200]); 
sgtitle('Figure S5C');                 

speed = speeds([1,2,4]); 
xticknames = [{['0.0 (',num2str(numDir(1)),')']},...
              {['0.25 (',num2str(numDir(2)),')']},...
              {['0.5 (',num2str(numDir(3)),')']},...
              {['0.75 (',num2str(numDir(4)),')']}]; 

legend_mtx1 = []; 
legend_mtx2 = []; 

subplot(1,1,1); 
for n=1:length(speed)
    plot(thr, perfD(:,n), 'o-',...
        'color', 0.4*(3-n)*ones(1,3), 'MarkerFaceColor',0.4*(3-n)*ones(1,3), 'LineWidth', 2); 
    hold on; 
    legend_mtx1 = [legend_mtx1; {['Motion speed: ',num2str(speed(n)),'°/s']}]; 

    plot(thr, perf45D(:,n), '^:',...
        'color', 0.4*(3-n)*ones(1,3), 'MarkerFaceColor',0.4*(3-n)*ones(1,3), 'LineWidth', 2); 
    legend_mtx2 = [legend_mtx2; {['Motion speed: ',num2str(speed(n)),'°/s']}]; 
   
end
xticks(thr);
xticklabels(xticknames);
%legend([legend_mtx1; legend_mtx2],'Location','southeast'); 
xlabel('DI threshold (#units)');
ylabel('Performance correct');
xlim([-0.05,0.8]);
ylim([0.1,0.85]); 
set(gca,'TickDir','out','box','off');


%% Figure S5D
% D. Low dimensional representation of responses to surface motion 
%    yields limited patterning. 

load ../dataFiles/cellData_NPX_FigS5DE.mat; 

figure; 
set(gcf,'position',[100,100,600,250]);  
sgtitle('Figure S5DE');                         


% make train dataset (Surface motion: dX = 5)
test = []; 
for dX =1:5  % for 5 stimulus types
    trainData = []; 
    DI_base = []; 
    for i = 1:length(cellData_NPX_FigS5DE) %this is the number of neurons
        % taking 8 columns corresponding to to the different stim conditions
        % taking 19 trials 
        col_start = (dX-1)*8 + 1; 
        col_end = dX*8; 
        d1 = cellData_NPX_FigS5DE(i).respMtx(1:19,col_start:col_end); 
        nreps = size(d1, 1); %this is the number of repeats    
    
        % accumulate trainData
        trainData(:, i) = reshape(d1, 1, []);
    
        % accumulate DI values for Surface motion
        DI_base = [DI_base; cellData_NPX_FigS5DE(i).DI_base]; 
    end
    test(dX,:,:) = trainData;    

    if dX==5  % this is for Figure S5D
        % direction selective neurons
        neuid = find(DI_base(:,dX) > 0.5);
        
        % run PCA
        [cfDir, sDir, ~,~, expl, mu] = pca(trainData(:,neuid));
        
        subplot(1,2,1);  hold on;
        
        for i = 1:8    
            if i ~= 7
                scatter(sDir((i-1)*nreps+1:i*nreps,1), sDir((i-1)*nreps+1:i*nreps,2), 20, clrs(i*nJ, :), 'filled');
            else
                scatter(sDir((i-1)*nreps+1:i*nreps,1), sDir((i-1)*nreps+1:i*nreps,2), 20, [0.2 0.7 0.2], 'filled');
            end
        end
        axis([-50 100 -50 50]);
        legend({'0°','45°','90°','135°','180°','225°','270°','315°'},'position',[0.41,0.3,0.05,0.15]); 
        set(gca,'tickdir','out','fontname','arial','xticklabels','','yticklabels','','xtick',[],'ytick',[]); 
        xlabel('PC1'); ylabel('PC2');
        pbaspect([1 1 1]);
        title(['Surface Motion: ',num2str(length(neuid)),' neurons with DI > 0.5']);
    end
end

% Figure S5E
% E. Population decoding performance for object and surface motion 
%    direction calculated for four different DI thresholds.

thr = [0.0 0.245 0.495 0.745];

for t1 = 1:4          % for threshold levels
    for dX = [3 5]    % two trainsets: object (RF/12), surface

        neuid = find(DI_base(:,dX)> thr(t1));
        idx = min(4, length(neuid));
        numPCs = idx;

        trainData = squeeze(test(dX,:,:));
        [cfDir, sDir, ~,~, expl, mu] = pca(trainData(:,neuid));
        mdl = fitctree(sDir(:,1:numPCs),[repmat(0,19,1); repmat(45,19,1); repmat(90,19,1); repmat(135,19,1); repmat(180,19,1); repmat(225,19,1); repmat(270,19,1); repmat(315,19,1)]);
        btruth = [repmat(0,19,1); repmat(45,19,1); repmat(90,19,1); repmat(135,19,1); repmat(180,19,1); repmat(225,19,1); repmat(270,19,1); repmat(315,19,1)]; 
        %project data
        
        for j = 1:5 % test all the conditions now for a fit
            testData = squeeze(test(j,:,:));
            testS1 = (testData(:,neuid)-mu)*cfDir(:,1:idx);   
            
            %%%then predict
            dir_pred1 = predict(mdl,testS1(:,1:numPCs));
            perC(t1,dX,j)= sum(dir_pred1 == btruth)/152; 
        
            %within 45 deg
            diff1 = abs(dir_pred1 - btruth); 
            perC45(t1,dX,j) = length(find(diff1 <= 45))/152; 
        end
    end
end

subplot(1,2,2); 
xticknames = [{'0.0'},{'0.25'},{'0.5'},{'0.75'}]; 
plot(thr, perC(:,3,4),'o-',...
        'color', [0,0,0], 'MarkerFaceColor',[0,0,0], 'LineWidth', 2);  
hold on;
plot(thr, perC(:,3,5),'o-',...
        'color', [1,0,0], 'MarkerFaceColor',[1,0,0], 'LineWidth', 2);  
plot(thr, perC(:,5,3),'o-',...
        'color', [0,0,1], 'MarkerFaceColor',[0,0,1], 'LineWidth', 2);  
plot(thr, perC(:,5,4),'o-',...
        'color', [0.5,0.5,1], 'MarkerFaceColor',[0.5,0.5,1], 'LineWidth', 2);  

xticks(thr);
xticklabels(xticknames);
xlabel('DI threshold');
ylabel('Performance correct');
title('Exact decoding');
xlim([-0.05,0.8]);
ylim([0,0.75]); 
set(gca,'TickDir','out','box','off');

clearvars -except cellData_NPX_*; 