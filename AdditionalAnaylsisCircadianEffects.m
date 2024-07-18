clear all
close all

% WBDataRawExp1 = readtable('Y:\niels\Jianfeng\WesternBlotResults\WB_raw_26wells12wells-1_AllAnimals.csv');
% WBDataRawExp1 = readtable('Y:\niels\Jianfeng\WesternBlotResults\WB_raw_26wells12wells-1.csv');
WBDataRawExp1 = readtable('Y:\niels\Jianfeng\WesternBlotResults\WB_raw_26wells12wells-1_NoNorm.csv');
WBAnimalNamesExp1 = table2array(WBDataRawExp1(:,1));
WBDataExp1.HypoThalamus = table2array(WBDataRawExp1(:,2:5));
WBDataExp1.Cortex = table2array(WBDataRawExp1(:,8:11));
WBMeasureNamesTmp = fieldnames(WBDataRawExp1);
WBMeasureNamesExp1.Hypo = WBMeasureNamesTmp(2:5,1);
WBMeasureNamesExp1.Cortex = WBMeasureNamesTmp(8:11,1);

% DataExp1 = WBDataExp1.HypoThalamus(1:16,:)./repmat(mean(WBDataExp1.HypoThalamus(17:32,:),1),[16,1]);
DataExp1 = WBDataExp1.HypoThalamus(1:10,:)./repmat(mean(WBDataExp1.HypoThalamus(11:20,:),1),[10,1]);
% DataExp1 = WBDataExp1.HypoThalamus(1:10,:)*repmat(mean(WBDataExp1.HypoThalamus(1:10,:),1),[10,1]);
% DataExp1 = WBDataExp1.HypoThalamus(1:10,:);

% strfind(WBAnimalNamesExp1,'Sleep')

% WBDataRawExp2 = readtable('Y:\niels\LunAMPAScaling\NewDataSet\RSD_WB_Data_Cort_Hypo_Merg.csv');
WBDataRawExp2 = readtable('Y:\niels\Jianfeng\WesternBlotResults\Exp2_WB_raw_26wells12wells-1_NoNorm.csv');
WBAnimalNamesExp2 = table2array(WBDataRawExp2(:,1));
WBDataExp2.HypoThalamus = table2array(WBDataRawExp2(:,2:5));
WBDataExp2.Cortex = table2array(WBDataRawExp2(:,8:11));
WBMeasureNamesTmp = fieldnames(WBDataRawExp2);
WBMeasureNamesExp2.Hypo = WBMeasureNamesTmp(2:5,1);
WBMeasureNamesExp2.Cortex = WBMeasureNamesTmp(8:11,1);

DataExp2 = WBDataExp2.HypoThalamus(1:8,:)./repmat(mean(WBDataExp2.HypoThalamus(9:16,:),1),[8,1]);
% DataExp2 = WBDataExp2.HypoThalamus(1:8,:);

hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 20 500 500])
set(gcf,'color','white')
set(gcf,'name',' Protein levels Exp1 vs Exp2 Hypothalamus','NumberTitle','off')
tiledlayout(1,4)
for iPlot = 1:4
    nexttile

    b = bar([1 2], [mean(DataExp1(:,iPlot),1), mean(DataExp2(:,iPlot),1)]);
    b.FaceColor = 'flat';
    b.CData([1],:) = [1 0 0];
    b.CData([2],:) = [0 0 1];

    ylabel('Sleep/Wake or TSD');
    Labels = {'Exp1','Exp2'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca, 'box', 'off')
    hold on
    StandardError = [];
    tmp =  [];
    %adding data points for distribution
    for barIndex = 1: size(b.YData,2)
        if barIndex ==1
            Data = DataExp1(:,iPlot);
        end
        if barIndex ==2
            Data = DataExp2(:,iPlot);
        end

        sorted_data_y = sort(Data, 'ascend');


        sorted_data_x = [];
        [N,edges] = histcounts(Data(:,1),100);
        for iBin = 1: size(N,2)
            if N(1,iBin) >0
                %             if mod(N(1,iBin),2)== 0 %gerade Zahlen
                XtmpLocation = [-(N(1,iBin)-1)/2:(N(1,iBin)-1)/2];
                %             end
                if N(1,iBin) == 1
                    XtmpLocation = 0;
                end
                %             if mod(N(1,iBin),2)== 1 && N(1,iBin) ~= 1
                %                 XtmpLocation = [-(N(1,iBin))/2:(N(1,iBin))/2];
                %             end
                sorted_data_x = [sorted_data_x XtmpLocation];
            end
        end
        sorted_data_x = sorted_data_x/10+barIndex;
        line(sorted_data_x, sorted_data_y(:,1), 'linestyle', 'none', 'marker', '.','MarkerSize', 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [0 0 0])

        StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
        tmp(barIndex) =  mean(Data);
    end
    e = errorbar([1:2],tmp,StandardError,'.');
    e.MarkerSize = 1;
    e.Color = 'black';
    e.CapSize = 10;
    clear tmp


    %adding significance
    [h,p,ci,stats] = ttest2(DataExp1(:,iPlot),DataExp2(:,iPlot));
    HighlySig = find(p<0.01);
    p(p<0.01)=1;
    Sig = find(p<0.05);
    tmp = ylim;
    %plot sig
    if sum(HighlySig,2)>0
        plot(1.4,tmp(end)+0.5,'k*','MarkerSize',5)
        plot(1.6,tmp(end)+0.5,'k*','MarkerSize',5)
        plot([1,2],[tmp(end) tmp(end)],'k')
    end
    if sum(Sig,2)>0
        plot(1.5,tmp(end)+0.5,'k*','MarkerSize',5)
        plot([1,2],[tmp(end) tmp(end)],'k')

    end
    title(WBDataRawExp1.Properties.VariableNames{iPlot+1}(end-4:end))

end
VectorExperiment = [ones(size(DataExp1,1)*4,1);ones(size(DataExp2,1)*4,1)*2];
VectorProtein = [ones(size(DataExp1,1),1);ones(size(DataExp1,1),1)*2;ones(size(DataExp1,1),1)*3;ones(size(DataExp1,1),1)*4;...
    ones(size(DataExp2,1),1);ones(size(DataExp2,1),1)*2;ones(size(DataExp2,1),1)*3;ones(size(DataExp2,1),1)*4];
ANOVAData = [DataExp1(:,1);DataExp1(:,2);DataExp1(:,3);DataExp1(:,4);...
    DataExp2(:,1);DataExp2(:,2);DataExp2(:,3);DataExp2(:,4)];
[p,tablestats,stats,terms] = anovan(ANOVAData',{VectorExperiment',VectorProtein'}, 'model','interaction', 'display','off');
txt = {'Hypothalamus:', strcat('Experiment: ',num2str(round(tablestats{2, 7},3)),' F: ',num2str(round(tablestats{2, 6},3))),...
    strcat('Protein',num2str(round(tablestats{3, 7},3)),' F: ',num2str(round(tablestats{3, 6},3))),...
    strcat('Interact',num2str(round(tablestats{4, 7},3)),' F: ',num2str(round(tablestats{4, 6},3)))}
DataExp1Hypo = DataExp1;
DataExp2Hypo = DataExp2;

% DataExp1 = WBDataExp1.Cortex(1:16,:)./repmat(mean(WBDataExp1.Cortex(17:32,:),1),[16,1]);
DataExp1 = WBDataExp1.Cortex(1:10,:)./repmat(mean(WBDataExp1.Cortex(11:20,:),1),[10,1]);

DataExp2 = WBDataExp2.Cortex(1:8,:)./repmat(mean(WBDataExp2.Cortex(9:16,:),1),[8,1]);

hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 20 500 500])
set(gcf,'color','white')
set(gcf,'name',' Protein levels Exp1 vs Exp2 Cortex','NumberTitle','off')
tiledlayout(1,4)
for iPlot = 1:4
    nexttile

    b = bar([1 2], [mean(DataExp1(:,iPlot),1), mean(DataExp2(:,iPlot),1)]);
    b.FaceColor = 'flat';
    b.CData([1],:) = [1 0 0];
    b.CData([2],:) = [0 0 1];

    ylabel('Sleep/Wake or TSD');
    Labels = {'Exp1','Exp2'};
    set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
    set(gca,'FontSize',8,'TickLength',[0.025 0.025])
    set(gca,'TickDir','out');
    set(gca, 'box', 'off')
    hold on
    StandardError = [];
    tmp =  [];
    %adding data points for distribution
    for barIndex = 1: size(b.YData,2)
        if barIndex ==1
            Data = DataExp1(:,iPlot);
        end
        if barIndex ==2
            Data = DataExp2(:,iPlot);
        end

        sorted_data_y = sort(Data, 'ascend');


        sorted_data_x = [];
        [N,edges] = histcounts(Data(:,1),100);
        for iBin = 1: size(N,2)
            if N(1,iBin) >0
                %             if mod(N(1,iBin),2)== 0 %gerade Zahlen
                XtmpLocation = [-(N(1,iBin)-1)/2:(N(1,iBin)-1)/2];
                %             end
                if N(1,iBin) == 1
                    XtmpLocation = 0;
                end
                %             if mod(N(1,iBin),2)== 1 && N(1,iBin) ~= 1
                %                 XtmpLocation = [-(N(1,iBin))/2:(N(1,iBin))/2];
                %             end
                sorted_data_x = [sorted_data_x XtmpLocation];
            end
        end
        sorted_data_x = sorted_data_x/10+barIndex;
        line(sorted_data_x, sorted_data_y(:,1), 'linestyle', 'none', 'marker', '.','MarkerSize', 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [0 0 0])

        StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
        tmp(barIndex) =  mean(Data);
    end
    e = errorbar([1:2],tmp,StandardError,'.');
    e.MarkerSize = 1;
    e.Color = 'black';
    e.CapSize = 10;
    clear tmp


    %adding significance
    [h,p,ci,stats] = ttest2(DataExp1(:,iPlot),DataExp2(:,iPlot));
    HighlySig = find(p<0.01);
    p(p<0.01)=1;
    Sig = find(p<0.05);
    tmp = ylim;
    %plot sig
    if sum(HighlySig,2)>0
        plot(1.4,tmp(end)+0.5,'k*','MarkerSize',5)
        plot(1.6,tmp(end)+0.5,'k*','MarkerSize',5)
        plot([1,2],[tmp(end) tmp(end)],'k')
    end
    if sum(Sig,2)>0
        plot(1.5,tmp(end)+0.5,'k*','MarkerSize',5)
        plot([1,2],[tmp(end) tmp(end)],'k')

    end
    title(WBDataRawExp1.Properties.VariableNames{iPlot+1}(end-4:end))

end

VectorExperiment = [ones(size(DataExp1,1)*4,1);ones(size(DataExp2,1)*4,1)*2];
VectorProtein = [ones(size(DataExp1,1),1);ones(size(DataExp1,1),1)*2;ones(size(DataExp1,1),1)*3;ones(size(DataExp1,1),1)*4;...
    ones(size(DataExp2,1),1);ones(size(DataExp2,1),1)*2;ones(size(DataExp2,1),1)*3;ones(size(DataExp2,1),1)*4];
ANOVAData = [DataExp1(:,1);DataExp1(:,2);DataExp1(:,3);DataExp1(:,4);...
    DataExp2(:,1);DataExp2(:,2);DataExp2(:,3);DataExp2(:,4)];
[p,tablestats,stats,terms] = anovan(ANOVAData,{VectorExperiment',VectorProtein'}, 'model','interaction', 'display','off');
txt = {'Cortex:', strcat('Experiment: ',num2str(round(tablestats{2, 7},3)),' F: ',num2str(round(tablestats{2, 6},3))),...
    strcat('Protein',num2str(round(tablestats{3, 7},3)),' F: ',num2str(round(tablestats{3, 6},3))),...
    strcat('Interact',num2str(round(tablestats{4, 7},3)),' F: ',num2str(round(tablestats{4, 6},3)))}


DataExp1Cortex = DataExp1;
DataExp2Cortex = DataExp2;
%repeated measures ANOVA FOR Exp1/Exp2 * Hypo/Cortex
for iProtein= 1:4
%     meas = [DataExp1Hypo(:,iProtein) DataExp2Hypo(:,iProtein)];
%     meas = [meas;[DataExp1Cortex(:,iProtein) DataExp2Cortex(:,iProtein)]];

    meas = [DataExp1Hypo(:,iProtein) DataExp1Cortex(:,iProtein)];
    meas = [meas;[DataExp2Hypo(:,iProtein) DataExp2Cortex(:,iProtein)]];

    clear Stage
    for iCells = 1: size(DataExp1Hypo(:,iProtein),1)
        Stage{iCells,1} = 'Exp1';
    end
    for iCells = size(DataExp1Hypo(:,iProtein))+1:size(DataExp1Hypo(:,iProtein),1)+size(DataExp2Hypo(:,iProtein),1)
        Stage{iCells,1} = 'Exp2';
    end
    t = table(Stage,meas(:,1),meas(:,2),...
        'VariableNames',{'Experiment','Hypo','Cortex'});
    Meas = table([{'Hypo'} {'Cortex'}]','VariableNames',{'Area'});
    rm = fitrm(t,'Hypo-Cortex~Experiment','WithinDesign',Meas);
    ranovatbl = ranova(rm,'WithinModel','Area');
    FValue = ranovatbl{5,4}
    PValue = ranovatbl{5,5}

%   Regions  = repmat([repmat({'cort'}, 1,1); repmat({'hypo'}, 1,1)],1,1);
%   Proteins = repmat({'p1'}, 2,1); %Proteins = repmat({'p1'; 'p2'; 'p3'; 'p4'}, 2,1);
%   within = table(Regions, Proteins, 'VariableNames', {'Regions', 'Proteins'});
%   rm = fitrm(t,'V1-V2 ~ Group','WithinDesign',within); %, 'WithinModel','orthogonalcontrasts')
%   output.exp_SleepWake.rm = rm;
%   [ranovatbl] = ranova(rm, 'WithinModel','Regions');
%   output.exp_SleepWake.ranovatbl = ranovatbl;

%   fp{iProtein,1} = protein_str{iProtein};
  fp{iProtein,2} = ['F_interaction(',num2str(ranovatbl.DF(5), '%.0f'),', ', num2str(ranovatbl.DF(5+1),'%.0f'),') = ',  num2str(ranovatbl.F(5),'%.4f'),...
    ', p = ',  num2str(ranovatbl.pValue(5),'%.4f')];
  fp{iProtein,3} = ['F_group(',num2str(ranovatbl.DF(2),'%.0f'),', ',  num2str(ranovatbl.DF(2+1),'%.0f'),') = ',  num2str(ranovatbl.F(2),'%.4f'),...
    ', p = ',  num2str(ranovatbl.pValue(2),'%.4f')];
  fp{iProtein,4} = ['F_region(',num2str(ranovatbl.DF(4),'%.0f'),', ',  num2str(ranovatbl.DF(4+2),'%.0f'),') = ',  num2str(ranovatbl.F(4),'%.4f'),...
    ', p = ',  num2str(ranovatbl.pValue(4),'%.4f')];

end



%%
WBDataRawExp1 = readtable('Y:\niels\Jianfeng\WesternBlotResults\Exp1_WB_raw_26wells12wells-1_DataManfred.csv');
WBAnimalNamesExp1 = table2array(WBDataRawExp1(:,1));
WBDataExp1.Cortex = table2array(WBDataRawExp1(:,2:5));
WBDataExp1.HypoThalamus = table2array(WBDataRawExp1(:,8:11));
WBMeasureNamesTmp = fieldnames(WBDataRawExp1);
WBMeasureNamesExp1.Cortex = WBMeasureNamesTmp(2:5,1);
WBMeasureNamesExp1.Hypo = WBMeasureNamesTmp(8:11,1);

% DataExp1 = WBDataExp1.HypoThalamus(1:16,:)./repmat(mean(WBDataExp1.HypoThalamus(17:32,:),1),[16,1]);
%repeated measures ANOVA FOR Sleep/WAK * Hypo/Cortex
for iProtein= 1:4
%     meas = [WBDataExp1.HypoThalamus(1:10,iProtein) WBDataExp1.HypoThalamus(11:20,iProtein)];
%     meas = [meas;[WBDataExp1.Cortex(1:10,iProtein) WBDataExp1.Cortex(11:20,iProtein)]];

    meas = [WBDataExp1.HypoThalamus(1:10,iProtein) WBDataExp1.Cortex(1:10,iProtein)];
    meas = [meas;[WBDataExp1.HypoThalamus(11:20,iProtein) WBDataExp1.Cortex(11:20,iProtein)]];

    clear Stage
    for iCells = 1: size(WBDataExp1.HypoThalamus(1:10,iProtein),1)
        Stage{iCells,1} = 'Sleep';
    end
    for iCells = size(WBDataExp1.Cortex(1:10,iProtein))+1:size(WBDataExp1.HypoThalamus(1:10,iProtein),1)+size(WBDataExp1.Cortex(1:10,iProtein),1)
        Stage{iCells,1} = 'Wake';
    end
    t = table(Stage,meas(:,1),meas(:,2),...
        'VariableNames',{'Stage','Hypo','Cortex'});
    Meas = table([{'Hypo'} {'Cortex'}]','VariableNames',{'Area'});
    rm = fitrm(t,'Hypo-Cortex~Stage','WithinDesign',Meas);
    ranovatbl = ranova(rm,'WithinModel','Area');
    FValue = ranovatbl{5,4}
    PValue = ranovatbl{5,5}

%   Regions  = repmat([repmat({'cort'}, 1,1); repmat({'hypo'}, 1,1)],1,1);
%   Proteins = repmat({'p1'}, 2,1); %Proteins = repmat({'p1'; 'p2'; 'p3'; 'p4'}, 2,1);
%   within = table(Regions, Proteins, 'VariableNames', {'Regions', 'Proteins'});
%   rm = fitrm(t,'V1-V2 ~ Group','WithinDesign',within); %, 'WithinModel','orthogonalcontrasts')
%   output.exp_SleepWake.rm = rm;
%   [ranovatbl] = ranova(rm, 'WithinModel','Regions');
%   output.exp_SleepWake.ranovatbl = ranovatbl;

%   fp{iProtein,1} = protein_str{iProtein};
  fp{iProtein,2} = ['F_interaction(',num2str(ranovatbl.DF(5), '%.0f'),', ', num2str(ranovatbl.DF(5+1),'%.0f'),') = ',  num2str(ranovatbl.F(5),'%.4f'),...
    ', p = ',  num2str(ranovatbl.pValue(5),'%.4f')];
  fp{iProtein,3} = ['F_group(',num2str(ranovatbl.DF(2),'%.0f'),', ',  num2str(ranovatbl.DF(2+1),'%.0f'),') = ',  num2str(ranovatbl.F(2),'%.4f'),...
    ', p = ',  num2str(ranovatbl.pValue(2),'%.4f')];
  fp{iProtein,4} = ['F_region(',num2str(ranovatbl.DF(4),'%.0f'),', ',  num2str(ranovatbl.DF(4+2),'%.0f'),') = ',  num2str(ranovatbl.F(4),'%.4f'),...
    ', p = ',  num2str(ranovatbl.pValue(4),'%.4f')];

end





%%
%Experiment 1
DirData = 'Y:\niels\Jianfeng\Results\PerHour\';
DirSave = 'Y:\niels\Jianfeng\Results\PerHourAllAnimals\';
WBDataRaw = readtable('Y:\niels\Jianfeng\WesternBlotResults\WB_raw_26wells12wells-1.csv');
WBAnimalNames = table2array(WBDataRaw(:,1));
WBData.HypoThalamus = table2array(WBDataRaw(:,2:5));
WBData.Cortex = table2array(WBDataRaw(:,8:11));
WBMeasureNamesTmp = fieldnames(WBDataRaw);
WBMeasureNames.Hypo = WBMeasureNamesTmp(2:5,1);
WBMeasureNames.Cortex = WBMeasureNamesTmp(8:11,1);

% DirData = '/gpfs01/born/animal/niels/Jianfeng/EEGDataScored/Matlab/';
% DirSave = '/gpfs01/born/animal/niels/Jianfeng/Results/';
Files = dir(DirData);
clear AnimalName
for iFile = 3: size(Files,1)
    AnimalName{iFile-2,1} = Files(iFile,1).name;
end
iWake =1;
iSleep = 1;
for iAnimal = 1:size(AnimalName,1)
    load(strcat(DirData,AnimalName{iAnimal,1}));
    ResultNames = fieldnames(ResultsPerHour{1});
    
    for iCh = 1:3
        clear CurrData
        for iGraph = 1: numel(fieldnames(ResultsPerHour{1}))-6 %last six variables are alredy mean or sum
            for ihour = 1:6
                if ~isempty(ResultsPerHour{ihour}.(ResultNames{iGraph})) && size(ResultsPerHour{ihour}.(ResultNames{iGraph}),1)>1
                    CurrData(1,ihour) = nanmean(ResultsPerHour{ihour}.(ResultNames{iGraph}){iCh,1});
                else
                    CurrData(1,ihour) = 0;
                end
            end
            
            if isempty(strfind(AnimalName{iAnimal,1},'Sleep')) %check if animal is from wake group or sleep group
                Exp1.ResultsAllWAKAnimals.(ResultNames{iGraph}){iCh,1}(iWake,:) = CurrData;
                
            else
                Exp1.ResultsAllSLEEPAnimals.(ResultNames{iGraph}){iCh,1}(iSleep,:) = CurrData;
                
            end
            
            
        end
        clear CurrData
        for iGraph = numel(fieldnames(ResultsPerHour{1}))-5 : numel(fieldnames(ResultsPerHour{1}))
            for ihour = 1:6
                CurrData(1,ihour) = nanmean(ResultsPerHour{ihour}.(ResultNames{iGraph}));
            end
            if isempty(strfind(AnimalName{iAnimal,1},'Sleep')) %check if animals is from wake group or sleep group
                Exp1.ResultsAllWAKAnimals.(ResultNames{iGraph}){iCh,1}(iWake,:) = CurrData;

            else
                Exp1.ResultsAllSLEEPAnimals.(ResultNames{iGraph}){iCh,1}(iSleep,:) = CurrData;

            end
        end
    end
    
    if isempty(strfind(AnimalName{iAnimal,1},'Sleep')) %check if animals is from wake group or sleep group
        Exp1.ResultsAllWAKAnimals.WBResults.Hypo(iWake,:) =  WBData.HypoThalamus(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);        
        Exp1.ResultsAllWAKAnimals.WBResults.Cortex(iWake,:) =  WBData.Cortex(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);        
    else
        Exp1.ResultsAllSLEEPAnimals.WBResults.Hypo(iSleep,:) = WBData.HypoThalamus(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);        
        Exp1.ResultsAllSLEEPAnimals.WBResults.Cortex(iSleep,:) =  WBData.Cortex(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);        
    end
    
    
    if isempty(strfind(AnimalName{iAnimal,1},'Sleep'))
        iWake = iWake + 1;
    else
        iSleep = iSleep + 1;
    end
    
end



%Experiment 2
% restoredefaultpath
addpath 'Y:\niels\software'
DirData = 'Y:\niels\LunAMPAScaling\NewDataSet\Results\PerHour\';
DirSave = 'Y:\niels\LunAMPAScaling\NewDataSet\Results\PerHourAllAnimals\';
WBDataRaw = readtable('Y:\niels\LunAMPAScaling\NewDataSet\RSD_WB_Data_Cort_Hypo_Merg.csv');
WBAnimalNames = table2array(WBDataRaw(:,1));
WBData.HypoThalamus = table2array(WBDataRaw(:,2:5));
WBData.Cortex = table2array(WBDataRaw(:,8:11));
WBMeasureNamesTmp = fieldnames(WBDataRaw);
WBMeasureNames.Hypo = WBMeasureNamesTmp(2:5,1);
WBMeasureNames.Cortex = WBMeasureNamesTmp(8:11,1);

Files = dir(DirData);
clear AnimalName
for iFile = 4: size(Files,1)
    AnimalName{iFile-3,1} = Files(iFile,1).name;
end
iTSD =1;
iSleep = 1;
iRSD =1;

for iAnimal = 1:size(AnimalName,1)
    load(strcat(DirData,AnimalName{iAnimal,1}));
    ResultNames = fieldnames(ResultsPerHour{1});

    for iCh = 1:3
        clear CurrData
        for iGraph = 1: numel(fieldnames(ResultsPerHour{1}))-6 %last six variables are alredy mean or sum
            for ihour = 1:6
                if ~isempty(ResultsPerHour{ihour}.(ResultNames{iGraph})) && size(ResultsPerHour{ihour}.(ResultNames{iGraph}),1)>1
                    CurrData(1,ihour) = nanmean(ResultsPerHour{ihour}.(ResultNames{iGraph}){iCh,1});
                else
                    CurrData(1,ihour) = 0;
                end
            end

            if ~isempty(strfind(AnimalName{iAnimal,1},'TSD')) %check if animal is from TSD group
                Exp2.ResultsAllTSDAnimals.(ResultNames{iGraph}){iCh,1}(iTSD,:) = CurrData;
            end
            if ~isempty(strfind(AnimalName{iAnimal,1},'RSD')) %check if animal is from RSD group
                Exp2.ResultsAllRSDAnimals.(ResultNames{iGraph}){iCh,1}(iRSD,:) = CurrData;
            end
            if ~isempty(strfind(AnimalName{iAnimal,1},'sleep')) %check if animal is from sleep group
                Exp2.ResultsAllSLEEPAnimals.(ResultNames{iGraph}){iCh,1}(iSleep,:) = CurrData;

            end


        end
        clear CurrData
        for iGraph = numel(fieldnames(ResultsPerHour{1}))-5 : numel(fieldnames(ResultsPerHour{1}))
            for ihour = 1:6
                CurrData(1,ihour) = nanmean(ResultsPerHour{ihour}.(ResultNames{iGraph}));
            end
            if ~isempty(strfind(AnimalName{iAnimal,1},'TSD')) %check if animal is from TSD group
                Exp2.ResultsAllTSDAnimals.(ResultNames{iGraph}){iCh,1}(iTSD,:) = CurrData;
            end
            if ~isempty(strfind(AnimalName{iAnimal,1},'RSD')) %check if animal is from RSD group
                Exp2.ResultsAllRSDAnimals.(ResultNames{iGraph}){iCh,1}(iRSD,:) = CurrData;
            end
            if ~isempty(strfind(AnimalName{iAnimal,1},'sleep')) %check if animal is from sleep group
                Exp2.ResultsAllSLEEPAnimals.(ResultNames{iGraph}){iCh,1}(iSleep,:) = CurrData;

            end
        end
    end

    if ~isempty(strfind(AnimalName{iAnimal,1},'TSD')) %check if animals is from wake group or sleep group
        Exp2.ResultsAllTSDAnimals.WBResults.Hypo(iTSD,:) =  WBData.HypoThalamus(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);
        Exp2.ResultsAllTSDAnimals.WBResults.Cortex(iTSD,:) =  WBData.Cortex(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);
    end
    if ~isempty(strfind(AnimalName{iAnimal,1},'RSD')) %check if animals is from wake group or sleep group
        Exp2.ResultsAllRSDAnimals.WBResults.Hypo(iRSD,:) =  WBData.HypoThalamus(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);
        Exp2.ResultsAllRSDAnimals.WBResults.Cortex(iRSD,:) =  WBData.Cortex(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);
    end
    if ~isempty(strfind(AnimalName{iAnimal,1},'sleep'))
        Exp2.ResultsAllSLEEPAnimals.WBResults.Hypo(iSleep,:) = WBData.HypoThalamus(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);
        Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(iSleep,:) =  WBData.Cortex(find(strncmp(WBAnimalNames,AnimalName{iAnimal,1}(15:end),6)),:);
    end


    if ~isempty(strfind(AnimalName{iAnimal,1},'TSD')) %check if animal is from TSD group
        iTSD = iTSD +1;
    end
    if ~isempty(strfind(AnimalName{iAnimal,1},'RSD')) %check if animal is from RSD group
        iRSD =iRSD + 1;
    end
    if ~isempty(strfind(AnimalName{iAnimal,1},'sleep')) %check if animal is from sleep group
        iSleep = iSleep + 1;
    end

end


%%
%prepare figures


hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [(iCh-1)*500 50 500 900])
set(gcf,'color','white')
set(gca,'box','off')
tiledlayout(4,2)
nexttile
b = bar([1 2 3], [mean(sum(Exp1.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,4:6),2),1)/60, mean(sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,4:6),2),1)/60, mean(sum(Exp2.ResultsAllRSDAnimals.AbsNREMDur{1,1}(:,4:6),2),1)/60]);
b.FaceColor = 'flat';
b.CData([1],:) = [1 1 1];
b.CData([2],:) = [1 1 1];
b.CData([3],:) = [1 1 1];
b.FaceAlpha = 1;


ylabel('SWS Duration last 3 hours in min');
Labels = {'Exp1 Sleep','Exp2 Sleep','Exp2 RSD'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on
%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    if barIndex ==1
        Data = sum(Exp1.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,4:6),2)/60;
    end
    if barIndex ==2
        Data = sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,4:6),2)/60;
    end
    if barIndex ==3
        Data = sum(Exp2.ResultsAllRSDAnimals.AbsNREMDur{1,1}(:,4:6),2)/60;
    end
    
    sorted_data_y = sort(Data, 'ascend');


    sorted_data_x = [];
    [N,edges] = histcounts(Data(:,1),100);
    for iBin = 1: size(N,2)
        if N(1,iBin) >0
            %             if mod(N(1,iBin),2)== 0 %gerade Zahlen
            XtmpLocation = [-(N(1,iBin)-1)/2:(N(1,iBin)-1)/2];
            %             end
            if N(1,iBin) == 1
                XtmpLocation = 0;
            end
            %             if mod(N(1,iBin),2)== 1 && N(1,iBin) ~= 1
            %                 XtmpLocation = [-(N(1,iBin))/2:(N(1,iBin))/2];
            %             end
            sorted_data_x = [sorted_data_x XtmpLocation];
        end
    end
    sorted_data_x = sorted_data_x/10+barIndex;
    line(sorted_data_x, sorted_data_y(:,1), 'linestyle', 'none', 'marker', '.','MarkerSize', 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [0 0 0])

    StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
    tmp(barIndex) =  mean(Data);
end
e = errorbar([1:3],tmp,StandardError,'.');
e.MarkerSize = 1;
e.Color = 'black';
e.CapSize = 10;
clear tmp
title('Variability SWS Duration last 3 hours');

nexttile

ylabel('SWS in min');
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on
plot([1:6],mean(Exp1.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 1 1])
hold all
plot([1:6],mean(Exp1.ResultsAllWAKAnimals.AbsNREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 0 0])
xticks([1 : 6]);
ylim([0 50])
xlim([1 6])
xticklabels ({'1h', '2h', '3h', '4h', '5h', '6h'});
for iGroup = 1:2
    if iGroup ==1
        b = bar([1:6], mean(Exp1.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    else
        b = bar([1:6], mean(Exp1.ResultsAllWAKAnimals.AbsNREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    end
b.FaceColor = 'none';
%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    if iGroup ==1
        Data = (Exp1.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,barIndex))/60;
    end
    if iGroup ==2
        Data = (Exp1.ResultsAllWAKAnimals.AbsNREMDur{1,1}(:,barIndex))/60;
    end
   
    StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
    tmp(barIndex) =  mean(Data);
end
e = errorbar([1:6],tmp,StandardError,'.');
e.MarkerSize = 1;
e.Color = 'black';
e.CapSize = 10;
clear tmp
end
title('SWS amount Exp1');

nexttile

ylabel('SWS in min');
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on
plot([1:6],mean(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 1 1])
hold all
plot([1:6],mean(Exp2.ResultsAllRSDAnimals.AbsNREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 0.5])
xticks([1 : 6]);
ylim([0 50])
xlim([1 6])
xticklabels ({'1h', '2h', '3h', '4h', '5h', '6h'});
for iGroup = 1:2
    if iGroup ==1
        b = bar([1:6], mean(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    else
        b = bar([1:6], mean(Exp2.ResultsAllRSDAnimals.AbsNREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    end
b.FaceColor = 'none';
%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    if iGroup ==1
        Data = (Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,barIndex))/60;
    end
    if iGroup ==2
        Data = (Exp2.ResultsAllRSDAnimals.AbsNREMDur{1,1}(:,barIndex))/60;
    end
   
    StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
    tmp(barIndex) =  mean(Data);
end
e = errorbar([1:6],tmp,StandardError,'.');
e.MarkerSize = 1;
e.Color = 'black';
e.CapSize = 10;
clear tmp
end
title('SWS amount Exp2');

nexttile
b = bar([1 2 3], [mean(sum(Exp1.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,4:6),2),1)/60, mean(sum(Exp2.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,4:6),2),1)/60, mean(sum(Exp2.ResultsAllRSDAnimals.AbsREMDur{1,1}(:,4:6),2),1)/60]);
b.FaceColor = 'flat';
b.CData([1],:) = [1 1 1];
b.CData([2],:) = [1 1 1];
b.CData([3],:) = [1 1 1];
b.FaceAlpha = 1;

ylabel('REM Duration last 3 hours in min');
Labels = {'Exp1 Sleep','Exp2 Sleep','Exp2 RSD'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on
StandardError = [];
tmp =  [];
%adding data points for distribution
for barIndex = 1: size(b.YData,2)
    if barIndex ==1
        Data = sum(Exp1.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,4:6),2)/60;
    end
    if barIndex ==2
        Data = sum(Exp2.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,4:6),2)/60;
    end
    if barIndex ==3
        Data = sum(Exp2.ResultsAllRSDAnimals.AbsREMDur{1,1}(:,4:6),2)/60;
    end

    sorted_data_y = sort(Data, 'ascend');


    sorted_data_x = [];
    [N,edges] = histcounts(Data(:,1),100);
    for iBin = 1: size(N,2)
        if N(1,iBin) >0
            %             if mod(N(1,iBin),2)== 0 %gerade Zahlen
            XtmpLocation = [-(N(1,iBin)-1)/2:(N(1,iBin)-1)/2];
            %             end
            if N(1,iBin) == 1
                XtmpLocation = 0;
            end
            %             if mod(N(1,iBin),2)== 1 && N(1,iBin) ~= 1
            %                 XtmpLocation = [-(N(1,iBin))/2:(N(1,iBin))/2];
            %             end
            sorted_data_x = [sorted_data_x XtmpLocation];
        end
    end
    sorted_data_x = sorted_data_x/10+barIndex;
    line(sorted_data_x, sorted_data_y(:,1), 'linestyle', 'none', 'marker', '.','MarkerSize', 15, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [0 0 0])

    StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
    tmp(barIndex) =  mean(Data);
end
e = errorbar([1:3],tmp,StandardError,'.');
e.MarkerSize = 1;
e.Color = 'black';
e.CapSize = 10;
clear tmp
title('Variability REM Duration last 3 hours');

nexttile


ylabel('REM in min');
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on
plot([1:6],mean(Exp1.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 1 1])
hold all
plot([1:6],mean(Exp1.ResultsAllWAKAnimals.AbsREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 0 0])
xticks([1 : 6]);
ylim([0 6])
xlim([1 6])
xticklabels ({'1h', '2h', '3h', '4h', '5h', '6h'});
for iGroup = 1:2
    if iGroup ==1
        b = bar([1:6], mean(Exp1.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    else
        b = bar([1:6], mean(Exp1.ResultsAllWAKAnimals.AbsREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    end
    b.FaceColor = 'none';
    %adding data points for distribution
    for barIndex = 1: size(b.YData,2)
        if iGroup ==1
            Data = (Exp1.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,barIndex))/60;
        end
        if iGroup ==2
            Data = (Exp1.ResultsAllWAKAnimals.AbsREMDur{1,1}(:,barIndex))/60;
        end

        StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
        tmp(barIndex) =  mean(Data);
    end
    e = errorbar([1:6],tmp,StandardError,'.');
    e.MarkerSize = 1;
    e.Color = 'black';
    e.CapSize = 10;
    clear tmp
end
title('REM amount Exp1');

nexttile

ylabel('REM in min');
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on
plot([1:6],mean(Exp2.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 1 1])
hold all
plot([1:6],mean(Exp2.ResultsAllRSDAnimals.AbsREMDur{1,1}(:,:)/60),'-o','Color','k','LineWidth',1,...
    'MarkerSize',5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 0.5])
xticks([1 : 6]);
ylim([0 6])
xlim([1 6])
xticklabels ({'1h', '2h', '3h', '4h', '5h', '6h'});
for iGroup = 1:2
    if iGroup ==1
        b = bar([1:6], mean(Exp2.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    else
        b = bar([1:6], mean(Exp2.ResultsAllRSDAnimals.AbsREMDur{1,1}(:,:)/60,1),'EdgeColor','none');
    end
    b.FaceColor = 'none';
    %adding data points for distribution
    for barIndex = 1: size(b.YData,2)
        if iGroup ==1
            Data = (Exp2.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,barIndex))/60;
        end
        if iGroup ==2
            Data = (Exp2.ResultsAllRSDAnimals.AbsREMDur{1,1}(:,barIndex))/60;
        end

        StandardError(barIndex) = std(Data)./repmat(sqrt(length(Data)),1,1);
        tmp(barIndex) =  mean(Data);
    end
    e = errorbar([1:6],tmp,StandardError,'.');
    e.MarkerSize = 1;
    e.Color = 'black';
    e.CapSize = 10;
    clear tmp
end
title('REM amount Exp2');

nexttile


b = bar([1 2 3], [var(sum(Exp1.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,4:6),2),1), var(sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1,1}(:,4:6),2),1), var(sum(Exp2.ResultsAllRSDAnimals.AbsNREMDur{1,1}(:,4:6),2),1)]);
b.FaceColor = 'flat';
b.CData([1],:) = [1 1 1];
b.CData([2],:) = [1 1 1];
b.CData([3],:) = [1 1 1];
b.FaceAlpha = 1;

ylabel('Variance SWS Duration last 3 hours');
Labels = {'Exp1 Sleep','Exp2 Sleep','Exp2 RSD'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on

title('Variability NREM Duration last 3 hours');

nexttile
b = bar([1 2 3], [var(sum(Exp1.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,4:6),2),1), var(sum(Exp2.ResultsAllSLEEPAnimals.AbsREMDur{1,1}(:,4:6),2),1), var(sum(Exp2.ResultsAllRSDAnimals.AbsREMDur{1,1}(:,4:6),2),1)]);
b.FaceColor = 'flat';
b.CData([1],:) = [1 1 1];
b.CData([2],:) = [1 1 1];
b.CData([3],:) = [1 1 1];
b.FaceAlpha = 1;

ylabel('Variance SWS Duration last 3 hours');
Labels = {'Exp1 Sleep','Exp2 Sleep','Exp2 RSD'};
set(gca, 'XTick', 1:3, 'XTickLabel', Labels);
set(gca,'FontSize',8,'TickLength',[0.025 0.025])
set(gca,'TickDir','out');
set(gca, 'box', 'off')
hold on

title('Variability NREM Duration last 3 hours');

hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [1000 50 800 400])
set(gcf,'color','white')
set(gca,'box','off')
iCh = 1;
tiledlayout(1,2)
CurrentResultsEphys = [sum(Exp2.ResultsAllSLEEPAnimals.TotalNumberOfSpi{iCh}(:,4:6),2)./sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{iCh}(:,4:6),2);...
    sum(Exp2.ResultsAllRSDAnimals.TotalNumberOfSpi{iCh}(:,4:6),2)./sum(Exp2.ResultsAllRSDAnimals.AbsNREMDur{iCh}(:,4:6),2)];
CurrentResultsWB = [Exp2.ResultsAllSLEEPAnimals.WBResults.Hypo(:,1); Exp2.ResultsAllRSDAnimals.WBResults.Hypo(:,1)];
nexttile
scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1); repmat([0 1 0],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('Spi Density');
ylabel('GluA1 Hypo');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');

nexttile
CurrentResultsEphys = [sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{iCh}(:,4:6),2);...
    sum(Exp2.ResultsAllRSDAnimals.AbsNREMDur{iCh}(:,4:6),2)];
CurrentResultsWB = [Exp2.ResultsAllSLEEPAnimals.WBResults.Hypo(:,1); Exp2.ResultsAllRSDAnimals.WBResults.Hypo(:,1)];

scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1); repmat([0 1 0],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('SWS Duration');
ylabel('GluA1 Hypo');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');
title('correlations for Sleep and RSD animals')

%%
% Hypo and Cortex pooled
hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [500 50 1200 400])
set(gcf,'color','white')
set(gca,'box','off')
iCh = 1;
tiledlayout(1,3)
CurrentResultsEphys = sum(Exp2.ResultsAllSLEEPAnimals.TotalNumberOfSpi{iCh}(:,4:6),2)./sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{iCh}(:,4:6),2);
CurrentResultsWB = mean([Exp2.ResultsAllSLEEPAnimals.WBResults.Hypo(:,1), Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(:,1)],2);
nexttile
scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('Spi Density');
ylabel('Mean GluA1 Hypo + Cortex');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');

nexttile
CurrentResultsEphys = sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{iCh}(:,4:6),2);
CurrentResultsWB = mean([Exp2.ResultsAllSLEEPAnimals.WBResults.Hypo(:,1), Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(:,1)],2);

scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('SWS Duration');
ylabel('Mean GluA1 Hypo + Cortex');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');
title('correlations Merged Hypo and Cortex GluA1')
nexttile
CurrentResultsEphys = sum(Exp2.ResultsAllSLEEPAnimals.REMThetaEnergy{iCh}(:,4:6),2);
CurrentResultsWB = mean([Exp2.ResultsAllSLEEPAnimals.WBResults.Hypo(:,1), Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(:,1)],2);

scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('REM Theta Energy');
ylabel('Mean GluA1 Hypo + Cortex');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');

%%
% correlation spindle density across all channels for Cortex
hFig = figure;
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [500 50 1600 400])
set(gcf,'color','white')
set(gca,'box','off')
tiledlayout(1,4)
CurrentResultsEphys = mean([sum(Exp2.ResultsAllSLEEPAnimals.TotalNumberOfSpi{1}(:,4:6),2)./sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1}(:,4:6),2),...
    sum(Exp2.ResultsAllSLEEPAnimals.TotalNumberOfSpi{2}(:,4:6),2)./sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{2}(:,4:6),2),...
    sum(Exp2.ResultsAllSLEEPAnimals.TotalNumberOfSpi{3}(:,4:6),2)./sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{3}(:,4:6),2)],2);
CurrentResultsWB = mean([Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(:,1)],2);
nexttile
scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('Spi Density');
ylabel('GluA1 Cortex + Spi density across all Channels');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');
title('GluA1')

CurrentResultsWB = mean([Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(:,2)],2);
nexttile
scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('Spi Density');
ylabel('GluA2 Cortex + Spi density across all Channels');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');
title('GluA2')

CurrentResultsWB = mean([Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(:,3)],2);
nexttile
scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('Spi Density');
ylabel('Ser 845 Cortex + Spi density across all Channels');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');
title('845')

CurrentResultsWB = mean([Exp2.ResultsAllSLEEPAnimals.WBResults.Cortex(:,4)],2);
nexttile
scatter(CurrentResultsEphys,CurrentResultsWB,...
    [],[repmat([0 1 1],8,1)],'filled');

[r,p] = corr(CurrentResultsEphys,CurrentResultsWB);
xlabel('Spi Density');
ylabel('Ser 831 Cortex + Spi density across all Channels');
ylim([0 2.5])
xLimits = get(gca,'XLim');  % Get the range of the x axis
h2 = lsline;
text(xLimits(1,1),0.25,strcat('r=',num2str(r),' p =',num2str(p)));
set(h2,'color','black');
title('831')


%% stepwise regression
%Hypo GluA1
CurrentResultsEphys = [sum(Exp2.ResultsAllSLEEPAnimals.TotalNumberOfSpi{1}(:,4:6),2)./sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{1}(:,4:6),2),...
    sum(Exp2.ResultsAllSLEEPAnimals.AbsNREMDur{iCh}(:,4:6),2),...
    sum(Exp2.ResultsAllSLEEPAnimals.REMThetaEnergy{iCh}(:,4:6),2)];

CurrentResultsWB = Exp2.ResultsAllSLEEPAnimals.WBResults.Hypo(:,1);

mdl = stepwiselm(tmpData,'ResponseVar','AMPAreceptorDensity','PredictorVars',{'SpiDensity','NREMDuration','REMtheta'})

mdl = stepwiselm(tmpData,'linear','AMPAreceptorDensity ~ SpiDensity + NREMDuration + REMtheta','ResponseVar','AMPAreceptorDensity','PredictorVars',{'SpiDensity','NREMDuration','REMtheta'},'Criterion','AIC','Verbose',2)

mdl = stepwiselm(tmpData,'linear','ResponseVar','AMPAreceptorDensity','PredictorVars',{'SpiDensity','NREMDuration','REMtheta'},'Criterion','AIC','Verbose',2)
mdl = stepwiselm(tmpData,'interactions','ResponseVar','AMPAreceptorDensity','PredictorVars',{'SpiDensity','NREMDuration','REMtheta'},'Criterion','AIC','Verbose',2)
mdl = stepwiselm(tmpData,'constant','ResponseVar','AMPAreceptorDensity','PredictorVars',{'SpiDensity','NREMDuration','REMtheta'},'Criterion','AIC','Verbose',2)

mdl2 = fitlm([tmpData.SpiDensity tmpData.NREMDuration tmpData.REMtheta],tmpData.AMPAreceptorDensity)
step(mdl2)