clear
close all
% 
% DirData = 'Y:\niels\Jianfeng\EEGDataScored\Matlab\';
% DirSave = 'Y:\niels\Jianfeng\Results\DetectedEvents\';
DirData = '/gpfs01/born/animal/niels/Jianfeng/EEGDataScored/Matlab/';
DirSave = '/gpfs01/born/animal/niels/Jianfeng/Results/DetectedEventsSpiNew/';

% DirData = 'Y:\niels\LunAMPAScaling\MatlabData\';
% DirSave = 'Y:\niels\LunAMPAScaling\Results\DetectedEvents\';
% DirData = '/gpfs01/born/animal/niels/LunAMPAScaling/MatlabData/';
% DirSave = '/gpfs01/born/animal/niels/LunAMPAScaling/Results/DetectedEventsSpiNew/';
Files = dir(DirData);
for iFile = 3: size(Files,1)
    AnimalName{iFile-2,1} = Files(iFile,1).name;
end
numCh = 5;
ChannelNames = cell(numCh,1);
ChannelNames{1,1} = 'EEG1';
ChannelNames{2,1} = 'EEG2';
ChannelNames{3,1} = 'EEG3';
ChannelNames{4,1} = 'EMG1';
ChannelNames{5,1} = 'EMG2';
%Minimum Spindle Duration in ms
MinSpiDuration = 500;
MaxSpiDuration = 2500;
SpindleThreshold(1,1) = 1.00;
SpindleThreshold(2,1) = 1.5;
SpindleThreshold(3,1) = 2.0;
Fs = 1000;
%design filters
[SpiFilterHigh1,SpiFilterHigh2] = butter(3,2*11/Fs,'high');
[SpiFilterLow1,SpiFilterLow2] = butter(3,2*16/Fs,'low');
% SpiFilter = designfilt('bandpassiir','FilterOrder',20, ...
%     'HalfPowerFrequency1',9,'HalfPowerFrequency2',17, ...
%     'SampleRate',Fs);

[SOFilterHigh1,SOFilterHigh2] = butter(3,2*0.1/Fs,'high');
[SOFilterLow1,SOFilterLow2] = butter(3,2*4/Fs,'low');
% SOFilter = designfilt('bandpassiir','FilterOrder',20, ...
%     'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',4, ...
%     'SampleRate',Fs);

[ThetaFilterHigh1,ThetaFilterHigh2] = butter(3,2*5/Fs,'high');
[ThetaFilterLow1,ThetaFilterLow2] = butter(3,2*10/Fs,'low');
%%
for iAnimal = 6:size(AnimalName,1) 
    RawData = load(strcat(DirData,AnimalName{iAnimal,1}));
    RawChannelNames = fieldnames(RawData);
    RawScoring = RawData.SlStNew.codes(:,1);
    OrigScoringLength = length(RawScoring);

    
    %find NREMEpisodes
    NREMBegEpisode = strfind((RawScoring==2)',[0 1]);
    NREMEndEpisode = strfind((RawScoring==2)',[1 0]);
    if RawScoring(1,1) ==2
        NREMBegEpisode = [1 NREMBegEpisode];
    end
    if RawScoring(end,1) ==2
        NREMEndEpisode = [NREMEndEpisode length(RawScoring)];
    end
    NREMEpisodes = [NREMBegEpisode*10; NREMEndEpisode*10]; %create Matrix with NRem on and offset time in sec
    
    %find REMEpisodes
    REMBegEpisode = strfind((RawScoring==3)',[0 1]);
    REMEndEpisode = strfind((RawScoring==3)',[1 0]);
    if RawScoring(1,1) ==3
        REMBegEpisode = [1 REMBegEpisode];
    end
    if RawScoring(end,1) ==3
        REMEndEpisode = [REMEndEpisode length(RawScoring)];
    end
    REMEpisodes = [REMBegEpisode*10; REMEndEpisode*10]; %create Matrix with NRem on and offset time in sec
    
    %find WakeEpisodes
    WAKBegEpisode = strfind((RawScoring==1)',[0 1]);
    WAKEndEpisode = strfind((RawScoring==1)',[1 0]);
    if RawScoring(1,1) ==1
        WAKBegEpisode = [1 WAKBegEpisode];
    end
    if RawScoring(end,1) ==1
        WAKEndEpisode = [WAKEndEpisode length(RawScoring)];
    end
    WAKEpisodes = [WAKBegEpisode*10; WAKEndEpisode*10]; %create Matrix with NRem on and offset time in sec
    
    
    RawData = struct2cell(RawData);
    Data = zeros(RawData{1, 1}.length,numCh);
    for iCh = 1:numCh
        TempIdx = find(~cellfun('isempty', (strfind(RawChannelNames,ChannelNames{iCh,1})))); %find index of specific EEG/EMG Channel
        if RawData{TempIdx,1}.length < RawData{1, 1}.length %compensate that some channels might start one sample later
            Tmp = RawData{1, 1}.length - RawData{TempIdx,1}.length;
            Data(1+Tmp:end,iCh)= RawData{TempIdx,1}.values; 
        else
            Data(:,iCh)= RawData{TempIdx,1}.values; 
        end
    end
    %compenstae that scoring is shorter than Data
    DiffScoringData = OrigScoringLength*10*Fs - length(Data);
    Data(end-DiffScoringData:end,:) = [];
       
    %create upsampled scoring vector
    SleepScoring = zeros(size(Data,1),1);
    for iEpoch = 1:size(NREMEpisodes,2)
        SleepScoring(NREMEpisodes(1,iEpoch)*Fs:NREMEpisodes(2,iEpoch)*Fs) = 2;
    end
    %%
    %detect spindles
    SpiBand = zeros(size(Data,1),3);
    for iCh = 1:3
        SpiBand(:,iCh) = filtfilt(SpiFilterHigh1,SpiFilterHigh2,Data(:,iCh));
        SpiBand(:,iCh) = filtfilt(SpiFilterLow1,SpiFilterLow2,SpiBand(:,iCh));
    end
    SpiBandAmp = abs(hilbert(SpiBand));
    %calculate std for thresholds
    SpiAmpSTD = std(SpiBandAmp(SleepScoring==2,:));
    SpiAmpMean = mean(SpiBandAmp(SleepScoring==2,:));
        
    %Detect Spindles
    AllSpindles.FastSpi = cell(size(NREMEpisodes,2),3);
    for iEpoch = 1:size(NREMEpisodes,2)
        DataTmp = SpiBandAmp(NREMEpisodes(1,iEpoch)*Fs : NREMEpisodes(2,iEpoch)*Fs,:);
        for iCh = 1:3
            %%
            FastSpiAmplitudeTmp = smooth(DataTmp(:,iCh),40);%get smoothed instantaneous amplitude
            %detect timepoints above threshold
            above_threshold = FastSpiAmplitudeTmp > SpindleThreshold(1,1)*SpiAmpSTD(iCh);
            isLongEnough = bwareafilt(above_threshold, [MinSpiDuration, MaxSpiDuration]); %find spindle within duration range
            isLongEnough = [0; isLongEnough]; %compensate that spindle might start in the beginning
            SpiBeginning =  strfind(isLongEnough',[0 1]); %find spindle Beginning line before compensates that it find last 0
            SpiEnd = strfind(isLongEnough',[1 0])-1; %find spindle Ending subtract 1 because of added 0 in the beginning
            %check if Spinde was detected
            if ~isempty(SpiBeginning) || ~isempty(SpiEnd)
                %check if spindle started in the end of the Epoch
                if length(SpiEnd)<length(SpiBeginning)
                    SpiBeginning(:,end)=[];
                end
                if ~isempty(SpiBeginning) || ~isempty(SpiEnd) && SpiBeginning(1,1)==1
                    SpiBeginning(:,1) = [];
                    SpiEnd(:,1) = [];
                end
                FastSpindles = [SpiBeginning;SpiEnd];
                AllSpindles.FastSpi{iEpoch,iCh} = FastSpindles+(NREMEpisodes(1,iEpoch)*Fs);%include beginning of NREMEpoch
                
                
                
            else
                AllSpindles.FastSpi{iEpoch,iCh} = [];
            end
            if ~isempty(AllSpindles.FastSpi{iEpoch,iCh})
                AllSpindles.FastSpi{iEpoch,iCh}(:,AllSpindles.FastSpi{iEpoch,iCh}(2,:) >length(FastSpiAmplitudeTmp)-5000)=[];
            end
            
            %%
            %second threshold criteria
            CurrentSpindles = AllSpindles.FastSpi{iEpoch,iCh};
            TempIdx = [];
            for iSpi = 1: size (CurrentSpindles,2)
                
                DataTmpSpi = SpiBand(CurrentSpindles(1,iSpi)-5000:CurrentSpindles(2,iSpi)+5000,iCh); %get filteres spindle signal for eachspindle + - 5sec
                FastSpiAmplitudeTmp = smooth(abs(hilbert(DataTmpSpi)),40);%get smoothed instantaneous amplitude
                
                above_threshold = FastSpiAmplitudeTmp(5000:end-5000) > SpindleThreshold(2,1)*SpiAmpSTD(iCh);
                isLongEnough = bwareafilt(above_threshold, [250, MaxSpiDuration]); %find spindle within duration range
                %third threshold criteria
                above_Max = FastSpiAmplitudeTmp(5000:end-5000) > SpindleThreshold(3,1)*SpiAmpSTD(iCh);
                MaxIsThere = bwareafilt(above_Max, [1, MaxSpiDuration]); %find spindle within duration range
                [pks,locs] = findpeaks(DataTmpSpi(5000:end-5000,1),'MinPeakProminence',SpindleThreshold(1,1)*SpiAmpSTD(iCh));
                if sum(double(isLongEnough))>1 && sum(double(MaxIsThere))>1 && max(diff(locs))<100 %check if long enough spindle is present and check that no peak to peak distance is more than 125ms
                    
                else %if criteria not fullfilled store index of Spindles
                    TempIdx = [TempIdx iSpi];
                end
                
            end
            AllSpindles.FastSpi{iEpoch,iCh}(:,TempIdx)=[];%if not criteriy fullfilled delete detected spindle
        end
    end
    %calculated Spindel density
    TotalNumberOfSpi = 0;
    EpisodeDurations = 0;
    for iEpoch = 1:size(AllSpindles.FastSpi,1)
        CurrentSpindles = AllSpindles.FastSpi{iEpoch,iCh};
        TotalNumberOfSpi = TotalNumberOfSpi +size(CurrentSpindles,2);
        EpisodeDurations = EpisodeDurations + NREMEpisodes(2,iEpoch)-NREMEpisodes(1,iEpoch);
    end
    SpindleDensity = TotalNumberOfSpi/(EpisodeDurations/60);%spindle density in spindles per minute
    %%
    %detect SOs
    SOBand = zeros(size(Data,1),3);
    for iCh = 1:3
        SOBand(:,iCh) = filtfilt(SOFilterHigh1,SOFilterHigh2,Data(:,iCh));
        SOBand(:,iCh) = filtfilt(SOFilterLow1,SOFilterLow2,SOBand(:,iCh));
    end
    %calculate std for thresholds
    SOSTD = std(SOBand(SleepScoring==2,:));
    SOMean = mean(abs(SOBand(SleepScoring==2,:)));
    %find negative amplitudes bigger than Threshold
    SOEpisodes = cell(3,1);
    for iCh = 1:3
        SoThreshold = SOMean(iCh) + 0.66* SOSTD(iCh);
        for iEpoch = 1:size(NREMEpisodes,2)
            DataTmp = SOBand(NREMEpisodes(1,iEpoch)*Fs : NREMEpisodes(2,iEpoch)*Fs,iCh);
            
            SOBegEpisode = strfind((DataTmp<-SoThreshold)',[0 1])-1;
            SOEndEpisode = strfind((DataTmp<-SoThreshold)',[1 0]);
            if size(SOEndEpisode,1)>0
                if SOEndEpisode(1,1) < SOBegEpisode(1,1)
                    SOEndEpisode(:,1) = [];
                end
                if length(SOBegEpisode) > length(SOEndEpisode)
                    SOBegEpisode(:,end) = [];
                end
                SOEpisodes{iCh,1} = [SOEpisodes{iCh,1} [SOBegEpisode+NREMEpisodes(1,iEpoch)*Fs; SOEndEpisode+NREMEpisodes(1,iEpoch)*Fs]];
            end
        end
    end
    
    % find zero crossings
    ZeroCrossings = cell(3,1);
    for iCh = 1:3
        ZeroCrossings{iCh,1} = zeros(3,size(SOEpisodes{iCh,1},2));
        for iEvent = 1: size(SOEpisodes{iCh,1},2)
            X = 0;
            Y = 0;
            for iSearchCrossing = 1:2000
                if SOBand(SOEpisodes{iCh,1}(1,iEvent)-iSearchCrossing,iCh)>0 && X == 0
                    ZeroCrossings{iCh,1} (1,iEvent) = SOEpisodes{iCh,1}(1,iEvent)-iSearchCrossing;
                    X = 1;
                end
                if SOBand(SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing,iCh)>0 && Y == 0
                    ZeroCrossings{iCh,1} (2,iEvent) = SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing;
                    Y = 1;
                end
                if Y ==1 && SOBand(SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing,iCh)<0
                    ZeroCrossings{iCh,1}(3,iEvent) = SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing;
                    Y = 2;
                end
            end
        end
        ZeroCrossings{iCh,1}(:,find(ZeroCrossings{iCh,1}(1,:)==0))=[];
        ZeroCrossings{iCh,1}(:,find(ZeroCrossings{iCh,1}(2,:)==0))=[];
        ZeroCrossings{iCh,1}(:,find(ZeroCrossings{iCh,1}(3,:)==0))=[];
        %compensate if several peaks are below threshold but arent above 0 in
        %between and zerocrossings are the same
        tmp1 = unique(ZeroCrossings{iCh,1}(1,:));
        tmp2 = unique(ZeroCrossings{iCh,1}(2,:));
        tmp3 = unique(ZeroCrossings{iCh,1}(3,:));
        ZeroCrossings{iCh,1} = [tmp1; tmp2; tmp3];
        
        %remove SO with to long downstate
        ZeroCrossings{iCh,1}(:,(ZeroCrossings{iCh,1}(2,:)-ZeroCrossings{iCh,1}(1,:))>2000)=[];
        %remove SOs with to short duration
        ZeroCrossings{iCh,1}(:,(ZeroCrossings{iCh,1}(3,:)-ZeroCrossings{iCh,1}(1,:))<500)=[];
        
        %remove SOs with to small peak to peak amplitude
        NegPeakValue{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
        PosPeakValue{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
		Peak2PeakAmp{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
		for iEvent = 1:size(ZeroCrossings{iCh,1},2)
			NegPeakValue{iCh,1}(iEvent,1) = min(SOBand(ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent),iCh)',[],2);
			PosPeakValue{iCh,1}(iEvent,1) = max(SOBand(ZeroCrossings{iCh,1}(2,iEvent):ZeroCrossings{iCh,1}(3,iEvent),iCh)',[],2);
			Peak2PeakAmp{iCh,1}(iEvent,1) = abs(NegPeakValue{iCh,1}(iEvent,1))+PosPeakValue{iCh,1}(iEvent,1);
        end
        
        cfg.slo_rel_thr = 33;

        if isfield(cfg, 'slo_rel_thr')
            ZeroCrossings{iCh,1}(:,NegPeakValue{iCh,1}>prctile(NegPeakValue{iCh,1},cfg.slo_rel_thr)) = []; 
            PosPeakValue{iCh,1}(NegPeakValue{iCh,1}>prctile(NegPeakValue{iCh,1},cfg.slo_rel_thr),:) = [];
            Peak2PeakAmp{iCh,1}(NegPeakValue{iCh,1}>prctile(NegPeakValue{iCh,1},cfg.slo_rel_thr),:) = [];
            NegPeakValue{iCh,1}(NegPeakValue{iCh,1}>prctile(NegPeakValue{iCh,1},cfg.slo_rel_thr),:) = [];
            output.slo.thr =prctile(NegPeakValue{iCh,1},cfg.slo_rel_thr);
        else
            output.slo.thr = SoThreshold;
        end
		% Find negative peaks
		NegativePeaks{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
		for iEvent = 1: size(ZeroCrossings{iCh,1},2)
			[M,I] = min(SOBand(ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent),iCh)',[],2);
			NegativePeaks{iCh,1}(iEvent,1) = ZeroCrossings{iCh,1}(1,iEvent)+I;
        end
        
        
%         %remove SOs with to small peak to peak amplitude
%         Peak2PeakAmp{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
%         for iEvent = 1:size(ZeroCrossings{iCh,1},2)
%             NegPeakValue = min(SOBand(ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent),iCh),[],1);
%             PosPeakValue = max(SOBand(ZeroCrossings{iCh,1}(2,iEvent):ZeroCrossings{iCh,1}(3,iEvent),iCh),[],1);
%             Peak2PeakAmp{iCh,1}(iEvent,1) = abs(NegPeakValue)+PosPeakValue;
%         end
%         ZeroCrossings{iCh,1}(:,Peak2PeakAmp{iCh,1}<0.07) = [];
%         Peak2PeakAmp{iCh,1}(Peak2PeakAmp{iCh,1}<0.07,:) = [];
%         
%         %find negative peaks
%         NegativePeaks{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
%         for iEvent = 1: size(ZeroCrossings{iCh,1},2)
%             [M,I] = min(SOBand(ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent),iCh),[],1);
%             NegativePeaks{iCh,1}(iEvent,1) = ZeroCrossings{iCh,1}(1,iEvent)+I;
%         end
        %%
        %extract SO phase
        SOSpiCoupling{iCh,1} = zeros(size(NegativePeaks{iCh,1},1),1);
        SOGA{iCh,1} = zeros(size(NegativePeaks{iCh,1},1),5001);
        SOPhase{iCh,1} = zeros(size(NegativePeaks{iCh,1},1),5001);
        for iSO= 1:size(NegativePeaks{iCh,1},1)
            SOGA{iCh,1}(iSO,:) =  SOBand(NegativePeaks{iCh,1}(iSO,1)-2500:NegativePeaks{iCh,1}(iSO,1)+2500);
            SOPhase{iCh,1}(iSO,:) = rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:))));
            [TmpAmp, SpiAmpIndex] = max(SpiBandAmp(NegativePeaks{iCh,1}(iSO,1)-2500:NegativePeaks{iCh,1}(iSO,1)+2500));
            SOSpiCoupling{iCh,1}(iSO,1) = SOPhase{iCh,1}(iSO,SpiAmpIndex);
        end
        %find spindle band max
    end
    %%
    %calculated theta power during REM
    ThetaBand = zeros(size(Data,1),3);
    ThetaAmp = zeros(size(Data,1),3);
    for iCh = 1:3
        ThetaBand(:,iCh) = filtfilt(ThetaFilterHigh1,ThetaFilterHigh2,Data(:,iCh));
        ThetaBand(:,iCh) = filtfilt(ThetaFilterLow1,ThetaFilterLow2,ThetaBand(:,iCh));
        ThetaAmp(:,iCh) = abs(hilbert(ThetaBand(:,iCh)));
    end
    for iEpoch = 1:size(REMEpisodes,2)
        for iCh = 1:3
        REMThetaMeanAmp{iCh,1}(iEpoch,1) = mean(ThetaAmp(REMEpisodes(1,iEpoch):REMEpisodes(2,iEpoch),iCh));
        REMThetaEnergy{iCh,1}(iEpoch,1) = sum(ThetaAmp(REMEpisodes(1,iEpoch):REMEpisodes(2,iEpoch),iCh));
        end
    end
    if isempty(REMEpisodes)
        REMThetaMeanAmp{iCh,1}(iEpoch,1) = 0;
        REMThetaEnergy{iCh,1}(iEpoch,1) = 0;
    end

    save(strcat(DirSave,'DetectedEvents',AnimalName{iAnimal,1}),'AllSpindles','SOGA','SOSpiCoupling','NREMEpisodes','REMEpisodes','SleepScoring','REMThetaEnergy','REMThetaMeanAmp','-v7.3');
    %        save(strcat(DirSave,'DetectedEvents',AnimalName{iAnimal,1}),'','-v7.3');
    
    clear SOGA SOPhase SOSpiCoupling REMThetaMeanAmp REMThetaEnergy
    
end

%%
SOGA = zeros(size(NegativePeaks{1,1},1),251);
for iSO= 1:size(NegativePeaks{1,1},1)
   SOGA(iSO,:) =  downsample(Data(NegativePeaks{1,1}(iSO,1)-2500:NegativePeaks{1,1}(iSO,1)+2500),20);
%      plot(SOGA(iSO,:)');
%      waitforbuttonpress
end


iCh = 3;
GoodSpiIndex = 0;
for iEpoch = 1:size(AllSpindles.FastSpi,1)
    CurrentSpindles = AllSpindles.FastSpi{iEpoch,iCh};
    for iSpi = 1: size (CurrentSpindles,2)
        
        DataTmpSpi = SpiBand(CurrentSpindles(1,iSpi)-5000:CurrentSpindles(2,iSpi)+5000,iCh); %get filteres spindle signal for eachspindle + - 5sec
        FastSpiAmplitudeTmp = smooth(abs(hilbert(DataTmpSpi)),40);%get smoothed instantaneous amplitude
        
        above_threshold = FastSpiAmplitudeTmp(5000:end-5000) > SpindleThreshold(2,1)*SpiAmpSTD(iCh);
        isLongEnough = bwareafilt(above_threshold, [250, MaxSpiDuration]); %find spindle within duration range
        above_Max = FastSpiAmplitudeTmp(5000:end-5000) > SpindleThreshold(3,1)*SpiAmpSTD(iCh);
        MaxIsThere = bwareafilt(above_Max, [1, MaxSpiDuration]); %find spindle within duration range
        [pks,locs] = findpeaks(DataTmpSpi(5000:end-5000,1),'MinPeakProminence',1.5*SpiAmpSTD(iCh));
        if sum(double(isLongEnough))>1 && sum(double(MaxIsThere))>1 && max(diff(locs))<100 %check if long enough spindle is present and check that no peak to peak distance is more than 125ms
            GoodSpiIndex = GoodSpiIndex+1;
            subplot(2,1,1)
            DataTmp = Data(CurrentSpindles(1,iSpi)-5000:CurrentSpindles(1,iSpi)+5000,iCh);
            plot(DataTmp);
            ylim([-0.5 0.5]);
            subplot(2,1,2)
            DataTmp = filtfilt(SpiFilterHigh1,SpiFilterHigh2,DataTmp);
            DataTmp = filtfilt(SpiFilterLow1,SpiFilterLow2,DataTmp);
            plot(DataTmp);
            hold all
            DataTmp = smooth(abs(hilbert(DataTmp)),40);
            plot(DataTmp);
            %yline(5*FastSpiSTD);
            %yline(2.5*FastSpiSTD);
            yline(2.0*SpiAmpSTD(iCh));
            yline(1.5*SpiAmpSTD(iCh));
            yline(3*SpiAmpSTD(iCh));
            ylim([-0.3 0.3]);
            hold off
            waitforbuttonpress
        end
        
    end
    
end


for iCh=1:3
gradBin     = 360/10;
figure
obj2                    = CircHist(SOSpiCoupling{iCh,1},gradBin);
obj2.colorBar           = 'k';  % change color of bars
obj2.avgAngH.LineStyle  = '--'; % make average-angle line dashed
obj2.avgAngH.LineWidth  = 3; % make average-angle line thinner
obj2.colorAvgAng        = [1 0.5 0.5]; % change average-angle line color
end