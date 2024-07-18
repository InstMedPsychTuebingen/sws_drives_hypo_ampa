clear
close all

% DirData = 'Y:\niels\Jianfeng\EEGDataScored\Matlab\';
% DirSave = 'Y:\niels\Jianfeng\Results\PerHour\';
DirData = '/gpfs01/born/animal/niels/LunAMPAScaling/NewDataSet/MatlabEEGData/';
DirSave = '/gpfs01/born/animal/niels/LunAMPAScaling/NewDataSet/Results/PerHour/';
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
[FastSpiFilterHigh1,FastSpiFilterHigh2] = butter(3,2*10/Fs,'high');
[FastSpiFilterLow1,FastSpiFilterLow2] = butter(3,2*16/Fs,'low');

[SlowSpiFilterHigh1,SlowSpiFilterHigh2] = butter(3,2*7/Fs,'high');
[SlowSpiFilterLow1,SlowSpiFilterLow2] = butter(3,2*10/Fs,'low');

[SOFilterHigh1,SOFilterHigh2] = butter(3,2*0.1/Fs,'high');
[SOFilterLow1,SOFilterLow2] = butter(3,2*4/Fs,'low');
% SOFilter = designfilt('bandpassiir','FilterOrder',20, ...
%     'HalfPowerFrequency1',0.1,'HalfPowerFrequency2',4, ...
%     'SampleRate',Fs);

[ThetaFilterHigh1,ThetaFilterHigh2] = butter(3,2*5/Fs,'high');
[ThetaFilterLow1,ThetaFilterLow2] = butter(3,2*10/Fs,'low');
%%
for iAnimal = 1:size(AnimalName,1)
    iAnimal
    RawData = load(strcat(DirData,AnimalName{iAnimal,1}));
    RawChannelNames = fieldnames(RawData);
    RawScoring = RawData.SlStNew.codes(:,1);
    OrigScoringLength = length(RawScoring);
    %cut RawScoring to 6h
    if length(RawScoring) > (6*60*60)/10+1
        RawScoring = RawScoring(end-(6*60*60)/10+1:end,1);
    else
        Tmp = ones(((6*60*60)/10+1) - length(RawScoring),1,1);
        RawScoring = [Tmp; RawScoring];%compensate that this recording was stopped to early
    end
    
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
    Data(1,:) = [];%compensate that channels do not start simoultaneously
    %remove Data that exists without scoring
    DiffScoringData = OrigScoringLength*10*Fs - length(Data);
    Data(end-DiffScoringData:end,:) = [];
    %cut Data to the last 6h
    if length(RawScoring) > (6*60*60)/10+1
        Data = Data(end-(6*60*60*Fs)+1:end,:);
    else
        Tmp = zeros(((6*60*60*Fs)+1) - length(Data),5);
        Data = [Tmp; Data];
    end
    
     if NREMEpisodes(2,end)*Fs > size(Data,1)
        NREMEpisodes(2,end) = NREMEpisodes(2,end)-10;
    end
    
    %create upsampled scoring vector
    SleepScoring = zeros(size(Data,1),1);
    for iEpoch = 1:size(NREMEpisodes,2)
        SleepScoring(NREMEpisodes(1,iEpoch)*Fs:NREMEpisodes(2,iEpoch)*Fs) = 2;
    end
    if length(SleepScoring)>length(Data)
        SleepScoring = SleepScoring(1:length(Data),:);
    end
    %%
    %detect spindles
    SpiBand = zeros(size(Data,1),3);
    for iCh = 1:3
        SpiBand(:,iCh) = filtfilt(FastSpiFilterHigh1,FastSpiFilterHigh2,Data(:,iCh));
        SpiBand(:,iCh) = filtfilt(FastSpiFilterLow1,FastSpiFilterLow2,SpiBand(:,iCh));
    end
    SpiBandAmp = abs(hilbert(SpiBand));
    %calculate std for thresholds
    SpiAmpSTD = std(SpiBandAmp(SleepScoring==2,:));
    SpiAmpMean = mean(SpiBandAmp(SleepScoring==2,:));
    
    SOBand = zeros(size(Data,1),3);
    DeltaAmp = zeros(size(Data,1),3);
    for iCh = 1:3
        SOBand(:,iCh) = filtfilt(SOFilterHigh1,SOFilterHigh2,Data(:,iCh));
        SOBand(:,iCh) = filtfilt(SOFilterLow1,SOFilterLow2,SOBand(:,iCh));
        DeltaAmp(:,iCh) = abs(hilbert(SOBand(:,iCh)));
    end
    
    ThetaBand = zeros(size(Data,1),3);
    ThetaAmp = zeros(size(Data,1),3);
    for iCh = 1:3
        ThetaBand(:,iCh) = filtfilt(ThetaFilterHigh1,ThetaFilterHigh2,Data(:,iCh));
        ThetaBand(:,iCh) = filtfilt(ThetaFilterLow1,ThetaFilterLow2,ThetaBand(:,iCh));
        ThetaAmp(:,iCh) = abs(hilbert(ThetaBand(:,iCh)));
        
    end
    
    AllNREMEpisodes = NREMEpisodes;
    AllREMEpisodes = REMEpisodes;
    AllWAKEpisodes = WAKEpisodes;
    
    clear  ResultsPerHour
    ResultsPerHour =[];
    
    for ihour = 1:6 %repeat analysis for every hour
        clear NREMEpisodes REMEpisodes WAKEpisodes
        BegEpisode = AllNREMEpisodes(1,(AllNREMEpisodes(1,:)<ihour*60*60 & AllNREMEpisodes(1,:)>(ihour-1)*60*60));
        EndEpisode = AllNREMEpisodes(2,(AllNREMEpisodes(2,:)<ihour*60*60 & AllNREMEpisodes(2,:)>(ihour-1)*60*60));
        if length(BegEpisode)<length(EndEpisode)
            BegEpisode = [(ihour-1)*60*60+1 BegEpisode];
        end
        if length(BegEpisode)>length(EndEpisode)
            EndEpisode = [EndEpisode ihour*60*60];
        end
        if ~isempty(BegEpisode) || ~isempty(EndEpisode)
            if BegEpisode(1,end)> EndEpisode(1,end)
                EndEpisode = [EndEpisode ihour*60*60];
            end
            if BegEpisode(1,1)> EndEpisode(1,1)
                BegEpisode = [(ihour-1)*60*60+1 BegEpisode];
            end
        end
        NREMEpisodes = [BegEpisode; EndEpisode];
        if ~isempty(AllREMEpisodes)
        BegEpisode = AllREMEpisodes(1,(AllREMEpisodes(1,:)<ihour*60*60 & AllREMEpisodes(1,:)>(ihour-1)*60*60));
        EndEpisode = AllREMEpisodes(2,(AllREMEpisodes(2,:)<ihour*60*60 & AllREMEpisodes(2,:)>(ihour-1)*60*60));
        if length(BegEpisode)<length(EndEpisode)
            BegEpisode = [(ihour-1)*60*60+1 BegEpisode];
        end
        if length(BegEpisode)>length(EndEpisode)
            EndEpisode = [EndEpisode ihour*60*60];
        end
        if ~isempty(BegEpisode) || ~isempty(EndEpisode)
            if BegEpisode(1,end)> EndEpisode(1,end)
                EndEpisode = [EndEpisode ihour*60*60];
            end
            if BegEpisode(1,1)> EndEpisode(1,1)
                BegEpisode = [(ihour-1)*60*60+1 BegEpisode];
            end
        end
        REMEpisodes = [BegEpisode; EndEpisode];
        else
        REMEpisodes = [];
        end
        BegEpisode = AllWAKEpisodes(1,(AllWAKEpisodes(1,:)<ihour*60*60 & AllWAKEpisodes(1,:)>(ihour-1)*60*60));
        EndEpisode = AllWAKEpisodes(2,(AllWAKEpisodes(2,:)<ihour*60*60 & AllWAKEpisodes(2,:)>(ihour-1)*60*60));
        if length(BegEpisode)<length(EndEpisode)
            BegEpisode = [(ihour-1)*60*60+1 BegEpisode];
        end
        if length(BegEpisode)>length(EndEpisode)
            EndEpisode = [EndEpisode ihour*60*60];
        end
        if ~isempty(BegEpisode) || ~isempty(EndEpisode)
            if BegEpisode(1,end)> EndEpisode(1,end)
                EndEpisode = [EndEpisode ihour*60*60];
            end
            if BegEpisode(1,1)> EndEpisode(1,1)
                BegEpisode = [(ihour-1)*60*60+1 BegEpisode];
            end
        end
        WAKEpisodes = [BegEpisode; EndEpisode];
        
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
                %%
                %second threshold criteria
                if size(AllSpindles.FastSpi{iEpoch,iCh},1) >0
                    AllSpindles.FastSpi{iEpoch,iCh}(:,AllSpindles.FastSpi{iEpoch,iCh}(1,:)>length(Data)-5000) =[]; %delete events to close to ending
                end
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
        for iCh=1:3
            TotalNumberOfSpi{iCh} = 0;
            EpisodeDurations{iCh} = 0;
            for iEpoch = 1:size(AllSpindles.FastSpi,1)
                CurrentSpindles = AllSpindles.FastSpi{iEpoch,iCh};
                TotalNumberOfSpi{iCh} = TotalNumberOfSpi{iCh} +size(CurrentSpindles,2);
                EpisodeDurations{iCh} = EpisodeDurations{iCh} + NREMEpisodes(2,iEpoch)-NREMEpisodes(1,iEpoch);
            end
            SpindleDensity{iCh} = TotalNumberOfSpi{iCh}/(EpisodeDurations{iCh}/60);%spindle density in spindles per minute
            ResultsPerHour{ihour}.SpindleDensity{iCh,1} = SpindleDensity{iCh};
            ResultsPerHour{ihour}.TotalNumberOfSpi{iCh,1} = TotalNumberOfSpi{iCh};
        end
        %%
        %detect SOs 
        for iEpoch = 1:size(NREMEpisodes,2)
            for iCh = 1:3
                DeltaMeanAmp{iCh,1}(iEpoch,1) = mean(DeltaAmp(NREMEpisodes(1,iEpoch)*Fs:NREMEpisodes(2,iEpoch)*Fs,iCh));
                DeltaEnergy{iCh,1}(iEpoch,1) = sum(DeltaAmp(NREMEpisodes(1,iEpoch)*Fs:NREMEpisodes(2,iEpoch)*Fs,iCh));
            end
        end
        
        ResultsPerHour{ihour}.DeltaMeanAmp = DeltaMeanAmp;
        ResultsPerHour{ihour}.DeltaEnergy = DeltaEnergy;
        DeltaMeanAmp = [];
        DeltaEnergy = [];
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
            
            %remove SO with to to long downstate
            ZeroCrossings{iCh,1}(:,(ZeroCrossings{iCh,1}(2,:)-ZeroCrossings{iCh,1}(1,:))>2000)=[];
            %remove SOs with two short duration
            ZeroCrossings{iCh,1}(:,(ZeroCrossings{iCh,1}(3,:)-ZeroCrossings{iCh,1}(1,:))<500)=[];
            %remove SOs with to small peak to peak amplitude
            Peak2PeakAmp{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
            for iEvent = 1:size(ZeroCrossings{iCh,1},2)
                [NegPeakValue, IMin] = min(SOBand(ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent),iCh),[],1);
                [PosPeakValue, IMax] = max(SOBand(ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(3,iEvent),iCh),[],1);
                SOSlope{iCh,1}(iEvent,1) = (abs(NegPeakValue)+PosPeakValue)/(IMax - IMin);
                Peak2PeakAmp{iCh,1}(iEvent,1) = abs(NegPeakValue)+PosPeakValue;
            end
            ZeroCrossings{iCh,1}(:,Peak2PeakAmp{iCh,1}<0.07) = [];
            Peak2PeakAmp{iCh,1}(Peak2PeakAmp{iCh,1}<0.07,:) = [];
            SOSlope{iCh,1}(Peak2PeakAmp{iCh,1}<0.07,:) = [];
            ResultsPerHour{ihour}.Peak2PeakAmp{iCh,1} = Peak2PeakAmp{iCh,1};
            ResultsPerHour{ihour}.SOSlope{iCh,1} = SOSlope{iCh,1};
            %find negative peaks
            NegativePeaks{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
            for iEvent = 1: size(ZeroCrossings{iCh,1},2)
                [M,I] = min(SOBand(ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent),iCh),[],1);
                NegativePeaks{iCh,1}(iEvent,1) = ZeroCrossings{iCh,1}(1,iEvent)+I;
            end            
            %%
            %extract SO phase
%             SOSpiCoupling{iCh,1} = zeros(size(NegativePeaks{iCh,1},1),1);
            SOGA{iCh,1} = zeros(size(NegativePeaks{iCh,1},1),5001);
%             SOPhase{iCh,1} = zeros(size(NegativePeaks{iCh,1},1),5001);
            iSOTmp = 0;
            
            SloSpiDetCoupling{iCh,1}	= [];
            SOSpiCoupling{iCh,1}        = [];
			SOPhase						= [];
			cnt							= 1;
            twindow = 2.5;
			for iSO = 1:size(NegativePeaks{iCh,1},1)
				spi_cur = [];
				% If there is a spindle fully inside SO plus minus time
				% window
                for iEpoch = 1: size(AllSpindles.FastSpi,1)
                    if ~isempty(AllSpindles.FastSpi{iEpoch,iCh})
                        spi_ind = find(AllSpindles.FastSpi{iEpoch,iCh}(1,:) > NegativePeaks{iCh,1}(iSO,1)-round(twindow*Fs) & AllSpindles.FastSpi{iEpoch,iCh}(2,:) < NegativePeaks{iCh,1}(iSO,1)+round(twindow*Fs));
                        if size(spi_ind,2) > 0 % if at least one spindle was found
                            spi_cur = [spi_cur AllSpindles.FastSpi{iEpoch,iCh}(:,spi_ind)];
                        end
                    end
                end
                SOGA{iCh,1}(iSO,:) =  SOBand(NegativePeaks{iCh,1}(iSO,1)-2500:NegativePeaks{iCh,1}(iSO,1)+2500,iCh)';
				for iSpi = 1:size(spi_cur,2)
					SOPhase = rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:)))); % SO phase along entire window
					[~, SpiAmpIndex] = max(SpiBandAmp(spi_cur(1,iSpi):spi_cur(2,iSpi),iCh)); % find spindle maximum amp (samples from spindle start)
					tmp = spi_cur(1,iSpi) + SpiAmpIndex - 1; % spindle maximum amp in global samples
					tmp = tmp - (NegativePeaks{iCh,1}(iSO,1)-round(twindow*Fs)); % spindle maximum amp in samples from SO window start
					SOSpiCoupling{iCh,1}(cnt,1) = SOPhase(round(tmp)); % note down phase there
					cnt = cnt + 1;
                end
                
			end
            if size(SOSpiCoupling,1)<iCh
                SOSpiCoupling{iCh,1} = [];
            end
%             
%             for iSO= 1:size(NegativePeaks{iCh,1},1)
%                 y = 0;
%                 CurrentSpindle = [];
%                 for iEpoch = 1: size(AllSpindles.FastSpi,1)
%                     if size(AllSpindles.FastSpi{iEpoch,iCh},2) >0
%                         if size(find(AllSpindles.FastSpi{iEpoch,iCh}(1,:)*((ihour-1)*60*60*Fs+1) > NegativePeaks{iCh,1}(iSO,1)-2500 & AllSpindles.FastSpi{iEpoch,iCh}(1,:)*((ihour-1)*60*60*Fs+1) < NegativePeaks{iCh,1}(iSO,1)+2500),2) > 0 && y == 0
%                             y=1;
%                             CurrentSpindle = AllSpindles.FastSpi{iEpoch,iCh}(:,(AllSpindles.FastSpi{iEpoch,iCh}(1,:)*((ihour-1)*60*60*Fs+1) > NegativePeaks{iCh,1}(iSO,1)-2500 & AllSpindles.FastSpi{iEpoch,iCh}(1,:)*((ihour-1)*60*60*Fs+1) < NegativePeaks{iCh,1}(iSO,1)+2500));                           
%                         end
%                     end
%                 end
%                 %                 SOGA{iCh,1}(iSO,:) =  SOBand(NegativePeaks{iCh,1}(iSO,1)-2500:NegativePeaks{iCh,1}(iSO,1)+2500);
%                 %                 SOPhase{iCh,1}(iSO,:) = rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:))));
%                 %                 [TmpAmp, SpiAmpIndex] = max(SpiBandAmp(NegativePeaks{iCh,1}(iSO,1)-2500:NegativePeaks{iCh,1}(iSO,1)+2500));
%                 %                 SOSpiCoupling{iCh,1}(iSO,1) = SOPhase{iCh,1}(iSO,SpiAmpIndex);
%                 if size(CurrentSpindle,2)>0
%                     for iSpi = 1:size(CurrentSpindle,2)
%                         iSOTmp = iSOTmp + 1;
%                         SOGA{iCh,1}(iSO,:) =  SOBand(NegativePeaks{iCh,1}(iSO,1)-2500:NegativePeaks{iCh,1}(iSO,1)+2500);
%                         SOPhase{iCh,1}(iSO,:) = rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:))));
%                         [TmpAmp, SpiAmpIndex] = max(SpiBandAmp(CurrentSpindle(1,iSpi):CurrentSpindle(2,iSpi)));
%                         SOSpiCoupling{iCh,1}(iSOTmp,1) = SOPhase{iCh,1}(iSO,SpiAmpIndex);
%                     end
%                 end
%             end
            ResultsPerHour{ihour}.SOSpiCoupling{iCh,1} = SOSpiCoupling{iCh,1};
            SOSpiCoupling{iCh,1} = [];
            ResultsPerHour{ihour}.SONumber{iCh,1} = size(NegativePeaks{iCh,1},1);
        end
        %%
        %calculated theta power during REM
        for iCh = 1:3
            REMThetaMeanAmp{iCh,1} = zeros(size(REMEpisodes,2),1);
            REMThetaEnergy{iCh,1} = zeros(size(REMEpisodes,2),1);
        end
        
        for iEpoch = 1:size(REMEpisodes,2)
            for iCh = 1:3
                REMThetaMeanAmp{iCh,1}(iEpoch,1) = mean(ThetaAmp(REMEpisodes(1,iEpoch)*Fs:REMEpisodes(2,iEpoch)*Fs,iCh));
                REMThetaEnergy{iCh,1}(iEpoch,1) = sum(ThetaAmp(REMEpisodes(1,iEpoch)*Fs:REMEpisodes(2,iEpoch)*Fs,iCh));
            end
        end
        ResultsPerHour{ihour}.REMThetaMeanAmp = REMThetaMeanAmp;
        ResultsPerHour{ihour}.REMThetaEnergy = REMThetaEnergy;
        
        %%
        %calculated theta power during NREM
        for iCh = 1:3
            NREMThetaMeanAmp{iCh,1} = zeros(size(NREMEpisodes,2),1);
            NREMThetaEnergy{iCh,1} = zeros(size(NREMEpisodes,2),1);
        end
        
        for iEpoch = 1:size(NREMEpisodes,2)
            for iCh = 1:3
                NREMThetaMeanAmp{iCh,1}(iEpoch,1) = mean(ThetaAmp(NREMEpisodes(1,iEpoch)*Fs:NREMEpisodes(2,iEpoch)*Fs,iCh));
                NREMThetaEnergy{iCh,1}(iEpoch,1) = sum(ThetaAmp(NREMEpisodes(1,iEpoch)*Fs:NREMEpisodes(2,iEpoch)*Fs,iCh));
            end
        end
        ResultsPerHour{ihour}.NREMThetaMeanAmp = NREMThetaMeanAmp;
        ResultsPerHour{ihour}.NREMThetaEnergy = NREMThetaEnergy;
        %%
        if ~isnan(mean(NREMEpisodes(2,:)'-NREMEpisodes(1,:)'))
            ResultsPerHour{ihour}.MeanNREMDur = mean(NREMEpisodes(2,:)'-NREMEpisodes(1,:)');
        else
            ResultsPerHour{ihour}.MeanNREMDur = 0;
        end
        if ~isempty(REMEpisodes)
        if ~isnan(mean(REMEpisodes(2,:)'-REMEpisodes(1,:)'))
            ResultsPerHour{ihour}.MeanREMDur = mean(REMEpisodes(2,:)'-REMEpisodes(1,:)');
        else
            ResultsPerHour{ihour}.MeanREMDur = 0;
        end
        else
        ResultsPerHour{ihour}.MeanREMDur = 0;
        end
        if ~isnan(mean(WAKEpisodes(2,:)'-WAKEpisodes(1,:)'))
            ResultsPerHour{ihour}.MeanWAKDur = mean(WAKEpisodes(2,:)'-WAKEpisodes(1,:)');
        else
            ResultsPerHour{ihour}.MeanWAKDur =  0;
        end
        
        if ~isnan(sum(NREMEpisodes(2,:)'-NREMEpisodes(1,:)'))
            ResultsPerHour{ihour}.AbsNREMDur = sum(NREMEpisodes(2,:)'-NREMEpisodes(1,:)');
        else
            ResultsPerHour{ihour}.AbsNREMDur = 0;
        end
        if ~isempty(REMEpisodes)
        if ~isnan(sum(REMEpisodes(2,:)'-REMEpisodes(1,:)'))
            ResultsPerHour{ihour}.AbsREMDur = sum(REMEpisodes(2,:)'-REMEpisodes(1,:)');
        else
            ResultsPerHour{ihour}.AbsREMDur = 0;
        end
        else
            ResultsPerHour{ihour}.AbsREMDur = 0;
        end
        if ~isnan(sum(WAKEpisodes(2,:)'-WAKEpisodes(1,:)'))
            ResultsPerHour{ihour}.AbsWAKDur = sum(WAKEpisodes(2,:)'-WAKEpisodes(1,:)');
        else
            ResultsPerHour{ihour}.AbsWAKDur =  0;
        end        
    end
    
    save(strcat(DirSave,'ResultsPerHour',AnimalName{iAnimal,1}),'ResultsPerHour','-v7.3');
    
    clear SOGA SOPhase SOSpiCoupling REMThetaMeanAmp REMThetaEnergy
    
end

