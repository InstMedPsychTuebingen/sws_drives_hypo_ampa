 
clear; close all;

% dir and load data 
cd add_path
load('output_events_detected.mat');

%% global config
exp       = fieldnames(Dat);
params    = {'AbsNREMDur','SpindleDensity','DeltaMeanAmp','DeltaEnergy',...
  'SODensity','AbsREMDur','REMThetaMeanAmp','REMThetaEnergy'};
para_corr ={'SWS time (min)','Spindle density (#/min)','SWS power','SWA energy',...
  'SO density (#/min)','REM time (min)','REM theta power','REM theta energy'};
hours     = 4 : 6; % 2nd 3 hours
gs        = fieldnames(Dat.(exp{1}));
igs    = find(contains(fieldnames(Dat.(exp{1})),'SLEEP'));
reg       = fieldnames(Dat.(exp{1}).(gs{igs}).WBResults);
reg_idx   = 1;
corr_type = 'Pearson'; % options: 'Pearson';% 'Spearman';

%% start for loops
clear res
for ipara = 1: numel(params)
  for iexp=1 : 2  
    for reg_idx = 1 : 2 

      if contains(params{ipara},{'Mean','Density','Peak2Peak','Slope','Coupling'})
        if contains(params{ipara},'Peak2Peak')
          tmp_para = [...
            nanmean(Dat.(exp{iexp}).(gs{igs}).(params{ipara}){1,1}(:,hours(1,1)),2) - ...
            nanmean(Dat.(exp{iexp}).(gs{igs}).(params{ipara}){1,1}(:,hours(1,end)),2); ];
        else
          tmp_para = [...
            nanmean(Dat.(exp{iexp}).(gs{igs}).(params{ipara}){1,1}(:,hours),2) ];
        end
      elseif contains(params{ipara},{'Energy','Number','Abs'})
        if contains(params{ipara},'Peak2Peak')
          tmp_para = [...
            nansum(Dat.(exp{iexp}).(gs{igs}).(params{ipara}){1,1}(:,hours(1,1)),2) - ...
            nansum(Dat.(exp{iexp}).(gs{igs}).(params{ipara}){1,1}(:,hours(1,end)),2) ];
        else
          tmp_para = [...
            nansum(Dat.(exp{iexp}).(gs{igs}).(params{ipara}){1,1}(:,hours),2) ];
        end
      end % if contains

      tmp_wb =  Dat.(exp{iexp}).(gs{igs}).WBResults.(reg{reg_idx})(:,1)  ;

      [r1,p1] = corr(tmp_para,tmp_wb,'Type',corr_type);
      res.(exp{iexp}).(gs{igs}).(params{ipara}){reg_idx,1} = reg{reg_idx};
      res.(exp{iexp}).(gs{igs}).(params{ipara}){reg_idx,2} = ['r=',num2str(r1,3),'; p=',num2str(p1,3)];
      res.(exp{iexp}).(gs{igs}).(params{ipara}){reg_idx,3} = para_corr{ipara};

    end % for reg_idx
  end % for iexp
end % for ipara
%% save
save('output_corr_sleep_param.mat','res','-v7.3')

