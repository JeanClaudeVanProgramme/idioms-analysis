function[exclude_trial_list, errorType] = ExclusionTrialEval(resp_1)
%List of stimuli from Sobh that are unreliable
stimuli_to_exclude={'~bio1-';
    '~cul9-';
    '~ins14-';
    '~ins16-';
    '~ins4-';
    '~ins7-';
    '~ins8-';
    '~ins9-';
    '~met1-';
    '~met13-';
    '~met14-';
    '~met17-';
    '~met2-';
    '~met3-';
    '~met7-'; };


%do a string comparison
for qq=1:size(resp_1,1)
    for nn=1:length(stimuli_to_exclude)
        tmp1=stimuli_to_exclude{nn};
        tmp2=resp_1{qq};
        if(~isempty(tmp2))
            exclude_trial_list(qq,nn)=strcmp(tmp1,tmp2(1:length(tmp1)));
        else %no response? experiment error
            errorType=2;
            exclude_trial_list(qq,nn)=1;
        end
    end
    
end

exclude_trial_list=sum(exclude_trial_list,2)>0;

return
