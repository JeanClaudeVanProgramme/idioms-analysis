clear all; clc;

visualize=0;
xy_res=[1280 1024];
smp_rate=1/120;
trialColorCode='rgkm';
errorType=0;

AOI_Hit_Threshold=0.4; %above what proportion of data must have an AOI hit?
tloss_threshold=0.7; % above what percent will we accept track losses

%define AOIS
box_sz=256*2; %border in & outside of screen
%rectangle uses bottom left + width & height
AOI(1,:)=[-box_sz*.25 xy_res(2)-box_sz*.75 box_sz box_sz];
AOI(2,:)=[xy_res(1)-box_sz*.75 xy_res(2)-box_sz*.75 box_sz box_sz];
AOI(3,:)=[xy_res(1)-box_sz*.75 -box_sz*.25 box_sz box_sz];
AOI(4,:)=[-box_sz*.25 -box_sz*.25 box_sz box_sz];

%old AOI using 256 boxes
%AOI(1,:)=[0 xy_res(2)-box_sz box_sz box_sz];
%AOI(2,:)=[xy_res(1)-box_sz xy_res(2)-box_sz box_sz box_sz];
%AOI(3,:)=[xy_res(1)-box_sz 0 box_sz box_sz];
%AOI(4,:)=[0 0 box_sz box_sz];

%load the master subject list
%subject_summary=readtable('AllSubjectSummary.txt','Delimiter','\t');
subject_summary=readtable('table idioms_3.xlsx');
%inclusionList;

datadir='C:\Users\bsun\Documents\GitHub\idioms-data\data\';
outputdir='\output\';

%read in the data directories and sort
filedirArray=dir(datadir);
idx=1;
for zz=1:length(filedirArray)
    %matlab insists on reading in hiddne directories, ignore them
    if( ~strcmp(filedirArray(zz).name,'.') & ...
            ~strcmp(filedirArray(zz).name,'..')  )
        fileDirList(idx)=str2num(filedirArray(zz).name);
        idx=idx+1;
    end
end

%fileDirList=sort(fileDirList);
fileDirList=unique(subject_summary.SubjectID);

AllSub_Details=zeros(length(fileDirList),2)-1;
AllSub_ave_AOI_per_idiom=zeros(4,10001,4,length(fileDirList)) +nan; % change 1st Dimension from 4 to 8 if using text & audio

for j=1:length(fileDirList) %per subject %first subjbect has weird eye data - skip
    n=0;
    fdir= num2str(fileDirList(j));
    
    sub_rng=find(subject_summary.SubjectID==str2num(fdir));
    
    filelistArray=dir([ datadir fdir '\*.csv']);
    [~,idx] = sort([filelistArray.datenum]);
    for zz= 1:length(idx)
        filelistCSV{zz}=filelistArray(idx(zz)).name;
    end
    
    %load mouse data from both blocks and combine together
    LoadMouseData;
    
    %open 2 files one for results the other for possible errors
    startDatafile;
    
    %load list of unreliable stimuli --  is this still valid list?  
    %exclude_trial_list=ExclusionTrialEval(resp_1);

    %evaluate if the trial belongs on inclusion list
    %instead of using inclusion list from inclusionList.m, use Sobh's excel
    %list
    inclusion_list=subject_summary.Item(sub_rng);
    keepTrial=zeros(length(condition),1);

    for bb=1:length(inclusion_list)
        tmp=strcmp(condition, inclusion_list(bb))>0;
        idx=find(tmp);
        if(isempty(idx))
           % keyboard
           disp('Entry in Sobhs list is not in data file!')
        else
            keepTrial(idx)=1;
        end
    end
    % a
    exclude_trial_list=~keepTrial'; %exclude_trial_list | ~keepTrial'; %exclude list is old
    sub_details=[-1 -1 -1 -1];

    
    for trial=1:length(filelistCSV) %eaxch csv trial file
        
        %test if this is a excluded trial
        if(exclude_trial_list(trial)) % | RT(trial)>=18000 ) RT's should already be screened
            %writeData; for now dont do anything, maybe have some entry for
            %skupped trials later?
            continue; %skip to next trial
        end
        
        fn=filelistCSV{trial}
        
        %test what stimili were in what location
        %resp 1-4 correspond to top left, top right, bottom right, bottom left
        r_tmp1=resp_1{trial}; r_tmp1=r_tmp1(end-4);
        r_tmp2=resp_2{trial}; r_tmp2=r_tmp2(end-4);
        r_tmp3=resp_3{trial}; r_tmp3=r_tmp3(end-4);
        r_tmp4=resp_4{trial}; r_tmp4=r_tmp4(end-4);
        
        %stimulus order: each entry corresponds to 1 to 4 onscreen
        stimulus_order=[r_tmp1 r_tmp2 r_tmp3 r_tmp4];
        %we want to read out as a code A=1, B=2, could do a if& then but faster
        %with sorting
        
        %stimulus_location corresponds to a standard ABCD order and indicates
        %which spot onscreen that stimulus appeared at
        % for instance 2 3 4 1 means A was in spot 2, B in spot 3, etc.
        %later we use stimulus location to give AOIs semantic labels
        [~,stimulus_location_tmp]=sort(stimulus_order);
        %sort a 2nd time to get in in terms of slot 1 is spatial and the entry
        %is stumulus
        
        [~,AOI_label]=sort(stimulus_location_tmp);
        
        %remap mouse Y data as it is upside down
        mouse(trial).Y=xy_res(2)-mouse(trial).Y;
        
        %AOI detection
        [mouse_AOI_PositionList, mouse_AOI_StimulusList]=...
            IsInAOI(mouse(trial).X,mouse(trial).Y,AOI,AOI_label);
        
        %mouse(trial).AOI_PositionList=mouse_AOI_PositionList;
        mouse(trial).AOI_StimulusList=mouse_AOI_StimulusList;
        
        %remove last header column is empty (why?)
        fid = fopen([ datadir fdir  '\' fn], 'r');
        str = fgetl(fid);
        fclose(fid);
        vars = regexp(str, '\t', 'split');
        vars = vars(1:end-1);
        eyedat=readtable([ datadir fdir  '\' fn] ,'ReadVariableNames',false,'HeaderLines',1,'Delimiter','\t');
        validVars = matlab.lang.makeValidName(vars);
        eyedat.Properties.VariableNames = validVars;
        
        
        %if we assume that RT has the duration from picture appearance to mouse
        %click we can trim the eye tracking data to the essential part
        rng=(size(eyedat,1)-round( (RT(trial)/1000)/smp_rate )):size(eyedat,1);
        if(rng(1)<1) %negative numbers...not sure why this would happen
            fprintf(fileID2, '%s\t%d Trial too short', fn);
            continue; %skip to next trial
        end
        
        %---Valid trial tests
        RL_eye_X_empty(trial,:)=[sum(eyedat.GazePoint2d_Right_x==-1) ...
            sum(eyedat.GazePoint2d_Left_x==-1)]/size(eyedat,1);
        RL_eye_Y_empty(trial,:)=[sum(eyedat.GazePoint2d_Right_y==0) ...
            sum(eyedat.GazePoint2d_Left_y==0)]/size(eyedat,1);
        %---End Valid trial tests
        
        %remap gaze data to screen resolution
        eyeL_X=eyedat.GazePoint2d_Left_x(rng)*xy_res(1);
        %Why do we need to transform the Y-data? Should be 0 to 1, not -100 to 100
        %also flip upside down as Tobii coords are 0,0 == Top Left
        eyeL_Y=xy_res(2)-((eyedat.GazePoint2d_Left_y(rng)+100)/200)*xy_res(2);
        
        eyeR_X=eyedat.GazePoint2d_Right_x(rng)*xy_res(1);
        eyeR_Y=xy_res(2)-((eyedat.GazePoint2d_Right_y(rng)+100)/200)*xy_res(2);
        
        %apply median filter to eye position
        eyeL_X=medfilter1(eyeL_X,4);
        eyeR_X=medfilter1(eyeR_X,4);
        eyeL_Y=medfilter1(eyeL_Y,4);
        eyeR_Y=medfilter1(eyeR_Y,4);
        eyeBinocXY = [(eyeL_X+eyeL_X)/2; (eyeL_Y+eyeR_Y)/2]';
        
        % [eye_AOI_PositionList, eye_AOI_StimulusList]=...
        [~, eye_AOI_StimulusList]=...
            IsInAOI(eyeBinocXY(:,1),eyeBinocXY(:,2),AOI, AOI_label);
        
        
        t=(eyedat.EyeTrackerTimestamp(rng)-eyedat.EyeTrackerTimestamp(rng(1)))/1000000;
        
        if(visualize)
            figure(1); clf;
            subplot(3,2,1);  hold on;
            
            %plot(eyeL_X, eyeL_Y,'b*')
            %plot(eyeR_X, eyeR_Y,'r*')
            plot(eyeBinocXY(:,1), eyeBinocXY(:,2),'*')
            DrawAOIs(AOI,xy_res);
            %legend('Left','Right');
            %color eye according to AOI detection
            
            for jj=1:4
                %idx=find(eye_AOI_PositionList==jj);
                idx=find(eye_AOI_StimulusList==jj);
                
                plot(eyeBinocXY(idx,1),eyeBinocXY(idx,2),[trialColorCode(jj) '*'])
            end
            xlim([0-box_sz*.25 xy_res(1)+box_sz*.25])
            ylim([0-box_sz*.25 xy_res(2)+box_sz*.25])
            xlabel('Horizontal Position')
            ylabel('Verical Position')
            title('Eye XY')
            legend('None','A','B','C','D')
            
            subplot(3,2,3);  hold on;
            %Why do we need to transform the Y-data? Should be 0 to 1, not -100 to 100
            plot(t,eyeL_X,'b*')
            plot(t,eyeR_X,'r*')
            %legend('Left','Right')
            title('Eye Horizontal')
            ylabel('Horizontal Position')
            xlabel('Time(s)')
            
            subplot(3,2,5); hold on;
            %Why do we need to transform the Y-data? Should be 0 to 1, not -100 to 100
            plot(t,eyeL_Y,'b*')
            plot(t,eyeR_Y,'r*')
            %legend('Left','Right')
            title('Eye Vertical')
            ylabel('Verical Position')
            xlabel('Time (s)')
            
            %---Plot mouse data---
            subplot(3,2,2); hold on; plot(mouse(trial).X,mouse(trial).Y,'*')
            DrawAOIs(AOI,xy_res);
            
            %color mouse according to AOI detection
            for jj=1:4
                %idx=find( mouse(trial).AOI_PositionList==jj);
                idx=find( mouse(trial).AOI_StimulusList==jj);
                
                plot(mouse(trial).X(idx),mouse(trial).Y(idx),[trialColorCode(jj) '*'])
            end
            
            
            xlim([0 xy_res(1)])
            ylim([0 xy_res(2)])
            title('Mouse XY Position')
            xlabel('Horizontal Position')
            ylabel('Verical Position')
            
            subplot(3,2,4); plot(mouse(trial).T/1000,mouse(trial).X)
            title('Mouse X Position vs. Time')
            ylabel('Horizontal Position')
            xlabel('Time(s)')
            
            subplot(3,2,6); plot(mouse(trial).T/1000,mouse(trial).Y)
            title('Mouse Y Position vs. Time')
            xlabel('Verical Position')
            ylabel('Time (s)')
            
            %animate in real time
            if(0)
                figure(1); subplot(3,2,1);
                for jj=1:length(eyeBinocXY)
                    DrawAOIs(AOI,xy_res);
                    plot(eyeBinocXY(jj,1),eyeBinocXY(jj,2),'r*')
                    xlim([0 xy_res(1)])
                    ylim([0 xy_res(2)])
                    xlabel('Horizontal Position')
                    ylabel('Verical Position')
                    title('Eye XY')
                    pause(0.008)
                end
            end
            
            figure(2); clf;
            subplot(2,1,1); hold on
            plot(t,eyeBinocXY(:,2),'b*')
            plot(t,eyeBinocXY(:,2),'b*')
            subplot(2,1,2); hold on
            plot(mouse(trial).T/1000,mouse(trial).X,'b*')
            plot(mouse(trial).T/1000,mouse(trial).Y,'b*')
            
            for jj=1:4
                %idx=find(eye_AOI_PositionList==jj);
                idx=find(eye_AOI_StimulusList==jj);
                
                subplot(2,1,1);
                title('Mouse vs. Eye AOI v Time Sync, Colors=AOIs, blue=no label')
                plot(t(idx),eyeBinocXY(idx,2),[trialColorCode(jj) '*'])
                plot(t(idx),eyeBinocXY(idx,2),[trialColorCode(jj) '*'])
                ylabel('Eye X Position')
                xlabel('Time (s)')
                
                subplot(2,1,2);
                %idx=find(mouse(trial).AOI_PositionList==jj);
                idx=find(mouse(trial).AOI_StimulusList==jj);
                
                plot(mouse(trial).T(idx)/1000,mouse(trial).X(idx),[trialColorCode(jj) '*'])
                plot(mouse(trial).T(idx)/1000,mouse(trial).Y(idx),[trialColorCode(jj) '*'])
                ylabel('Mouse Position')
                xlabel('Time (s)')
            end
        end
        
        time_from_samples(trial)=size(eyedat,1)*smp_rate;
        time_from_tracker(trial)=(eyedat.EyeTrackerTimestamp(end)-...
            eyedat.EyeTrackerTimestamp(1))/1000000;
        
        [mouse_distance_traveled,mouse_num_zero_crossings,...
            mouse_num_AOIs_visited, mouse_peakVel,...
            mouse_AOIs_visited,errorType]=mouseMetrics(mouse(trial),init_time(trial));
        
        [AOI_segments,AOI_total_dur,AOI_num_visits, AOI_mean_dur]=...
            eyeMetrics(eye_AOI_StimulusList,smp_rate);
        
        
        AOI_hit_perc(trial)=sum(eye_AOI_StimulusList>0)/length(eye_AOI_StimulusList);
        %if AOIs aren't looked at trial is probably junk - consider other
        %ways to filter bad trials
        if(AOI_hit_perc(trial)<AOI_Hit_Threshold);
            errorType=2;
        end
        
        %problem - Validity does not appear to be recorded correctly
        %contains odd values
        %tloss(trial)=sum(eyedat.Validity_Left==0)/length(eyedat.Validity_Left);
        %if( tloss(trial) < tloss_threshold)
        %    tloss(trial)
        %    errorType=3;
        %    keyboard
        %end 
        
       %sum(eyedat.Validity_Right==0)./length(eyedat.Validity_Right)
        
        percent_of_trial_with_trackloss(trial,:)=...
            [RL_eye_X_empty(trial,:) RL_eye_Y_empty(trial,:)]*100;
        
        if(errorType==0)
            %writeData;
            
            %put subject detail into numeric code
            if(strcmp( subject_summary.Group{sub_rng(1)}, 'control' ))%CN
                sub_type=0;
            else
                sub_type=1;
            end
            
            if(strcmp(  subject_summary.Age{sub_rng(1)}, 'young' ))
                sub_age=0;
            else %old
                sub_age=1;
            end
            
            
            tmp=condition{trial}; %subject_summary.Item{sub_rng(trial)};
            if(strcmp( tmp(1), 'a' )) %a for audio
                sub_presentType=0;
            else % t for text
                sub_presentType=1;
            end
            
            tmp= resp_1{trial}; %subject_summary.TargetPos1{sub_rng(trial)};
            %cul, ins, bio met
            
            if( strcmp( tmp(2), 'c'))
                sub_idiomType=0;
            elseif( strcmp( tmp(2), 'i'))
                sub_idiomType=1;
            elseif( strcmp( tmp(2), 'b'))
                sub_idiomType=2;
            elseif( strcmp( tmp(2), 'm'))
                sub_idiomType=3;
            end
            
            
            n=n+1;
            sub_details(n,:)=[sub_type sub_age sub_presentType sub_idiomType];
            subAOI(n).trial=eye_AOI_StimulusList;
            subAOI(n).norm=interp1(0:( 1/(length(eye_AOI_StimulusList)-1) ):1,eye_AOI_StimulusList,0:0.0001:1);
        end
        
        if(visualize)
            keyboard
        end
    end
    
    if(errorType==0)
        
        %take one subjects data and average over the 4 idiom types
        ave_AOI_per_idiom=[]; cntr=0;
        %for kk=0:1 %for now do not use text vs. audio
            for zz=0:3
                %for each type of idiom
                idx1=sub_details(:,4)==zz;
                %also was it text or audio?
                %idx2=sub_details(:,3)==kk; %for now do not use text vs. audio
                
                idx=find(idx1); %find(idx1&idx2); %for now do not use text vs. audio
                cntr=cntr+1;
                for ll=1:4 %types of pictures
                    
                    AOI_list=zeros(1,10001);
                    for qq=1:length(idx)
                        
                        %tmp=subAOI(idx(qq)).norm;
                        tmp=subAOI(idx(qq)).trial;
                        
                        %which AOI do we want? A,B,C,D?  A is target
                        tmp = tmp==ll;
                        %AOI_list(qq,:,ll)=tmp;
                        AOI_list(qq,1:length(tmp),ll)=tmp;
                    end
                    
                    if(~isempty(idx))
                        ave_AOI_per_idiom(cntr,:,ll)=mean(AOI_list(:,:,ll),1);
                    else
                        ave_AOI_per_idiom(cntr,:,ll)=zeros(1,10001)+nan;
                    end
                end
                
            %end %for now do not use text vs. audio
            end
        
        if(visualize)
            plot(ave_AOI_per_idiom(:,:,3)')
            ylim([0 1])
            ylabel('Normalized Time')
            xlabel('P(Fixate Target)')
            legend('cul-aud', 'ins-aud', 'bio-aud', 'met-aud',...
                'cul-txt', 'ins-txt', 'bio-txt', 'met-txt')
        end
        
        %now save these per subject with sub details
        AllSub_Details(j,:)=sub_details(1,1:2);
        AllSub_ave_AOI_per_idiom(:,:,:,j)=ave_AOI_per_idiom;
        
        fclose('all');
    end

    
    clear subjID RT  error mouse_distance_traveled mouse_num_zero_crossings ...
        mouse_peakVel mouse_num_AOIs_visited AOI_mouse_visits AOI_total_dur ...
        AOI_num_visits AOI_mean_dur exclude_trial_list RL_eye_X_empty ...
        RL_eye_Y_empty filelistCSV AOI_hit_perc percent_of_trial_with_trackloss
    

end %per subject file

% figure(101); clf;
% subplot(2,1,1)
% plot( mean(AllSub_ave_AOI_per_idiom(1:4,:,:),3)' )
% ylim([0 0.5])
% xlabel('Normalized Time')
% ylabel('P(Fixate Target)')
% legend('cul-aud', 'ins-aud', 'bio-aud', 'met-aud')
%
% subplot(2,1,2)
% plot( mean(AllSub_ave_AOI_per_idiom(5:8,:,:),3)' )
% ylim([0 0.5])
% xlabel('Normalized Time')
% ylabel('P(Fixate Target)')
% legend('cul-txt', 'ins-txt', 'bio-txt', 'met-txt')



%--generate a plot with all subject groups for each type of metaphor test--
%young_control
sub(1).list=find(AllSub_Details(:,1)==0 & AllSub_Details(:,2)==0);
%old_control=
sub(2).list=find(AllSub_Details(:,1)==0 & AllSub_Details(:,2)==1);
%young_asd=
sub(3).list=find(AllSub_Details(:,1)==1 & AllSub_Details(:,2)==0);
%old_asd=
sub(4).list=find(AllSub_Details(:,1)==1 & AllSub_Details(:,2)==1);

disp('subject types')
disp('YoungControl OldControl YoungASD OldASD')
[length(sub(1).list) length(sub(2).list) length(sub(3).list) length(sub(4).list)]

%use for text v audio
%ttl={'Cultural-Audio', 'Instructive-Audio', 'Biological-Audio', ...
%    'Novel Metaphor-Audio', 'Cultural-Text', 'Instructive-Text', ...
%    'Biological-Text', 'Novel Metaphor-Text'};

ttl={'Cultural', 'Instructive', 'Biological', ...
    'Novel Metaphor'};

picName={'A-Target'; 'B '; 'C '; 'D '};
c=50;

clr='rgcb';
t=(1:10001)/120;
cntr=0;

keyboard

for pictureNum=1:4
    for pxp=1:size(AllSub_ave_AOI_per_idiom,1)
        cntr=cntr+1;
        
        figure(151+cntr); clf; hold on
        set(gca,'FontSize',20)
        %plot first just so legend corresponds as shplot won't give proper
        %legend
        for szs=1:4
            plot(t, smooth(mean(AllSub_ave_AOI_per_idiom(pxp,:,pictureNum,sub(szs).list),4)',c),...
                clr(szs),'linewidth',3 )
        end
        legend('Control-Young','Control-Old','ASD-Young','ASD-Old')
        
        for szs=1:4
            plot(t, smooth(nanmean(AllSub_ave_AOI_per_idiom(pxp,:,pictureNum,sub(szs).list),4)',c),...
                clr(szs),'linewidth',3 )
            sem=smooth(nanstd(AllSub_ave_AOI_per_idiom(pxp,:,pictureNum,sub(szs).list),...
                [],4)./sqrt(length(sub(szs).list)),c)';
            shplot( t, ...
                smooth(nanmean(AllSub_ave_AOI_per_idiom(pxp,:,pictureNum,sub(szs).list),4)',c), ...
                sem, clr(szs) );
        end
        
        ylim([0 0.55])
        xlim([0 15])
        %xlabel('Normalized Time')
        xlabel('Time (s)')
        
        ylabel('P(Fixate Target)')
        title([picName(pictureNum) ttl{pxp}]);
        saveas(gcf,[picName{pictureNum} ttl{pxp} '.jpg'])
        ax = gca;
    end
end

