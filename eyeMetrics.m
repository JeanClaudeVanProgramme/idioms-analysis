function [AOI_segments,AOI_total_dur,AOI_num_visits,...
    AOI_mean_dur] = eyeMetrics(eyeBinocXY_AOI,smp_rate)


%minumum of 100ms is about 12 frames at 120hz
t_thresh=12;
AOI_segments=AOI_segmentation(eyeBinocXY_AOI);

td=AOI_segments(:,2)-AOI_segments(:,1);
AOI_segments=AOI_segments(find(td>=t_thresh),:);
AOI_sequence=AOI_segments(:,3)';

%find total duration
dur=(AOI_segments(:,2)-AOI_segments(:,1)+1)*smp_rate;

%step thru types of aoi 0, 1-4
for i=1:5
    AOI_total_dur(i)=sum(dur(find(AOI_sequence==(i-1))));
    AOI_num_visits(i)=sum(AOI_sequence==(i-1));
end

AOI_mean_dur=AOI_total_dur./AOI_num_visits;
%number of visits
%find mean duration
return