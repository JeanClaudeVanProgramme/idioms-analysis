function [distance_traveled, zero_crossings, num_AOIs_visited,...
    peakVel,AOI_segments,errorType] = mouseMetrics(mouse,init_time)

distance_traveled=[]; zero_crossings=[]; num_AOIs_visited=[];
peakVel=[]; AOI_segments=[]; errorType=[];
    
errorType=0;

vel_thresh=0.2; %pix/ms
X=mouse.X;
Y=mouse.Y;
T=mouse.T;
%AOI=mouse.AOI_PositionList;
AOI=mouse.AOI_StimulusList;



%total distance traveled'
X_vel=(X(2:end)-X(1:end-1))./(T(2:end)-T(1:end-1));
Y_vel=(Y(2:end)-Y(1:end-1))./(T(2:end)-T(1:end-1));

dist=sqrt( (X(2:end)-X(1:end-1)).^2 + ...
    (Y(2:end)-Y(1:end-1)).^2 );
distance_traveled=sum(dist);

%find zero crossings
vel = dist./( (T(2:end)-T(1:end-1))/1000) ; %pix/sec


[~,min_idx]=min(abs(mouse.T-init_time));

v_idx=min_idx:length(vel);

if(isempty(v_idx))
    errorType=3; %wrong timing for init?
    return;
end


Vel_segments=AOI_segmentation(vel(v_idx)>vel_thresh);
zero_crossings=length(find(Vel_segments(:,3)==0));

AOI_segments=AOI_segmentation(AOI);
num_AOIs_visited=size(AOI_segments,1);

peakVel=max(vel);

if(0)
    figure; hold on
    %plot(T(2:end),vel)
    %plot( T(v_idx-1:end-2),vel(v_idx:end),'r');
    plot(vel(v_idx))
    hold on; plot(Vel_segments(:,1), 0.5,'g*')
    hold on; plot(Vel_segments(:,2), 0.5,'r*')
end

%keyboard
return
