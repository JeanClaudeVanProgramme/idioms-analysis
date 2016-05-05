function [AOI_PositionList, AOI_StimulusList]= IsInAOI(x,y,AOI, AOI_label)
% inside = IsInRect(x,y,rect)
%
% Decide if location x,y is inside of the passed rect.
%
% 3/5/97  dhb  Wrote it., modified by bts 3-5-2015

%AOI label is from sorting the type of response
%entry slots 1:4 correspond to ABCD, numbers in the slot represent where
%on screen ABCD were located.
AOI_PositionList=zeros(length(x),1);
AOI_StimulusList=zeros(length(x),1);

%using matlab rect x,y,w,h
for i=1:size(AOI,1)
    RectLeft=AOI(i,1);
    RectRight=AOI(i,1)+AOI(i,3);
    RectBottom=AOI(i,2);
    RectTop=AOI(i,2)+AOI(i,4);
    
    inside(:,i)= (x >= RectLeft & x <= RectRight & ...
        y <= RectTop & y >= RectBottom);
    
    AOI_PositionList=AOI_PositionList + inside(:,i)*i; %give numeric location label
    AOI_StimulusList=AOI_StimulusList + inside(:,i)*AOI_label(i); %give numeric location label
end

return
