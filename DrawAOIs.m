function h = DrawAOIs(AOI,xy_res)
%drawAOIs


for i=1:size(AOI,1)
    rectangle('Position', AOI(i,:))
end

rectangle('Position',[0 0 xy_res],'linestyle','--')

%
% %top right
% rectangle('Position', topright)
%
% %bottom left
% rectangle('Position', bottomleft)
%
% %bottom right
% rectangle('Position', bottomright)

return
