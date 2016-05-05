function AOI_segments=AOI_segmentation(AOI)

start=1;
label=AOI(1);
t_i=1; t_f=[]; finalLabel=[];

for i=2:length(AOI)
    
    if(start)
        if(AOI(i)~=label)
            t_f(end+1)=i-1;
            finalLabel(end+1)=AOI(i-1);
            start=0;
            label=AOI(i);
        elseif(i==length(AOI))
            t_f(end+1)=i;
            finalLabel(end+1)=AOI(i);
        end
    elseif(~start & i~=length(AOI) )
        label=AOI(i);
        t_i(end+1)=i;
        start=1;
    elseif(i==length(AOI))
        label=AOI(i);
        t_i(end+1)=i;
        t_f(end+1)=i;
        finalLabel(end+1)=AOI(i);
        start=1;
    end
end

AOI_segments=[t_i' t_f' finalLabel'];

return