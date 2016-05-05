filelistArray=dir([pwd datadir fdir '\*.mt']);
[~,idx] = sort([filelistArray.datenum]);

for zz= 1:length(idx)
    filelistMT{zz}=filelistArray(idx(zz)).name;
    if(zz==1)
        [subjID_1,stim_1, order_1,condition_1,resp_1_1, resp_2_1, resp_3_1,...
            resp_4_1, response_1,error_1,resp_num_1, RT_1, init_time_1,...
            distractor_1, ideal_y_int_1, maxdev_1, real_time_1, comments_1,...
            timestamps_1,mouse1 ]=ReadMouseTrackerData([pwd datadir fdir '\' filelistMT{zz}]);
    elseif(zz==2)
        [subjID_2,stim_2, order_2,condition_2,resp_1_2, resp_2_2, resp_3_2,...
            resp_4_2, response_2,error_2,resp_num_2, RT_2, init_time_2,...
            distractor_2, ideal_y_int_2, maxdev_2, real_time_2, comments_2,...
            timestamps_2,mouse2 ]=ReadMouseTrackerData([pwd datadir fdir  '\' filelistMT{zz}]);
    end
end



subjID=[subjID_1; subjID_2];
stim=[stim_1; stim_2];
order=[order_1; order_2];
condition=[condition_1; condition_2];
resp_1=[resp_1_1; resp_1_2];
resp_2=[resp_2_1; resp_2_2];
resp_3=[resp_3_1; resp_3_2];
resp_4=[resp_4_1; resp_4_2];
response=[response_1; response_2];
error=[error_1; error_2];
resp_num=[resp_num_1; resp_num_2];
RT=[RT_1; RT_2];
init_time=[init_time_1; init_time_2];
distractor=[distractor_1; distractor_2];
%ideal_y_int=[ideal_y_int_1; ideal_y_int_2];
%maxdev=[maxdev_1; maxdev_2];
%real_time=[real_time_1; real_time_2];
comments=[comments_1; comments_2];
%timestamps=[timestamps_2; timestamps_2];

mouse=mouse1;

for zz=1:length(mouse2)
    mouse(length(mouse2)+zz).X=mouse2(zz).X;
    mouse(length(mouse2)+zz).Y=mouse2(zz).Y;
    mouse(length(mouse2)+zz).T=mouse2(zz).T;
    
end

clear subjID_1 stim_1  order_1 condition_1 resp_1_1  resp_2_1  resp_3_1 ...
    resp_4_1  response_1 error_1 resp_num_1  RT_1  init_time_1 ...
    distractor_1  ideal_y_int_1  maxdev_1  real_time_1  comments_1 ...
    timestamps_1  subjID_2 stim_2  order_2 condition_2 resp_1_2  resp_2_2  resp_3_2 ...
    resp_4_2  response_2 error_2 resp_num_2  RT_2  init_time_2 ...
    distractor_2  ideal_y_int_2  maxdev_2  real_time_2  comments_2 ...
    timestamps_2 zz mouse1 mouse2
