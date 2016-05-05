% function ans=smooth(data,smoothness)

% Smooth 'data' with a gaussian of stdev of 'smoothness'

% note: data must be a column vector

% 

% pilardi 5/02



function a=smooth(data,smoothness)

kernel=gauss(-5*smoothness:5*smoothness,0,smoothness)';


a=conv(data,kernel);
ksize = size(kernel,1);
datasize = length(a); %size(a,1);
a = a(floor(ksize/2):datasize-(floor(ksize/2)+1));


function ans=gauss(x,m,sd)

ans=exp(-(x-m).^2/(2*sd^2))/(sd*sqrt(2*pi));


function ans = myconv(data,kernel)
  