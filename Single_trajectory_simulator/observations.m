function [TY,Y]=observations(t,X,h,sigma,observation_period)

FT=t(size(t,2));
TY=[observation_period:observation_period:FT];

if size(sigma,2)==0
   Y=zeros(1,size(TY,2));
   return;
end

for i=1:size(TY,2)
    index=find(t<=i*observation_period,1,'last');
    Y(:,i)=h_function(h,X(:,index))+sigma'.*randn(size(sigma,2),1);
end