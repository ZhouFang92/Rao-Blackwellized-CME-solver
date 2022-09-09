function s=translate_propensities_for_multi_trajectory_simulation...
           (ps)

% This function translates the propensities into the format for
% multi-trajectory simulation. Here, ps is the propensity string. 

% X(i)==>X(i,:), and  *==>.*


ps=insertBefore(ps,"*","."); % *==>.*
ps=insertBefore(ps,"/","."); % /==>./
ps=insertBefore(ps,"^","."); % /==>./

ps=split(ps,"X("); % split at X(
index=strfind(ps,")"); % find the position of )
for i=2:size(ps,1)
    ps(i)=insertBefore(ps(i),index{i}(1),',:'); % it becomes i,:)
end
ps=join(ps,"X(");


% if the propensity is constant
if contains(ps,"X")==0
    %ps=strcat(ps,".*one_vector"); 
end


s=ps;