function observation_decomposition_F=...
         observations_belonging_to_each_follower_subsystem...
         (num_species,h,follower_subsystems)
     
% This algorithm returns a cell vector indicating which follower subsystem
% it belongs to. 

% species involved in each observation
v=species_measured_by_observation_channels(num_species,h);
observation_decomposition_F=zeros(1,max(size(h)));

for j=1:size(v,2)
     for i=1:size(follower_subsystems,2)
          if all(size(intersect(follower_subsystems{i},v{j})))>0  
             observation_decomposition_F(j)=i;
          end 
     end
end