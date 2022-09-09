function A_update_or_not=update_the_evolution_matrix_or_not(num_species,...
                            propensity,reaction_decomposition_F,...
                            follower_subsystems,Leader_subsystems,...
                            species_decomposition_LF)
                        
% this algorithm returns a matrix whose i,j elements tell if the evolution
% of the i-th follower-subsystem needs to be changed if we spot a change of
% the j-th leader subsystem. 

% 1 means that the matrix needs to be updated, otherwise 0;

A_update_or_not=zeros(size(follower_subsystems,2),size(Leader_subsystems,2));

leader_system_species=find(species_decomposition_LF==0);
v=species_influencing_reactions(num_species,propensity);

for i=1:size(follower_subsystems,2)
   %identify reactions influecing this subsystem
   reactions_involved=find(reaction_decomposition_F==i);
    
   for j=1:size(Leader_subsystems,2) 
       % check if j influences those reactions influcing i-th subsystems
       
       X=zeros(num_species,1);
       X(leader_system_species)=Leader_subsystems(:,j);
       species_changed=find(X~=0); % find the species whose dynamics are changed by the j-th leader subsystem
       
       %check if the species changed affects the reactions influcing the
       %follower subsystem
       for index=1:size(reactions_involved,2)
           if min(size(intersect(species_changed,v{reactions_involved(index)})))>0
               A_update_or_not(i,j)=1; % 
               break;
           end
       end
       
   end
end


