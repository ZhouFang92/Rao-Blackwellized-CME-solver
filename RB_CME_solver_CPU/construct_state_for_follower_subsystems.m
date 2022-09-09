function followerX_correction=construct_state_for_follower_subsystems...
                              (num_species,truncated_size,...
                               follower_subsystems)
                           
% This function returns the state for each follower subsystems. 

followerX_correction=cell(1,size(follower_subsystems,2));

for i=1:size(follower_subsystems,2) % consider the i-th follower subsystem
     truncated_size_follower_subsystem=truncated_size(follower_subsystems{i});
     size_of_the_matrix=prod(truncated_size_follower_subsystem);
     species_index=follower_subsystems{i};
     
     followerX_correction{i}=zeros(num_species,size_of_the_matrix);
    
     
     for LI=1:size_of_the_matrix % LI is the linear index
         MI=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI);
         if size(MI,1)==1
            MI=MI'; % make MI a column vector
         end
         followerX_correction{i}(species_index,LI)=MI-1;
     end
     
end