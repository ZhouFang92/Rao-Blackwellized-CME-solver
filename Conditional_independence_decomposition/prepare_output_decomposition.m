function  [follower_subsystems,reaction_decomposition_LF,...
           reaction_decomposition_F,Size_of_truncated_space]...
          =prepare_output_decomposition(num_reactions,species_decomposition_LF,...
           optimal_group,truncated_size,Reactions_in_leader_subsystems,...
           v_SpeciesInfluencingReactions,v_SpeciesInfluencedByReactions)
      
% species decomposition
num_species=size(species_decomposition_LF,2);
if max(species_decomposition_LF)==0
   follower_subsystems=[{}];
else
   index=1;
   for i=1:num_species
      if size(optimal_group{i},2)>0
         follower_subsystems{index}=optimal_group{i};
         index=index+1;
      end
   end
end

% Size_of_truncated_space
if size(follower_subsystems,2)==0
    Size_of_truncated_space=[];
else
    for i=1:size(follower_subsystems,2)
        TS=prod(truncated_size(follower_subsystems{i}));
        Size_of_truncated_space(i)=TS;
    end
end

% reaction decomposition
reaction_decomposition_LF=ones(1,num_reactions);
for j=1:size(Reactions_in_leader_subsystems,2)
    for i=1:size(Reactions_in_leader_subsystems{j},2)
        reaction_index=Reactions_in_leader_subsystems{j}(i);
        reaction_decomposition_LF(reaction_index)=0;
    end
end

% reaction decompostion follower subsystems
reaction_decomposition_F=zeros(1,num_reactions);
if size(follower_subsystems,2)>0 % there is follower systems
   for j=1:num_reactions
       species_involved= [v_SpeciesInfluencingReactions{j},v_SpeciesInfluencedByReactions{j}];
       for i=1:size(follower_subsystems,2)
           if all(size(intersect(species_involved,follower_subsystems{i})))>0
              reaction_decomposition_F(j)=i;
              break;
           end
       end
   end
end



