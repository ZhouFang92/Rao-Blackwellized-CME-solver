function [species_decomposition_LF,follower_subsystems,reaction_decomposition_LF,...
           reaction_decomposition_F,Size_of_truncated_space]...
          =conditional_independence_decomposition(num_species,reactant_vector,...
           product_vector,propensity,h,truncated_size,threshold_follower_subsystem)
       
% Initialization

% Preparation (show the interactions)
v_SpeciesInfluencingReactions=species_influencing_reactions(num_species,propensity);
v_SpeciesInfluencedByReactions=species_influenced_by_reactions(reactant_vector,product_vector);
v_SpeciesMeasuredByObservations=species_measured_by_observation_channels(num_species,h);

species_decomposition_LF=ones(1,num_species); % all species belongs to the follower system
optimal_group=cell(1,num_species);
optimal_group{1}=[1:1:num_species];   % only one follower subsystem and it contains all the species.
SL=prod(truncated_size); % size of the largest follower subsystem.
Reactions_in_leader_subsystems=[{}];

% Loop
while SL>threshold_follower_subsystem
    % save the decomposition obtained from the previous step
    old_species_decomposition_LF=species_decomposition_LF; 
    for i=1:num_species
        % First level decomposition
        if species_decomposition_LF(i)==1
            try_decomposition_LF=old_species_decomposition_LF;
            try_decomposition_LF(i)=0; %Flip the species
        else
            continue;
        end
        
        % second level decomposition
        [group,v_Reactions_in_leader_subsystems]...
         =second_level_decomposition(num_species,try_decomposition_LF,...
         reactant_vector,product_vector,h,...
         v_SpeciesInfluencingReactions,v_SpeciesInfluencedByReactions,...
         v_SpeciesMeasuredByObservations);
    
        % evaluate the size
        for j=1:num_species
         if  size(group{j},2)==0
             try_size_follower_subsystems(j)=0;
         else
             try_size_follower_subsystems(j)=1;
         end
         for i=1:size(group{j},2)
             species_index=group{j}(i);
             try_size_follower_subsystems(j)=try_size_follower_subsystems(j)*truncated_size(species_index);
         end
        end
        
        
        % update the optimal decomposition
        if max(try_size_follower_subsystems)<=SL
            SL=max(try_size_follower_subsystems);
            species_decomposition_LF=try_decomposition_LF;
            optimal_group=group;
            Reactions_in_leader_subsystems=v_Reactions_in_leader_subsystems;
        end
    
    
    
    end
end


% Determine the decomposition of the reactions
num_reactions=size(reactant_vector,2);
[follower_subsystems,reaction_decomposition_LF,...
           reaction_decomposition_F,Size_of_truncated_space]...
          =prepare_output_decomposition(num_reactions,species_decomposition_LF,...
           optimal_group,truncated_size,Reactions_in_leader_subsystems,...
           v_SpeciesInfluencingReactions,v_SpeciesInfluencedByReactions);