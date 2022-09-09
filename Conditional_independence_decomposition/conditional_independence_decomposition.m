function  [species_decomposition_LF,follower_subsystems,reaction_decomposition_LF,...
           reaction_decomposition_F,Size_of_truncated_space]...
          =conditional_independence_decomposition(num_species,reactant_vector,...
           product_vector,propensity,h,truncated_size,threshold_follower_subsystem)
       
% This function decompose the whole system.

% Preparation
v_SpeciesInfluencingReactions=species_influencing_reactions(num_species,propensity);
v_SpeciesInfluencedByReactions=species_influenced_by_reactions(reactant_vector,product_vector);
v_SpeciesMeasuredByObservations=species_measured_by_observation_channels(num_species,h);

% Decompose the network
species_decomposition_LF=zeros(1,num_species); % all species belongs to the leader system
size_follower_system=0;
optimal_group=cell(1,num_species);
Reactions_in_leader_subsystems=[{}];

for try_decomposition_index=1:2^num_species-1
    
    % First level decomposition
    try_decomposition_LF=num2binary(num_species,try_decomposition_index);
    
    % Second level decomposition
    for i=1:num_species
        if try_decomposition_LF(i)==1 %check if the i-th species belongs to the follower
          group{i}=[i];
        else
          group{i}=[];
        end
    end
    
    %merge interactive groups
    for j=1:size(reactant_vector,2)
        %find the species influencing or influenced by the j-th reaction
        species_involved= [v_SpeciesInfluencingReactions{j},v_SpeciesInfluencedByReactions{j}];
        % find the groups influencing or influenced by this reaction 
        Index=[];
        for i=1:num_species
            if all(size(intersect(species_involved,group{i})))>0
               Index=[Index,i]; 
            end
        end
        %merge these groups
        for i=2:size(Index,2)
            group{Index(1)}=[group{Index(1)},group{Index(i)}];
            group{Index(i)}=[];
        end
    end
    
    % merge those groups that influence the same leader subsystems
    [Leader_subsystems,v_Reactions_in_leader_subsystems]...
               =determine_leader_level_reactions(reactant_vector,...
                                product_vector,...
                                try_decomposition_LF);
                            % identify leader_systems and those species it
                            % influences. 
     for j=1:size(Leader_subsystems,2)
         %find the species influencing j-th leader subsystem
         species_index=[];
         for i=1:size(v_Reactions_in_leader_subsystems{j},2)
             reaction_index=v_Reactions_in_leader_subsystems{j}(i);
             species_index=[species_index,v_SpeciesInfluencingReactions{reaction_index},...
                            v_SpeciesInfluencedByReactions{reaction_index}];
         end
         % find the groups influencing this leader subsystem 
         Index=[];
         for i=1:num_species
            if all(size(intersect(species_index,group{i})))>0
               Index=[Index,i]; 
            end
         end
         %merge these groups
         for i=2:size(Index,2)
            group{Index(1)}=[group{Index(1)},group{Index(i)}];
            group{Index(i)}=[];
         end    
     end
     
     %Merge those groups measured by the same observations
     for j=1:size(h,2)
         %Identify the groups influencing this reaction channel
         Index=[];
         for i=1:num_species
            if all(size(intersect(v_SpeciesMeasuredByObservations{j},group{i})))>0
               Index=[Index,i]; 
            end
         end
         %merge these groups
         for i=2:size(Index,2)
            group{Index(1)}=[group{Index(1)},group{Index(i)}];
            group{Index(i)}=[];
         end   
     end
     
     % calculate the size
     try_size_follower_system=1;
     for i=1:num_species
         if try_decomposition_LF(i)==1
             try_size_follower_system=try_size_follower_system*truncated_size(i);
         end
     end
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
     if  try_size_follower_system>size_follower_system && ...
          all(try_size_follower_subsystems<=threshold_follower_subsystem)
         
         size_follower_system=try_size_follower_system;
         species_decomposition_LF=try_decomposition_LF;
         optimal_group=group;
         Reactions_in_leader_subsystems=v_Reactions_in_leader_subsystems;
     end
    
end


% Determine the decomposition of the reactions
num_reactions=size(reactant_vector,2);
[follower_subsystems,reaction_decomposition_LF,...
           reaction_decomposition_F,Size_of_truncated_space]...
          =prepare_output_decomposition(num_reactions,species_decomposition_LF,...
           optimal_group,truncated_size,Reactions_in_leader_subsystems,...
           v_SpeciesInfluencingReactions,v_SpeciesInfluencedByReactions);