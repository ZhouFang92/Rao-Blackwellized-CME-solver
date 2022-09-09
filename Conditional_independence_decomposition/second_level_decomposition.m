function [group,v_Reactions_in_leader_subsystems]...
         =second_level_decomposition(num_species,try_decomposition_LF,...
         reactant_vector,product_vector,h,...
         v_SpeciesInfluencingReactions,v_SpeciesInfluencedByReactions,...
         v_SpeciesMeasuredByObservations)


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