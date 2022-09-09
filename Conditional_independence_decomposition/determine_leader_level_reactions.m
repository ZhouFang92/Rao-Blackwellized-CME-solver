function   [Leader_subsystems,v]=determine_leader_level_reactions(reactant_vector,...
                                product_vector,...
                                species_decomposition_LF)
          
% This function returns the leader-subsystems and a cell array the j-th of
% which indicates the leader-level reactions corresponding to the j-th
% leader subsystem.

%species_decomposition_LF: 0 for leader level species, and 1 for follower
%level species.

reaction_vector=product_vector-reactant_vector; % reaction_vectors
index=find(species_decomposition_LF==0);
if size(index,2)==0
   Leader_subsystems=[];
   v=[{}];
end
reaction_vector_Leader=reaction_vector(index,:);
Leader_subsystems=(unique(reaction_vector_Leader','rows'))'; % find distinct vectors
Leader_subsystems(:,all(~Leader_subsystems,1))=[];           % Delete the zero vector

for j=1:size(Leader_subsystems,2) % traverse all the leader subsystems
    v{j}=[];
    for i=1:size(reaction_vector_Leader,2)
        if reaction_vector_Leader(:,i)==Leader_subsystems(:,j)
            v{j}=[v{j},i];
        end
    end
end

