function v=species_influenced_by_reactions(reactant_vector,product_vector)

% This function returns a cell array the j-th of which shows the species
% influenced by the j-th reaction. 

reaction_vector=product_vector-reactant_vector; % reaction_vectors
v{1}=[];
for j=1:size(reaction_vector,2)
    v{j}=find(reaction_vector(:,j)~=0)';
end