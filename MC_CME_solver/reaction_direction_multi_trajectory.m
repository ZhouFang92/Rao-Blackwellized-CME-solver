function zeta=reaction_direction_multi_trajectory(reaction_vector,mu)

num_species=size(reaction_vector,1);
%num_reactions=size(reaction_vector,2);
%zeta=zeros(num_species,size(mu,2),'single','gpuArray');


%for j=1:num_reactions
%    zeta= zeta+reaction_vector(:,j)*(mu==j); 
%end

zeta=reaction_vector(:,mu);
zeta=single(gpuArray(zeta));