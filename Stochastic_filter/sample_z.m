function X=sample_z(X,P,follower_subsystems,followerX_correction)

% this function samples z from the given distribution P

for j=1:size(X,2)
    for k=1:size(follower_subsystems,2) % consider z_k 
        species_follower_subsystem=follower_subsystems{k};
        d=cumsum(P{k,j});    
        a=rand();
        index=find(a<d,1,'first');
        % update z_k
        X(species_follower_subsystem,j)...
            =followerX_correction{k}(species_follower_subsystem,index);
    end
end