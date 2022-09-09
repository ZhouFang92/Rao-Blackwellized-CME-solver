function d=Distribution_RB(X,P,truncated_size,follower_subsystems,species_decomposition_LF)

% this function summarize the distribution of the RB method

leader_system_species=find(species_decomposition_LF==0);
follower_system_species=find(species_decomposition_LF~=0);
truncated_size_follower=truncated_size(follower_system_species);

sample_size=size(X,2); % number of particles
d=zeros(truncated_size'); % initialize the distribution
if size(truncated_size,1)==1
   d=zeros(truncated_size,1); 
end

parfor j=1:sample_size
    dj=zeros(truncated_size'); % the distribution by the j-th sample
    if size(truncated_size,1)==1
       dj=zeros(truncated_size,1); 
    end
    MI=zeros(size(X,1),1);
    MI(leader_system_species)=X(leader_system_species,j)+1;
    % calculate dj
    for LI_follower=1:prod(truncated_size_follower)
         MI_follower=Linear_index_2_matrix_index(truncated_size_follower,LI_follower);
         if size(MI_follower,1)==1
            MI_follower=MI_follower'; % make MI_follower a column vector
         end
         MI(follower_system_species)=MI_follower;
         LI=Matrix_index_2_linear_index(truncated_size,MI); % the linear index of dj

         value=1;
         for i=1:size(follower_subsystems,2)
             truncated_size_follower_subsystem=truncated_size(follower_subsystems{i});
             follower_subsystem_species=follower_subsystems{i};
             MI_sub=MI(follower_subsystem_species);
             LI_sub=Matrix_index_2_linear_index(truncated_size_follower_subsystem,MI_sub);
             value=value*P{i,j}(LI_sub);
         end
         dj(LI)=value;
    end % dj is created
    d=d+dj/sample_size;
end