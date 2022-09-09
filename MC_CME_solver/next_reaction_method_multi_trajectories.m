function [TX,X]=next_reaction_method_multi_trajectories...
                (k,phi,FT,reactant_vector,product_vector,...
                 IniDistribution_Mean,IniDistribution_Variance,...
                 sample_size)
             
% This function uses the next reaction method to generate
% multi_trajectories for the network. 


% Initialization
%num_species=size(reactant_vector,1); % the number of species
num_reactions=size(reactant_vector,2);
reaction_vector=product_vector-reactant_vector;

X=Initial_samples(IniDistribution_Mean,IniDistribution_Variance,...
                  sample_size);
X=single(gpuArray(X));
TX=zeros(1,sample_size,'single','gpuArray');     % the time vector

saved_number=0; % number of samples saved.
left_particles=gpuArray.linspace(1,sample_size,sample_size); % the index of the particles left
left_particles=single(left_particles);
tempX=X;
tempTX=TX;

% the time vector for the Poisson processes (associated with the random time changes)
T=zeros(num_reactions,sample_size,'single','gpuArray');     
% the jumping time of the Poisson processes
P=log(1./rand(num_reactions,sample_size,'single','gpuArray'));     

% updates

while saved_number<sample_size
      % deal with the left particles
      
      a=propensity_function_multi_trajectories(tempX,k,phi); % calculate the propensity
      delta_t=(P-T)./a;
      [delta,mu]=min(delta_t);    % find the next reaction time. 
      mu=single(mu);
      
      %update the prcoesses
      zeta=reaction_direction_multi_trajectory(reaction_vector,mu);
      try_tempTX=tempTX+delta;
      check_t=(try_tempTX<=FT); % whether the next jump time exceeds the final time
      tempX=tempX+zeta.*check_t;
      tempTX=min(try_tempTX,FT);
      
      %update T
      T=T+a.*delta;
      
      %update P
      r=log(1./rand(1,size(P,2),'single','gpuArray'));
      for j=1:size(P,1)
          P(j,:)= P(j,:)+r.*(mu==j);  % update the P value associated with the fired reaction
      end
      
      %save particles
      particles_to_save=find(tempTX>=FT); 
      particles_to_save=single(particles_to_save);
      if size(particles_to_save,2)>20000 || size(particles_to_save,2)==size(left_particles,2)
         X(:,saved_number+1:saved_number+size(particles_to_save,2))=tempX(:,particles_to_save);
         TX(:,saved_number+1:saved_number+size(particles_to_save,2))=tempTX(:,particles_to_save);
         saved_number=saved_number+size(particles_to_save,2);
         left_particles=find(tempTX<FT);
         left_particles=single(left_particles);
         %size(left_particles,2)
         tempX=tempX(:,left_particles);
         tempTX=tempTX(:,left_particles);
         T=T(:,left_particles);
         P=P(:,left_particles);
      end

end



% finished