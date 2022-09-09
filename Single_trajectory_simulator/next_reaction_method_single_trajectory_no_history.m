function [t,X]=next_reaction_method_single_trajectory_no_history...
               (X0,k,phi,FT,reactant_vector,product_vector)


% Initialization
reaction_vector=product_vector-reactant_vector;
num_reactions=size(reaction_vector,2);
X=X0;
t=0;
T=zeros(1,num_reactions);
r=rand(1,num_reactions);
P=log(1./r);
a=zeros(1,num_reactions);
for j=1:num_reactions
    a(j)=propensity_function(j,X(:,1),k,phi);
end
if max(a)>0
  DeltaT=(P-T)./a;
  [deltaT,mu]=min(DeltaT);
else
  deltaT=FT; % to avoid the all zero propensities.
end

% simulation
while t+deltaT<FT % update while the final time is not reached
    t=t+deltaT; %update the time
    DX=reaction_vector(:,mu);
    X=X+DX; %update the state
    T=T+deltaT*a; % update the used time for each reaction
    r(mu)=rand;
    P(mu)=P(mu)+log(1/r(mu)); % update the next reaction time for the mu-th reaction
    for j=1:num_reactions
        a(j)=propensity_function(j,X,k,phi); % calculate the propensity
    end  
    if max(a)>0
       DeltaT=(P-T)./a;
       [deltaT,mu]=min(DeltaT);
    else
       deltaT=FT; % to avoid the all zero propensities.
    end
end

t=FT;