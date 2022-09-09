function [t,X]=next_reaction_method_single_trajectory...
               (X0,k,phi,FT,reactant_vector,product_vector)


% Initialization
reaction_vector=product_vector-reactant_vector;
num_reactions=size(reaction_vector,2);
X(:,1)=X0;
t(1)=0;
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
num_step=1;

% simulation
while t(num_step)+deltaT<FT % update while the final time is not reached
    t(num_step+1)=t(num_step)+deltaT; %update the time
    DX=reaction_vector(:,mu);
    X(:,num_step+1)=X(:,num_step)+DX; %update the state
    num_step=num_step+1;
    T=T+deltaT*a; % update the used time for each reaction
    r(mu)=rand;
    P(mu)=P(mu)+log(1/r(mu)); % update the next reaction time for the mu-th reaction
    for j=1:num_reactions
        a(j)=propensity_function(j,X(:,num_step),k,phi); % calculate the propensity
    end  
    if max(a)>0
       DeltaT=(P-T)./a;
       [deltaT,mu]=min(DeltaT);
    else
       deltaT=FT; % to avoid the all zero propensities.
    end
end

t(num_step+1)=FT;
X(:,num_step+1)=X(:,num_step);
