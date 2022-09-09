function [w,p]=Correction_step_RBPF(w,X,p,y,follower_subsystems,...
                                    h,sigma,...
                                    observation_decomposition_F,...
                                    followerX_correction)

% this algorithm corresponds to the correction step in the RB-PF. 

% Compute the likelihood function
L=cell(1,size(follower_subsystems,2)+1); % the last one corresponds to the leader level likelihood function
for k=1:size(L,2)-1
    % construct the state matrix
    species_follower_subsystem=follower_subsystems{k};
    followerX=X.*ones(size(X,1),size(followerX_correction{k},2));
    followerX(species_follower_subsystem,:)...
        =followerX_correction{k}(species_follower_subsystem,:);
    
    total_L=likelihood_matrix_RB(y,followerX,h,sigma);
    observation_follower_subsystem=find(observation_decomposition_F==k);
    if min(size(observation_follower_subsystem))>0
       L{k}=prod(total_L(observation_follower_subsystem,:),1);
    else
       L{k}=ones(1,size(followerX,2)); 
    end
end
% the leader level observations
total_L=likelihood_matrix_RB(y,X,h,sigma);
observation_follower_subsystem=find(observation_decomposition_F==0);
if min(size(observation_follower_subsystem))>0
   L{end}=prod(total_L(observation_follower_subsystem,:),1);
else
   L{end}=1; 
end

% update weights and p
for k=1:size(L,2)-1
    p_temp=(L{k})'.*p{k};
    p_sum=sum(p_temp);
    if p_sum==0
        p{k}=0*p_temp;
    else
        p{k}=p_temp/p_sum;
    end
    w=w*p_sum;
end
w=w*L{end};
