function plot_RB_result(X,P,truncated_size,follower_subsystems,...
                        species_decomposition_LF)

sample_size=size(X,2);
num_species=size(X,1);
                    
% Initialization
d_marginal=cell(num_species,num_species); % the marginal distribution
for i=1:num_species
    for j=1:i
        d_marginal{i,j}=zeros(truncated_size(i),truncated_size(j));
    end
end

% calculate the marginal distribution
for i=1:num_species
     for j=1:i
          %Both are leader-level species
          if species_decomposition_LF(i)==0 && species_decomposition_LF(j)==0
              for index=1:sample_size
                  if X(i,index)+1<=truncated_size(i) && X(j,index)+1<=truncated_size(j)
                  d_marginal{i,j}(X(i,index)+1,X(j,index)+1)...
                   =d_marginal{i,j}(X(i,index)+1,X(j,index)+1)+1/sample_size; 
                  continue;
                  end
              end
          end
          
          % Both are follower-level species
          if species_decomposition_LF(i)~=0 && species_decomposition_LF(j)~=0
              d_marginal{i,j}=plot_RB_both_follower(i,j,X,P,truncated_size,follower_subsystems,...
                        species_decomposition_LF);
              continue;
          end
              
          % i is the leader, j is the follower
          if species_decomposition_LF(i)==0 && species_decomposition_LF(j)~=0
              d_marginal{i,j}=plot_RB_leader_plus_follower(i,j,X,P,truncated_size,follower_subsystems,...
                        species_decomposition_LF);
              continue;
          end
          
          % i is the follower, j is the leader
          if species_decomposition_LF(i)~=0 && species_decomposition_LF(j)==0
               d_marginal{i,j}=plot_RB_follower_plus_leader(i,j,X,P,truncated_size,follower_subsystems,...
                         species_decomposition_LF); 
              continue;
          end
     end
end


% plot
for i=1:num_species
    for j=1:i
        subplot(num_species,num_species,(i-1)*num_species+j);
        if i==j 
            % plot the marginal distribution
            value=diag(d_marginal{i,j});
            bar_plot=bar([0:1:size(value,1)-1],value);
            set(bar_plot,'edgecolor','blue');
        else
            % plot the joint distribution
            h=heatmap(d_marginal{i,j},'GridVisible','off');
            for index=1:size(d_marginal{i,j},2)
               h.XDisplayLabels{index}='';
            end
            for index=1:size(d_marginal{i,j},1)
               h.YDisplayLabels{index}='';
            end
            h.XDisplayLabels{1}='0';
            h.XDisplayLabels{size(d_marginal{i,j},2)}=string(size(d_marginal{i,j},2)-1);
            h.YDisplayLabels{1}='0';
            h.YDisplayLabels{size(d_marginal{i,j},1)}=string(size(d_marginal{i,j},1)-1);          
        end
    end
end