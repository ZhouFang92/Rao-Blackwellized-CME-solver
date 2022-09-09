function plot_result(d_MC,MC_samples)

% plot the result of the Monte Carlo method

% Initialization
if size(d_MC,1)>0
    truncated_size=(size(d_MC))';
    if size(d_MC,2)==1   % one dimensional distribution
       truncated_size=size(d_MC,1); 
    end
else
    truncated_size=max(MC_samples')+1;
    truncated_size=truncated_size';
end
num_species=size(truncated_size,1);
d_MC=gather(d_MC);
X=gather(MC_samples);

d_marginal=cell(num_species,num_species); % the marginal distribution
for i=1:num_species
    for j=1:i
        d_marginal{i,j}=zeros(truncated_size(i),truncated_size(j));
    end
end

% compute the marginal distribution
if size(d_MC,2)>0
    for index=1:prod(truncated_size)
        value=d_MC(index);
        MI=Linear_index_2_matrix_index(truncated_size,index); % get the matrix index
        for i=1:num_species
           for j=1:i
               d_marginal{i,j}(MI(i),MI(j))=d_marginal{i,j}(MI(i),MI(j))+value;
           end
        end
    end
else
    sample_size=size(X,2);
    for index=1:sample_size
       for i=1:num_species
           for j=1:i
                 d_marginal{i,j}(X(i,index)+1,X(j,index)+1)...
                   =d_marginal{i,j}(X(i,index)+1,X(j,index)+1)+1/sample_size;   
           end
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