angle=1;
figure(1)
ints_level_error={};
% for i=1:length(coc_mat)
for i=1:2
    coc_level=coc_mat{i,1};
    R_level=R_mat{i,1};
    ints_level=ints_mat{i,1};
    k=find(ints_level<=0);
    ints_level_error{i,1}=k;
    coc=coc_level(angle,:);
    R=R_level(angle,:);
    ints=ints_level(angle,:);
    for j=1:length(coc)
        if j~=1
            x=sum(ints(1:j))+sum(R(1:j-1))+R(j)/2;
        else
            x=sum(ints(1:j))+R(j)/2;
        end
        y=coc(j);
        r=R(j)/2;
        sita=-0.05:0.05:2*pi;
        subplot(1,2,1)
        plot(x+r*cos(sita),y+r*sin(sita),'k');
        xlim([4.9e5 5.01e5]);
        hold on;
        subplot(1,2,2)
        plot(x+r*cos(sita),y+r*sin(sita),'k');
        xlim([-0.01e5 0.1e5]);
        hold on;
    end
end
xlim([0 1000]);
ylim([0 1000]);