function [selected_fea,score] = selection_ttest(xa, y)
total_fea = size(xa,2);
no_of_data = size(xa,1);
idx_val = 0;

% for i=1: size(xa,2)
%    [r_val(i)] = relieff_test(xa(:,i), y,xa,7,'method','classification');
% end


[test,weight,hit,miss] = relieff_our(xa, y,7,'method','classification');

%[xa edges] = equal_width_quantization(xa, 5);
selected_fea = [];
score = [];
df = 1;


for i=1:size(xa,2)
    %df = size(unique(xa(:,i)),1)-1;
    p = zeros(1,size(unique(y),1));
    for jj = 1:size(p,2)
        t1 = miss(find(y==jj),:);
        t2 = hit(find(y==jj),:);
        [h1, p(jj)] = ttest2(t1(:,test(i)),t2(:,test(i)),'Tail','right');
    end
    %  try
    % [h1,p1]=ttest2(miss(:,test(i)),hit(:,test(i)),'Tail','right');
    %     catch
    %         h1  = 0;
    %     end
    %  p = [.002 ,.001];
    %
    %     if p <.05
    %         x=1;
    %     end
    
    %if weight(i) >(chi2inv(0.99,df)/(2*size(xa,1)))
    if p < 0.1
        selected_fea = [selected_fea test(i)];
        score = [score weight(i)];
    end
end
