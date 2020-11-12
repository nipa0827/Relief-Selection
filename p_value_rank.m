function [selected_fea,score] = p_value_rank(xa, y)
total_fea = size(xa,2);
no_of_data = size(xa,1);
idx_val = 0;

% for i=1: size(xa,2)
%    [r_val(i)] = relieff_test(xa(:,i), y,xa,7,'method','classification');
% end


[test,weight,hit,miss] = relieff_ttest(xa, y,7,'method','classification');

%[xa edges] = equal_width_quantization(xa, 5);
selected_fea = [];
score = [];
df = 1;

for i=1:size(xa,2)
%df = size(unique(xa(:,i)),1)-1;

[h1,p1(i)]=ttest2(miss(:,i),hit(:,i),'Tail','right');


end

[~,sorted] = sort(p1,'ascend');

for i=1:size(xa,2)
%df = size(unique(xa(:,i)),1)-1;


%if weight(i) >(chi2inv(0.99,df)/(2*size(xa,1)))
if p1(sorted(i)) < 0.01
selected_fea = [selected_fea sorted(i)];
score = [score weight(sorted(i))];
end
end