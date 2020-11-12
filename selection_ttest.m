function [selected_fea,score] = rank(xa, y)
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
    
    try
    [h1,p1]=ttest2(miss(:,test(i)),hit(:,test(i)),'Tail','right');
    catch
        h1  = 0;
    end
    
    %if weight(i) >(chi2inv(0.99,df)/(2*size(xa,1)))
    if p1 < 0.01
        selected_fea = [selected_fea test(i)];
        score = [score weight(i)];
    end
end
