function [selected_fea, score] = feature_rank_val_pval(xa, y,k)
total_fea = size(xa,2);
no_of_data = size(xa,1);
score = [];
% for i=1: size(xa,2)
%     [r_val(i)] = relieff_test(xa(:,i), y,xa,k,'method','classification');
% end

[test,r_val,hit,miss]=relieff_our(xa,y,k,'method','classification');
% 
% corr = corrcoef(xa);
% corr(corr==1)=0;
% reduncancy = sum(corr,1)/(size(corr,2)-1);
% for i=1: size(xa,2)
%     if reduncancy(i)>0.5
%         r_val(i) = r_val(i)-reduncancy(i);
%     end
% end

%[r_val, test] = sort(r_val,'descend');

[val1, idx1] = max(r_val);
score  = [score val1];
selected_fea = [];
selected_fea = [selected_fea test(1)];
remain1 = setdiff(1:size(xa,2),selected_fea);
for jj = 2:size(xa,2)
    zz = [selected_fea test(jj)];
    [r_val1, hit, miss] = relieff_complementary(xa(:,zz),y,xa,k,'method','classification');
    % [r_val1(j)] = relieff_test(xa(:,zz),y,7,'method','classification');
    
    [val, idx] = max(r_val1);
    
     p = zeros(1,size(unique(y),1));
    for tt = 1:size(p,2)
        t1 = miss(:,find(y==tt));
        t2 = hit(:,find(y==tt));
        [h1, p(tt)] = ttest2(t1,t2,'Tail','right');
    end
    
    if val>val1 && isempty(find(p>0.1))
        %pp = remain1(idx);
        selected_fea = [selected_fea test(jj)];
        score  = [score val];
        val1 = val;
        %     else
        %         break;
    end
end