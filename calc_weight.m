function [r_val,c] = calc_weight(f, y)
cnt1 =1;
%rc1 = 0;
rf = 0;
KNN = 0;
NO_OF_CLASS = length(unique(y));
distance = -100;
class_lbl = unique(y);
no_of_data = size(f,1); % check size 2/1
c=[];
for i=1:NO_OF_CLASS
    rc1 = 0;
    cnt1 =1;
    arr1=[];
    x = [];
    z = [];
    miss = [];
    miss_count = 1;
    for k=1:no_of_data
        if y(k) == class_lbl(i)
            arr1(cnt1,:) = f(k,:);
            z(cnt1) = y(k);
            cnt1 = cnt1+1;
        end
        if  ~(y(k) == class_lbl(i))
            x(miss_count) = y(k);
            miss(miss_count,:) = f(k,:);
            miss_count = miss_count +1;          
        end
        
    end
    [rc1,K,zz] = subCalcWeight(arr1,f,miss,z,y);
    rf = rf+rc1;
    rc1 = rc1/size(arr1,1);
    KNN = KNN +size(arr1,1)*K;
    c = [c ;zz];
end
%r_val = rf/(no_of_data*7);
r_val = rf;