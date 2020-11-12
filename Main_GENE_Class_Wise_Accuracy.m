clear all;
close all;
clc;
addpath('Libsvm/matlab');
addpath('D:\MS THESIS\code\Relief\data\gene data');
for dataset = ['5' ]
    cc = power(2,-5);
    number_neighbours=5;
    LOO=0;
    ss=0;
    similarity_ratio = 0;
    size_features = 0;
    pp = 0;
    stability=0;
    how_many =5;
    
    %paramteres
    max_qua_level = 50;
    no_of_fold=5;
    ts_num_max = 5000;
    
    if (isequal(dataset,'a'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('appendicitis.data');
        file = 'appendicitis.data';
        no_of_feature = 3;
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'!'))
        nclass = 2; %LOO
        clabel = [1 2];
        data = dlmread('merged_GSE106291.csv');
        no_of_feature=6;
        file = 'merged_GSE106291.csv';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'2'))
        nclass = 2; %LOO
        clabel = [1 2];
        data = dlmread('colon.txt');
        no_of_feature=5;
        file = 'colon.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'3'))
        nclass = 9; %LOO
        clabel = [1 2 3 4 5 6 7 8 9];
        data = dlmread('lymphoma.txt');
        no_of_feature=3;
        file = 'lymphoma.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'4'))
        nclass = 5;
        clabel = [1 2 3 4 5];
        data = dlmread('Lung.txt');
        no_of_feature=7;
        file = 'Lung.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'5'))
        dim=4 ;
        nclass = 2;
        clabel = [1 2 ];
        data = dlmread('CNS.txt');
        no_of_feature=5;
        file = 'CNS.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'6'))
        dim=8 ;
        nclass = 2;
        clabel = [1 2];
        data = dlmread('LEUKEMIA.txt');
        no_of_feature=2;
        file = 'LEUKEMIA.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'7'))
        nclass = 3;
        clabel = [1 2 3];
        data = dlmread('LEUKEMIA3C.txt');
        no_of_feature=4;
        file = 'LEUKEMIA3C.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'8'))
        dim=13 ;
        nclass = 4;
        clabel = [1 2 3 4];
        data = dlmread('LEUKEMIA4C.txt');
        no_of_feature=4;
        file = 'LEUKEMIA4C.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'9'))
        nclass = 3;
        clabel = [1 2 3];
        data = dlmread('MLL.txt');
        no_of_feature=4;
        file = 'MLL.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'@'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('OVARIAN.txt');
        no_of_feature=3;
        file = 'OVARIAN.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'1'))
        nclass = 3;
        clabel = [1 2 3];
        data = dlmread('SRBCT.txt');
        no_of_feature=5;
        file = 'SRBCT.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    
    fid = fopen('gene_result.txt', 'a');
    fprintf(fid,'\nDataset: %s\n', name);
    fprintf(fid,'\nNumber of selected feature: %d\n', no_of_feature);
    fclose(fid);
    
    
    xa = data(:, 2:end);
    
    Xmax = max(xa);
    Xmin = min(xa);
    Xdiff = Xmax-Xmin;
    
    xa = bsxfun(@rdivide,bsxfun(@minus,xa,mean(xa)),Xdiff);
    
    count=1;
    for ii=1:size(xa,2)
        if length(unique(xa(:,ii)))==1
            temp1(count)=ii;
            count=count+1;
        end
    end
    if count>1
        xa(:,temp1)=[];
    end
    clearvars temp1 count ;
    
    
    
    Sapp = data(:, 1);
    [m, dim] = size(xa);
    
    opts= struct;
    opts.att_1split= 2;
    opts.quaztization_level= 5;
    opts.dim= dim;
    
    
    if LOO==1
        acc=0;
        
        
        % %Leave-one-out cross-validation
        
        for i = 1 : size(xa,1)
            ind=1:size(xa,1);
            ind=ind(ind~=i);
            tr_fea = xa(ind,:);
            tr_label = Sapp(ind,:);
            ts_fea = xa(i,:);
            ts_label = Sapp(i,:);
            
            [tr_fea edges] = equal_width_quantization(xa(ind,:), 3);
            
            for ii=1:size(tr_fea,2)
                ts_fea(:,ii) = discretize(xa(i,ii), edges(:,ii));
            end
            %%%%%%%%%%%%%%%%%%%%%%% Equal Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  for j = 1:size(tr_fea,2)
            %             [newtr_fea(:, j), edges] = nowquantizeMI_random(tr_fea(:,j), 5);
            %             newts_fea(:, j) = testquantizeMI_random(ts_fea(:, j), edges);
            %  end
            % tr_fea = newtr_fea;
            % ts_fea = newts_fea;
            %%%%%%%%%%%%%%%%%%%%%%% Equal Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         [steps,sel_flag,rel,red,cond_red] = select_features(tr_fea,Sapp(ind,:),'direction','f','method','jmi','degree', 2);
            %         selectedFeatures=steps(1,:);
            %         aa{i}= selectedFeatures;
            %         save('mat\Lung_Cancer.mat','aa');
            pp = load('mat\wine');
            selectedFeatures = pp.aa{:,i};
            selectedFeatures=selectedFeatures(:,1:9);
            % aa{kk}= selectedFeatures;
            % if size(selectedFeatures,2)>50
            %         selectedFeatures=selectedFeatures(:,1:how_many);
            % end
            
            
            %         if i == 1
            %             temp = selectedFeatures;
            %             temp_fea = tr_fea;
            %         end
            %         if i >1
            %             %     C = intersect(temp,selectedFeatures);
            %             %     similarity_ratio = similarity_ratio + length(C)/length(temp);
            %             stability=infostability(unique(temp),unique(selectedFeatures),temp_fea,tr_fea);
            %             ss=ss+stability;
            %             temp = selectedFeatures;
            %             temp_fea = tr_fea;
            %         end
            tr_fea = tr_fea(:, selectedFeatures);
            ts_fea = ts_fea(:, selectedFeatures);
            % tree = fitctree(xa(ind,selectedFeatures),Sapp(ind));
            % label = predict(tree,xa(i,selectedFeatures) );
            % acc =( label==Sapp(i) ) + acc;
            disp(selectedFeatures);
            % % %KNN
            model2 = fitcknn(tr_fea,tr_label);
            pred_val2 = predict(model2,ts_fea);
            
            % % %     decision tree
            
            model1 = fitctree(tr_fea, tr_label);
            pred_val1 = predict(model1, ts_fea);
            
            
            
            options = [ '-t 0 '  ' -c 1' ];      % Libsvm parameter setting (linear SVM is used)
            
            model = svmtrain(double(Sapp(ind)), sparse(xa(ind,selectedFeatures)), options);
            clear tr_fea;
            
            [C, Acc, d2p] = svmpredict(double(Sapp(i)), sparse(xa(i,selectedFeatures)), model);                % svm test
            clear ts_fea;
            pred_val=C;
            
            acc = (pred_val==Sapp(i));
            acc1 = (pred_val1==Sapp(i));
            acc2 = (pred_val2==Sapp(i));
            
            acc_f(i)= acc;
            acc_f1(i)= acc1;
            acc_f2(i)= acc2;
            
            
            selected = size(unique(selectedFeatures),2);
            select(i) = selected ;
            
            fid = fopen('gene_result.txt.txt', 'a');
            
            fprintf('Iter---> %d:\t accu(svm): %f\t accu(knn): %f\t accu(dt): %f\t sel: %f\t stability : %f\n', i, acc_f(i),acc_f2(i),acc_f1(i),select(i),ss/9);
            fprintf(fid, 'Iter---> %d:\t accu(svm): %f\t accu(knn): %f\t accu(dt): %f\t sel: %f\t stability : %f\n', i, acc_f(i),acc_f2(i),acc_f1(i),select(i),ss/9);
            fclose(fid);
            
        end
        %     S = kuncheva_stability(featidx,dim);
        
        fid = fopen('gene_result.txt', 'a');
        %     fprintf('avg accu: %f\tsel: %f\tstability : %f\n', mean(accuracies),mean(selected),ss/9);
        fprintf('avg accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(acc_f),mean(acc_f2),mean(acc_f1),mean(select),ss/9);
        fprintf(fid,'avg accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(acc_f),mean(acc_f2),mean(acc_f1),mean(select),ss/9);
        fclose(fid);
    else
        
        for zz=1:10
            
            for kk=1:no_of_fold
                
                tr_idx = [];
                ts_idx = [];
                
                for jj = 1:nclass,
                    
                    idx_label = find(Sapp == clabel(jj));
                    num = length(idx_label);
                    rng(kk+jj+zz);
                    idx_rand = randperm(num);
                    tr_num=size(idx_label,1)-ceil(size(idx_label,1)/no_of_fold);
                    if num > tr_num + ts_num_max
                        tr_idx = [tr_idx; idx_label(idx_rand(1:tr_num))];
                        ts_idx = [ts_idx; idx_label(idx_rand(tr_num+1:tr_num+ts_num_max))];
                    else
                        tr_idx = [tr_idx; idx_label(idx_rand(1:tr_num))];
                        ts_idx = [ts_idx; idx_label(idx_rand(tr_num+1:end))];
                    end
                end
                tr_fea = zeros(length(tr_idx), dim);
                tr_label = zeros(length(tr_idx), 1);
                ts_fea = zeros(length(ts_idx), dim);
                ts_label = zeros(length(ts_idx), 1);
                
                
                tr_fea = xa(tr_idx,:);
                tr_label = Sapp(tr_idx);
                ts_fea = xa(ts_idx,:);
                ts_label = Sapp(ts_idx);
                % tic
                %%%%%%%%%%%%%%%%%%%%%%% Equal Width %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             %[tr_fea, edges] = equal_width_quantization(tr_fea, 3);
                %             for ii=1:size(tr_fea,2)
                %                discretized(:,ii) = discretize(ts_fea(:,ii), edges(:,ii));
                %             end
                %             ts_fea = discretized;
                %
                %             fprintf('%s\n',dataset);
                %%%%%%%%%%%%%%%%%%%%%%% Equal Width %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%% Equal Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %         for j = 1:size(tr_fea,2)
                %             [newtr_fea(:, j), edges] = nowquantizeMI_random(tr_fea(:,j), 5);
                %             newts_fea(:, j) = testquantizeMI_random(ts_fea(:, j), edges);
                %         end
                %         tr_fea = newtr_fea;
                %         ts_fea = newts_fea;
                %%%%%%%%%%%%%%%%%%%%%%% Equal Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % [selectedFeatures,weight]= relieff(tr_fea,tr_label,7,'method','classification');
                % no_of_feature = feature(kk);
                % [selectedFeatures,score] = relieff(tr_fea, tr_label,7,'method','classification');
                %   [selectedFeatures,score] = feature_rank(tr_fea, tr_label);
                
                matfile = strcat('mat\',name,'.mat');
                
                %         selectedFeatures=steps(1,:);
                %                                     aa{zz,kk}= selectedFeatures;
                %                                      tt{zz,kk} = score;
                %                                  save(matfile,'aa','tt');
                pp = load(matfile);
                selectedFeatures = pp.aa{:,kk};
                %                 %             selectedFeatures=selectedFeatures(:,1:9);
                selectedFeatures = selectedFeatures(1:no_of_feature);
                
                
                
                tr_fea = tr_fea(:, selectedFeatures);
                ts_fea = ts_fea(:, selectedFeatures);
                
                % % %KNN
                model2 = fitcknn(tr_fea,tr_label);
                pred_val2 = predict(model2,ts_fea);
                % % %     decision tree
                
                % model = fitctree(tr_fea, tr_label);
                %     pred_val = predict(model, ts_fea);
                % % %     decision tree
                model1 = fitctree(tr_fea, tr_label);
                pred_val1 = predict(model1, ts_fea);
                
                % % %     SVM
                
                
                
                % feature normalization
                %     f_min = min(tr_fea);    f_max = max(tr_fea);    f_tmp = f_max-f_min;
                %     r = 1./ (f_max - f_min);    r(f_tmp < 1e-10) = 1;
                %     tr_fea = (tr_fea - repmat(f_min,length(tr_label),1)).*repmat(r,length(tr_label),1);
                %
                c_chosen(1) = 1;
                %     options = [ '-s 0 -t 0 ' '-g ' num2str(power(2, -7)) ' -c ' num2str(cc(c_chosen(1)))];      % Libsvm parameter setting (linear SVM is used)
                options = [ '-t 0 '  ' -c 1' ];      % Libsvm parameter setting (linear SVM is used)
                
                model = svmtrain(double(tr_label), sparse(tr_fea), options);
                clear tr_fea;
                %
                %      ts_fea = (ts_fea - repmat(f_min,length(ts_label),1)).*repmat(r,length(ts_label),1);     % normalize the test feature
                [C, Acc, d2p] = svmpredict(double(ts_label), sparse(ts_fea), model);                % svm test
                clear ts_fea;
                pred_val=C;
                
                % % %     SVM
                
                
                
                cls = unique(ts_label);
                
                c1 = find(ts_label==1);
                class_1_svm(zz,kk) =  sum(eq(pred_val(c1,:),ts_label(c1,:)))/size(c1,1);
                class_1_knn(zz,kk) = sum(eq(pred_val2(c1,:),ts_label(c1,:)))/size(c1,1);
                class_1_DT(zz,kk) = sum(eq(pred_val1(c1,:),ts_label(c1,:)))/size(c1,1);
                
                c2 = find(ts_label==2);
                class_2_svm(zz,kk) =  sum(eq(pred_val(c2,:),ts_label(c2,:)))/size(c2,1);
                class_2_knn(zz,kk) = sum(eq(pred_val2(c2,:),ts_label(c2,:)))/size(c2,1);
                class_2_DT(zz,kk) = sum(eq(pred_val1(c2,:),ts_label(c2,:)))/size(c2,1);
                
                c3 = find(ts_label==3);
                class_3_svm(zz,kk) =  sum(eq(pred_val(c3,:),ts_label(c3,:)))/size(c3,1);
                class_3_knn(zz,kk) = sum(eq(pred_val2(c3,:),ts_label(c3,:)))/size(c3,1);
                class_3_DT(zz,kk) = sum(eq(pred_val1(c3,:),ts_label(c3,:)))/size(c3,1);
                
                c4 = find(ts_label==4);
                class_4_svm(zz,kk) =  sum(eq(pred_val(c4,:),ts_label(c4,:)))/size(c4,1);
                class_4_knn(zz,kk) = sum(eq(pred_val2(c4,:),ts_label(c4,:)))/size(c4,1);
                class_4_DT(zz,kk) = sum(eq(pred_val1(c4,:),ts_label(c4,:)))/size(c4,1);
                
                c5 = find(ts_label==5);
                class_5_svm(zz,kk) =  sum(eq(pred_val(c5,:),ts_label(c5,:)))/size(c5,1);
                class_5_knn(zz,kk) = sum(eq(pred_val2(c5,:),ts_label(c5,:)))/size(c5,1);
                class_5_DT(zz,kk) = sum(eq(pred_val1(c5,:),ts_label(c5,:)))/size(c5,1);
                
                c6 = find(ts_label==6);
                class_6_svm(zz,kk) =  sum(eq(pred_val(c6,:),ts_label(c6,:)))/size(c6,1);
                class_6_knn(zz,kk) = sum(eq(pred_val2(c6,:),ts_label(c6,:)))/size(c6,1);
                class_6_DT(zz,kk) = sum(eq(pred_val1(c6,:),ts_label(c6,:)))/size(c6,1);
                
                c7 = find(ts_label==7);
                class_7_svm(zz,kk) =  sum(eq(pred_val(c7,:),ts_label(c7,:)))/size(c7,1);
                class_7_knn(zz,kk) = sum(eq(pred_val2(c7,:),ts_label(c7,:)))/size(c7,1);
                class_7_DT(zz,kk) = sum(eq(pred_val1(c7,:),ts_label(c7,:)))/size(c7,1);
                
                c8 = find(ts_label==8);
                class_8_svm(zz,kk) =  sum(eq(pred_val(c8,:),ts_label(c8,:)))/size(c8,1);
                class_8_knn(zz,kk) = sum(eq(pred_val2(c8,:),ts_label(c8,:)))/size(c8,1);
                class_8_DT(zz,kk) = sum(eq(pred_val1(c8,:),ts_label(c8,:)))/size(c8,1);
                
                c9 = find(ts_label==9);
                class_9_svm(zz,kk) =  sum(eq(pred_val(c9,:),ts_label(c9,:)))/size(c9,1);
                class_9_knn(zz,kk) = sum(eq(pred_val2(c9,:),ts_label(c9,:)))/size(c9,1);
                class_9_DT(zz,kk) = sum(eq(pred_val1(c9,:),ts_label(c9,:)))/size(c9,1);
                
                accuracy1 = sum(eq(pred_val,ts_label(:,1)))/size(ts_label,1);
                accuracy2 = sum(eq(pred_val2,ts_label(:,1)))/size(ts_label,1);
                accuracy3 = sum(eq(pred_val1,ts_label(:,1)))/size(ts_label,1);
                EVAL = EvalMetric(ts_label,pred_val);
                % EVAL.balance
                accuracies(zz,kk)= accuracy1;
                accuracies2(zz,kk)= accuracy2;
                accuracies3(zz,kk)= accuracy3;
                selected(zz,kk) = size(unique(selectedFeatures),2);
                
                fid = fopen('gene_result.txt', 'a');
                fprintf('Iter---> %d: number----> %d:  accu(svm): %f\taccu(knn): %f\taccu(dt): %f\tsel: %d\tstability: %f\tc_avg_stabili: %f\n', zz, kk, accuracy1,accuracy2,accuracy3,size(unique(selectedFeatures),2), stability,ss/(kk-1));
                fprintf(fid, 'Iter---> %d: number----> %d:  accu: %f\taccu(knn): %f\taccu(dt): %f\tsel: %d\tstability: %f\tc_avg_stabili: %f\n', zz, kk, accuracy1,accuracy2,accuracy3,size(unique(selectedFeatures),2), stability,ss/(kk-1));
                fclose(fid);
            end
            % save('EF/EF_selected_iris_jmi.mat','aa');
            %  minsize = min(cellfun('size', featidx, 2));
            % newc = cellfun(@(x) x(1:minsize), featidx, 'uniformoutput',false);
            % featidx1=cell2mat(newc);
            
            % S = kuncheva_stability(featidx1,dim);
            % [assignment,cost] = munkres(featidx1)
            % Avg = getAvgEval(evals);
            
        end
        fid = fopen('gene_result.txt', 'a');
        fprintf('avg accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(accuracies(:)),mean(accuracies2(:)),mean(accuracies3(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'avg accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(accuracies(:)),mean(accuracies2(:)),mean(accuracies3(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-1 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_1_svm(:)),mean(class_1_knn(:)),mean(class_1_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-2 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_2_svm(:)),mean(class_2_knn(:)),mean(class_2_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-3 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_3_svm(:)),mean(class_3_knn(:)),mean(class_3_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-4 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_4_svm(:)),mean(class_4_knn(:)),mean(class_4_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-5 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_5_svm(:)),mean(class_5_knn(:)),mean(class_5_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-6 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_6_svm(:)),mean(class_6_knn(:)),mean(class_6_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-7 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_7_svm(:)),mean(class_7_knn(:)),mean(class_7_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-8 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_8_svm(:)),mean(class_8_knn(:)),mean(class_8_DT(:)),mean(selected(:)),ss/9);
        fprintf(fid, 'Class-9 accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(class_9_svm(:)),mean(class_9_knn(:)),mean(class_9_DT(:)),mean(selected(:)),ss/9);
        fclose(fid);
        
    end
    
    clear all;close all;
end
fclose('all');