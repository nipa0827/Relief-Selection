clear all;
close all;
clc;
addpath('Libsvm/matlab');
addpath('D:\MS THESIS\code\Relief\data\balanced data');
%addpath('mi');
% for dataset = ['p' 'a' 'r' '$' 'n' '-' '!' 'l' 'z' '#' 'o' 'g' 'm' 'k' 'f' 'q' 'j' 'y' 'd' 'x' 's' ')' 'c' '>' 'b' '4' 'v' '(' 'i' '?' ',' '&' '8' '@' '<' 'e' '3' '1' '%' '0' '7' 'w' '^' '6' '*' 'h' '2' 'u' 't' '9' '5']
% for dataset = ['a' 'r' '-' '!' 'l' 'o' 'g' 'k' 'q' 'j' 's' ')' 'c' 'b' '4' '?' ',' '&' '8' '@' '<' 'e' '3' '6' '*' 'u']
% for dataset = ['<' 'e' '3']
%for dataset = ['a' 'r' '-'  'l' 'o' 'g' 'k' 'q' 'j' 's' ')' 'c' 'b' '4' '?' '3' '#' '_' '@' '+' '!' ]
 for dataset = ['q' 'j' 's' ')' 'c' 'b' '4' '?' '3' '#' '_' '@' '+' '!' ]
    
    %for dataset = ['-']
    cc = power(2,-5);
    number_neighbours=5;
    no_of_fold=10;
    ts_num_max = 5000;
    size_features = 0;
    LOO = 1;
    stability = 0;
    similarity_ratio = 0;
    ss=0;
    c_acc =0;
    
    if (isequal(dataset,'a'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('appendicitis.dat');
        file = 'appendicitis.txt';
        no_of_feature = 5;
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'b'))
        dim=14 ;
        nclass = 2;
        clabel = [1 2];
        data = dlmread('australian.dat');
        no_of_feature=7;
        file = 'australian.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'c'))
        nclass = 3;
        clabel = [1 2 3];
        data = dlmread('balance.dat');
        no_of_feature = 4;
        file = 'balance.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'e'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('banana.dat');
        file = 'banana.txt';
        no_of_feature = 2;
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'g'))
        nclass = 5;
        clabel = [1 2 3 4 5];
        data = dlmread('cleveland.dat');
        file = 'cleveland.txt';
        no_of_feature = 4;
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'j'))
        nclass = 6;
        clabel = [1 2 3 4 5 6];
        data = dlmread('dermatology_formatted.txt');
        no_of_feature=6;
        file = 'dermatology_formatted.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'k'))
        nclass = 8;
        clabel = [1 2 3 4 5 6 7 8];
        data = dlmread('ecoli.dat');
        no_of_feature = 6;
        file = 'ecoli.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'l'))
        nclass = 6;
        clabel = [1 2 3 4 5 6];
        data = dlmread('glass.txt');
        no_of_feature=2;
        file = 'glass.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'o'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('heart.txt');
        no_of_feature=4;
        file = 'heart.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'q'))
        nclass = 2;
        clabel = [1 2];
        no_of_feature=14;
        data = dlmread('ionosphere.data');
        file = 'ionosphere.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'r'))
        dim=4 ;
        nclass = 3;
        clabel = [1 2 3];
        data = dlmread('iris.data');
        no_of_feature=2;
        file = 'iris.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    %     if (isequal(dataset,';'))
    %         dim=4 ;
    %         nclass = 3;
    %         clabel = [1 2 3];
    %         data = dlmread('test.txt');
    %         no_of_feature=2;
    %         file = 'test.txt';
    %         [pathstr,name,ext] = fileparts(file);
    %     end
    if (isequal(dataset,'s'))
        nclass = 10;
        clabel = [1 2 3 4 5 6 7 8 9 10];
        data = dlmread('led7digit.dat');
        no_of_feature = 5;
        file = 'led7digit.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'u'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('magic.dat');
        file = 'magic.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'3'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('phoneme.dat');
        file = 'phoneme.txt';
        no_of_feature = 5;
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'4'))
        dim=8 ;
        nclass = 2;
        clabel = [1 2];
        data = dlmread('pima-indians-diabetes.data');
        no_of_feature=4;
        file = 'pima-indians-diabetes.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'6'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('ring.dat');
        file = 'ring.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'8'))
        nclass = 7;
        clabel = [1 2 3 4 5 6 7];
        data = dlmread('Segment.txt');
        no_of_feature=16;
        file = 'Segment.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'!'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('sonar data lebel first10fold.txt');
        file = 'sonar data lebel first10fold.txt';
        no_of_feature=8;
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'%'))
        nclass = 11;
        clabel = [1 2 3 5 6 7 8 9 11 12 13];
        data = dlmread('texture.dat');
        file = 'texture.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'&'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('titanic.dat');
        file = 'titanic.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'*'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('twonorm.dat');
        file = 'twonorm.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'#'))
        nclass = 10;
        clabel = [1 2 3 4 5 6 7 8 9 10];
        data = dlmread('Cardio.txt');
        no_of_feature=10;
        file = 'Cardio.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'_'))
        nclass = 3;
        clabel = [1 2 3];
        data = dlmread('waveform.txt');
        no_of_feature=9;
        file = 'waveform.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'@'))
        dim=22 ;
        nclass = 2;
        clabel = [1 2];
        data = dlmread('Parkinsons.txt');
        no_of_feature=12;
        file = 'Parkinsons.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'^'))
        nclass = 7;
        clabel = [1 2 3 4 5 6 7];
        data = dlmread('steel.txt');
        no_of_feature=26;
        file = 'steel.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    if (isequal(dataset,'+'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('breast.txt');
        no_of_feature=6;
        file = 'breast.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,')'))
        nclass = 2;
        clabel = [1 2];
        data = dlmread('wdbc.dat');
        no_of_feature  = 5;
        file = 'wdbc.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'-'))
        dim=13 ;
        nclass = 3;
        clabel = [1 2 3];
        data = dlmread('wine.data');
        no_of_feature=4;
        file = 'wine.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,','))
        nclass = 6;
        clabel = [1 2 3 4 5 6];
        data = dlmread('winequality-red.dat');
        file = 'winequality-red.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'<'))
        nclass = 7;
        clabel = [1 2 3 4 5 6 7];
        data = dlmread('winequality-white.dat');
        file = 'winequality-white.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    if (isequal(dataset,'?'))
        nclass = 10;
        clabel = [1 2 3 4 5 6 7 8 9 10];
        data = dlmread('yeast.data');
        no_of_feature=7;
        file = 'yeast.txt';
        [pathstr,name,ext] = fileparts(file);
    end
    
    
    fid = fopen('rval_result.txt', 'a');
    fprintf(fid,'\nDataset: %s\n', name);
    fclose(fid);
    y = data(:,1); %class id
    xa = data(:, 2:end); %feature
    
    %     Xmax = max(xa);
    %     Xmin = min(xa);
    %     Xdiff = Xmax-Xmin;
    %
    %     fea_mean = mean(xa);
    %     %xa = bsxfun(@rdivide,bsxfun(@minus,xa,mean(xa)),Xdiff);
    
    %[xa edges] = equal_width_quantization(xa, 5);
    % feature normalization
    %         f_min = min(xa);    f_max = max(xa);    f_tmp = f_max-f_min;
    %         r = 1./ (f_max - f_min);    r(f_tmp < 1e-10) = 1;
    %         xa = (xa - repmat(f_min,length(y),1)).*repmat(r,length(y),1);
    %
    Sapp = data(:, 1);
    [m, dim] = size(xa);
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
    % xa = my_norm(xa);
    r_val = zeros(size(xa,2),1);
    %Puloma--------------------------------------------------------------
    
    for kk=1:no_of_fold
        %         tic
        tr_idx = [];
        ts_idx = [];
        
        for jj = 1:nclass,
            
            idx_label = find(Sapp == clabel(jj));
            num = length(idx_label);
            rng(kk+jj);
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
        
        tr_fea = xa(tr_idx, :);
        tr_label = Sapp(tr_idx);
        ts_fea = xa(ts_idx, :);
        ts_label = Sapp(ts_idx);
        
        Xmax = max(tr_fea);
        Xmin = min(tr_fea);
        Xdiff = Xmax-Xmin;
        
        fea_mean = mean(tr_fea);
        tr_fea = bsxfun(@rdivide,bsxfun(@minus,tr_fea,fea_mean),Xdiff);
        
        %%%%%%%%%%%%%%%%%%%%% Equal Width %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 [tr_fea edges] = equal_width_quantization(tr_fea, 5);
        %                 for ii=1:size(tr_fea,2)
        %                     discretized(:,ii) = discretize(ts_fea(:,ii), edges(:,ii));
        %                 end
        %                 ts_fea = discretized;
        %
        
        %%%%%%%%%%%%%%%%%%%%%%% Equal Frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %              for j = 1:size(tr_fea,2)
        %                         [newtr_fea(:, j), edges] = nowquantizeMI_random(tr_fea(:,j), 5);
        %                         newts_fea(:, j) = testquantizeMI_random(ts_fea(:, j), edges);
        %              end
        %             tr_fea = newtr_fea;
        %             ts_fea = newts_fea;
        %
        
        matfile = strcat('mat\',name,'.mat');
        
        
        [sel_feature] = selection_ttest(tr_fea, tr_label);
        aa{kk}= sel_feature;
        %         tt{kk} = score;
        save(matfile,'aa');
        
        %
        %[sel_feature, score] = relieff(tr_fea,tr_label,7,'method','classification');
        % % %                 pp = load(matfile);
        % % %         no_feature = size(pp.aa{:,kk},2);
        %   sel_feature = sel_feature(1:4);
        %         [tr_fea edges] = equal_width_quantization(tr_fea, 5);
        %         for ii=1:size(tr_fea,2)
        %             discretized(:,ii) = discretize(ts_fea(:,ii), edges(:,ii));
        %         end
        %         ts_fea = discretized;
        %
        
        %         if isempty(sel_feature)
        %             continue;
        %         end
        
        
        %                     pp = load(matfile);
        %                     sel_feature = pp.aa{:,kk};
        %              sel_feature=sel_feature(:,1:9);
        %       sel_feature = sel_feature(1:3);
        %
        
        
        
        ts_fea = bsxfun(@rdivide,bsxfun(@minus,ts_fea,fea_mean),Xdiff);
        
        tr_fea1 = tr_fea(:,sel_feature);
        tr_label = Sapp(tr_idx);
        ts_fea1 = ts_fea(:,sel_feature);
        ts_label = Sapp(ts_idx);
        
        
        tr_fea = tr_fea1;
        ts_fea = ts_fea1;
        
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
        
        accuracy1 = sum(eq(pred_val,ts_label(:,1)))/size(ts_label,1);
        accuracy2 = sum(eq(pred_val2,ts_label(:,1)))/size(ts_label,1);
        accuracy3 = sum(eq(pred_val1,ts_label(:,1)))/size(ts_label,1);
        EVAL = EvalMetric(ts_label,pred_val);
        % EVAL.balance
        accuracies(kk)= accuracy1;
        accuracies2(kk)= accuracy2;
        accuracies3(kk)= accuracy3;
        selected(kk) = size(unique(sel_feature),2);
        
        fid = fopen('rval_result.txt', 'a');
        fprintf('Iter---> %d: accu(svm): %f\taccu(knn): %f\taccu(dt): %f\tsel: %d\tstability: %f\tc_avg_stabili: %f\n', kk, accuracy1,accuracy2,accuracy3,size(unique(sel_feature),2), stability,ss/(kk-1));
        fprintf(fid, 'Iter---> %d: accu: %f\taccu(knn): %f\taccu(dt): %f\tsel: %d\tstability: %f\tc_avg_stabili: %f\n', kk, accuracy1,accuracy2,accuracy3,size(unique(sel_feature),2), stability,ss/(kk-1));
        fclose(fid);
    end
    
    fid = fopen('rval_result.txt', 'a');
    fprintf('avg accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(accuracies),mean(accuracies2),mean(accuracies3),mean(selected),ss/9);
    fprintf(fid, 'avg accu(svm): %f\tavg accu(knn): %f\taccu(dt): %f\tsel: %f\tstability : %f\n', mean(accuracies),mean(accuracies2),mean(accuracies3),mean(selected),ss/9);
    fclose(fid);
    clear all;close all;
end


fclose('all');