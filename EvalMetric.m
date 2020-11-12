function EVAL = EvalMetric(ACTUAL,PREDICTED)
% This fucntion evaluates the performance of a classification model by 
% calculating the common performance measures: Accuracy, Sensitivity, 
% Specificity, Precision, Recall, F-Measure, G-mean.
% Input: ACTUAL = Column matrix with actual class labels of the training
%                 examples
%        PREDICTED = Column matrix with predicted class labels by the
%                    classification model
% Output: EVAL = Row matrix with all the performance measures
EVAL = struct;

idx = (ACTUAL()==1);

p = length(ACTUAL(idx));
n = length(ACTUAL(~idx));
N = p+n;

tp = sum(ACTUAL(idx)==PREDICTED(idx));
tn = sum(ACTUAL(~idx)==PREDICTED(~idx));
fp = n-tn;
fn = p-tp;

tp_rate = tp/p;
tn_rate = tn/n;

accuracy = (tp+tn)/N;
sensitivity = tp_rate;
specificity = tn_rate;
precision = tp/(tp+fp);
recall = sensitivity;
f_measure = 2*((precision*recall)/(precision + recall));
gmean = sqrt(tp_rate*tn_rate);

pd=tp/(tp+fn);
pf=fp/(fp+tn);
balance =(1-sqrt((1- pd^2)+pf^2)/2^0.5)* 100;

EVAL.accuracy = accuracy;
EVAL.sensitivity = sensitivity;
EVAL.specificity = specificity;
EVAL.precision = precision;
EVAL.recall = recall;
EVAL.f_measure = f_measure;
EVAL.gmean = gmean;
EVAL.balance = balance;
% EVAL = [accuracy sensitivity specificity precision recall f_measure gmean];