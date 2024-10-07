function helperCommandWindowDisplay(MLtype,nAR,TrueNegative,FalseNegative,TruePositive,FalsePositive)
% Written by SE 165/265 (SP24) team.
% 
% helperCommandWindowDisplay displays a confusion chart, followed by:
% 
% True Positive Rate: the percent true positives relative to all 400 damage testing cases
% True Negative Rate: the percent true negatives relative to all 225 undamaged testing cases
% False Positive Rate: the percent of false negative relative to all 225 undamaged testing cases
% False Negative Rate: the percent of false negative relative to all 400 damaged testing cases 
% Overall Classification Accuracy: the percent of the sum of true positive plus true negatives relative to the sum of all 625 damaged and undamaged testing cases.

% 
%   INPUTS:
%       MLtype: boolean entry of 1 for 'Supervised Learning' or 0 for 'Unsupervised Learning'
%       nAR: number of auto-regressive coefficients (example: AR(5) has nAR = 5)
%       TrueNegative: number of true negative cases
%       FalseNegative: number of false negative cases
%       TruePositive: number of true positive cases
%       FalsePositive: number of FalsePositive negative cases
%
% Create Title
if MLtype == 1
    TitleString = '\t\t\t Supervised Learning';
else
    TitleString = '\t\t\t Unsupervised Learning';
end
disptitle = [TitleString,' Confusion Matrix AR(',num2str(nAR),') Model',repmat('\t',[1 10])];
fprintf([repmat('_',size(disptitle)),'\n\n',disptitle,'\n',repmat('_',size(disptitle)),'\n\n'])
% Print confusion chart
disp(array2table([TrueNegative,FalseNegative,TrueNegative+FalseNegative;...
                  FalsePositive,TruePositive,FalsePositive+TruePositive;...
                  TrueNegative+FalsePositive,FalseNegative+TruePositive,NaN],...
    'RowNames',{'Predicted undamaged';'Predicted damaged';'Total'},...
    'VariableNames',{'Actual Undamaged','Actual Damaged','Total'}))
% Print performance metrics
fprintf('\nTrue Positive (damage) Rate: %.1f%%\n' , 100*TruePositive/...
    (FalseNegative + TruePositive))
fprintf ('True Negative (undamaged) Rate: %.1f%%\n', 100*TrueNegative/...
    (TrueNegative + FalsePositive))
fprintf('False Positive Rate: %.1f%%\n', 100*FalsePositive/...
    (TrueNegative + FalsePositive))
fprintf('False Negative Rate: %.1f%%\n', 100*FalseNegative/...
    (FalseNegative + TruePositive))
fprintf('Overall Classification Accuracy: %.1f%%\n\n\n', 100*(TruePositive...
    + TrueNegative)/(FalseNegative + TruePositive + TrueNegative + FalsePositive))
%end function