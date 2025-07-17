function [DCSnum] = calDCS()

% This function loads the IEMbase list of IMDs along with their phenotype
% and calculates the "Diagnostic Confidence Score of each IMD (DCS)"
% DCS is defined as:
% DCS = number of feasible characteristic metabolite biomarkers * quantified effect size
%
% Input:
%       This function does not take any input and loads the data from excel
%       sheet in data folder
%
% Output:
%       DCSnum: is a double array, representing DCS number for all of the
%       IMDs in the Excel sheet
%
% Author: Farid Zare, July 2025

% set data folder
path = mfilename('fullpath');
mainFolder = fileparts(fileparts(path));
dataFolder = fullfile(mainFolder, 'data');

% Load IEMbase data list
iembaseData = readtable(fullfile(dataFolder, 'IMD_mapped.xlsx'), 'ReadVariableNames', true, 'ReadRowNames',false);
highlightedIndex = logical(iembaseData.Highlighted); % Characteristic biomarkers are shown as highlighted in the IEMbase website
nlt = numel(highlightedIndex);

% Get all the age specific effect sizes
startID = find(ismember(fieldnames(iembaseData), 'Neonatal_birth_1mth_'));
endID = find(ismember(fieldnames(iembaseData), 'Adulthood__16yrs_'));
effectSize = table2array(iembaseData(:, startID:endID));

% Averagre over effect sizes among different ages
effectSizeAbs = abs(effectSize);
effectSizeMean = mean(effectSizeAbs, 2);

% Lower bound for recruiting patients for each IMD
IMDLb = zeros(nlt, 1);
IMDUb = zeros(nlt, 1);
patientNum = cell(nlt, 1);
charStatus = cell(nlt, 1);
% charScore = characterized score for each IMD

uniqueIMD = unique(iembaseData.IEMNosologyCode, 'stable');
nlt = numel(uniqueIMD);

DCSnum = zeros(nlt, 1);
for i = 1: nlt
    % Find IEMbase codes from IEMbase data
    IMDid = ismember(iembaseData.IEMNosologyCode, uniqueIMD(i));
    % Only consider feasible metabolites
    DCSnum(i) = sum(IMDid .* logical(iembaseData.feasible) .* effectSizeMean);

    % Get the effect size of each betabolite biomarker
    if DCSnum(i) > 6
        patientNum(IMDid) = {'5 - 10'};
        IMDLb(IMDid) = 5;
        IMDUb(IMDid) = 10;
        charStatus(IMDid) = {'well characterized'};
    elseif 3 <= DCSnum(i) && DCSnum(i) <= 6
        patientNum(IMDid) = {'10 - 15'};
        IMDLb(IMDid) = 10;
        IMDUb(IMDid) = 15;
        charStatus(IMDid) = {'moderately characterized'};
    elseif DCSnum(i) < 3
        patientNum(IMDid) = {'15 - 20'};
        IMDLb(IMDid) = 15;
        IMDUb(IMDid) = 20;
        charStatus(IMDid) = {'poorly characterized'};
    end
end

end