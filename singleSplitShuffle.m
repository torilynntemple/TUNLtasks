%%SplitHalf reliability with random sorting of trial order 100 times.
%
%INPUT:
%   -trial: n X m X p, calcium traces. Where n is the cell number, m is
%   frames, p the number of trials.
%   -cl: Vector containing different trial lengths
%   -cut:Index of what kind of anchoring, 1 for Delay period, 2 for Front
%   Anchored, 3 for Front Anchored. Any excess transient is ignored.
%OUTPUT:
%   -out: Correlation matrix for each split half itteration.

function [out] = singleSplitShuffle(trial,cl,cut)
if ~isempty(trial)
    if ~(length(trial(1,1,:)) == 0)
    out = zeros(length(trial(:,1,1)),1);
    splitcorr = zeros(length(trial(:,1,1)),100);
    if ~isempty(trial) && length(trial(1,1,:))>1
        temp = permute(trial,[3 2 1]);
        if cut == 1 %Delay
            %do nothing
        elseif cut == 2 %Front Anchored: eliminate back end excess
            if ~(min(temp) == length(temp))
                temp = temp(:,1:min(cl),:);
            end
        elseif cut == 3 %Back Anchored: eliminate front end excess
            if ~(min(temp) == length(temp))
                temp = temp(:,end-min(cl)+1:end,:);
            end
        end
        count = 1;
        for i = 1 : length(temp(1,1,:))
            for j = 1: 100
                rearange = randperm(length(trial(1,1,:)));
                h1 = nanmean(temp(rearange(1:round(length(rearange)/2)),:,i),1);
                h2 = nanmean(temp(rearange(round(length(rearange)/2)+1:end),:,i),1);
                %             if ~isempty(find(isnan(h1),1))
                %                 h1(find(isnan(h1))) = 0;
                %             elseif ~isempty(find(isnan(h2),1))
                %                 h2(find(isnan(h2))) = 0;
                %             end
                splitcorr(i,j) = corr2(h1,h2);
            end
        end
        out = mean(splitcorr,2);
    end
    else 
        out = NaN;
    end
else
    out = NaN;
end
end
