%Shuffle Touchscreen data 100 times and compares to split half correlation
%
%INPUT:
%   -ms: miniscope structure, pertinant variables include:
%       *FiltTraces: n X m matrix, where n is the frame number and m the
%       cell number.
%   -events: string containing touchscreen event information
%OUTPUT:
%   - to be determined....
%

function [out] = TouchScreenShuffle(ms,events)
%seperate frame map from events file
for i = 1 : length(events)
    frameMap(i,1) = str2num(char(events(i,1)));
end
tic
out.ScorrDelaycorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrDelayincorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrDelayccor = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrDelayicor = zeros(length(ms.FiltTraces(1,:)),100);

out.DelaySplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.DelaySplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.DelaySplithalfccor = zeros(length(ms.FiltTraces(1,:)),1);
out.DelaySplithalficor = zeros(length(ms.FiltTraces(1,:)),1);

out.ScorrFrontcorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrFrontincorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrFrontccor = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrFronticor = zeros(length(ms.FiltTraces(1,:)),100);

out.FrontSplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.FrontSplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.FrontSplithalfccor = zeros(length(ms.FiltTraces(1,:)),1);
out.FrontSplithalficor = zeros(length(ms.FiltTraces(1,:)),1);

out.ScorrBackcorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrBackincorrect = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrBackccor = zeros(length(ms.FiltTraces(1,:)),100);
out.ScorrBackicor = zeros(length(ms.FiltTraces(1,:)),100);

out.BackSplithalfcorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.BackSplithalfincorrect = zeros(length(ms.FiltTraces(1,:)),1);
out.BackSplithalfccor = zeros(length(ms.FiltTraces(1,:)),1);
out.BackSplithalficor = zeros(length(ms.FiltTraces(1,:)),1);

for itteration = 1 : 101
    if itteration> 1
        msC.FiltTraces = CShuffle(ms.FiltTraces);
    end
    ctInd = [];                                                                 %Correct Trial Indices
    itInd = [];                                                                 %Incorrect Trial Indices
    ctCount = [];                                                               %Correct Trial count
    itCount = [];                                                               %Incorrect Trial count
    crcRow = [];                                                                %Correct correction trial row index
    criRow = [];                                                                %Incorrect correction trial row index
    crcInd = [];                                                                %Correct correction trial indices
    criInd = [];                                                                %Incorrect correction trial indices
    cri = 1;                                                                    %Correct correction trial index marker
    crc = 1;                                                                    %Incorrect correction trial index marker
    dc = [];
    di = [];
    dcc = [];
    dic = [];
    dccount = 1;
    dicount = 1;
    dcccount = 1;
    diccount = 1;
    Dcoef = zeros(length(ms.FiltTraces(1,:)),1);
    Acoef = Dcoef;
    
    %find all relevant events
    Trial = find(contains(events(:,5),'1'));
    Trial = [1 ; Trial];
    correct = find(contains(events(:,3),'Correct'));
    incorrect = find(contains(events(:,3),'Incorrect'));
    correction = find(contains(events(:,3),'correction'));
    delay = find(contains(events(:,6), '1'));
    delay2 = find(contains(events(:,6), '2'));
    if length(delay) > length(delay2)
        delay(end) = [];
    end
    timedelay = delay2-delay;
    dFrames = mode(timedelay);
    
    %find nearest events to correct/incorrect choices
    Cneighbr = knnsearch(Trial,correct,'k',2);
    Ineighbr = knnsearch(Trial,incorrect,'k',2);
    Cneighbr = sort(Cneighbr,2);
    Ineighbr = sort(Ineighbr,2);
    
    flagInf = 0;
    
    %Sort frames between correct/incorrect/correction trials
    %Correct trial sorting
    for i = 1 : length(Cneighbr)
        %Sorting lower bounds
        while min(Trial(Cneighbr(i,:)))> correct(i,1)
            temp = Cneighbr(find(Cneighbr == min(Cneighbr(i,:)),1)); %i == ic mod
            Cneighbr(find(Cneighbr == min(Cneighbr(i,1)),1)) = min(Cneighbr(i,:))-1;
            Cneighbr(find(Cneighbr(:,2) == max(Cneighbr(i,2)),1),2) = temp(1);
        end
        %Sorting upper bounds
        while ~any(isnan(Cneighbr(i,:)))&& flagInf <3
            if max(Trial(Cneighbr(i,:))) < correct(i,1)  %i == ic mod
                if(max(Cneighbr(i,:)))<= length(Trial)
                    Cneighbr(i,find(max(Cneighbr(i,2)) == Cneighbr(i,:),1)) = NaN;
                    Cneighbr(i,find(min(Cneighbr(i,1)) == Cneighbr(i,:),1)) = Cneighbr(i,find(min(Cneighbr(i,:)) == Cneighbr(i,1),1))+1;
                else
                    Cneighbr(find(min(Cneighbr == Cneighbr(i,:)),1)) = Cneighbr(find(min(Cneighbr == Cneighbr(i,:)),1)+1);
                end
            end
            flagInf = flagInf +1;
        end
        ctInd(i,1) = Trial(nanmin(Cneighbr(i,:)));
        
        %Identifying Correct correction trials
        if any(ctInd(i,1) == correction(:,1))
            crcRow(i,1) = i;
            crc = crc + 1;
        end
        %setting trial length
        if any(isnan(Cneighbr(i,:)))
            ctCount(i,1) = length(frameMap) - Trial(nanmax(Cneighbr(i,:)));
        else
            ctCount(i,1) = abs(diff(Trial(Cneighbr(i,:))));
        end
        flagInf = 0;
    end
    
    for i = 1: length(Ineighbr)
        while min(Trial(Ineighbr(i,:)))> incorrect(i,1)% && ~(min(Ineighbr(i,:) == Ineighbr(i-1,:)))
            temp = Ineighbr(find(Ineighbr == min(Ineighbr(i,:)),1));
            Ineighbr(find(Ineighbr == min(Ineighbr(i,:)),1)) = min(Ineighbr(i,:))-1;
            Ineighbr(find(Ineighbr(:,2) == max(Ineighbr(i,:)),1),2) = temp(1);
        end
        while ~any(isnan(Ineighbr(i,:)))&& flagInf <3
            if max(Trial(Ineighbr(i,:))) < incorrect(i,1)% && ~(min(Ineighbr(i,:) == Ineighbr(i-1,:)))
                if(max(Ineighbr(i,:)))>= length(Trial)
                    Ineighbr(i,find(max(Ineighbr(i,:)) == Ineighbr(i,:),1)) = NaN;
                    Ineighbr(i,find(min(Ineighbr(i,:)) == Ineighbr(i,:),1)) = Ineighbr(i,find(min(Ineighbr(i,:)) == Ineighbr(i,:),1))+1;
                else
                    Ineighbr(find(min(Ineighbr == Ineighbr(i,:)),1)) = Ineighbr(find(min(Ineighbr == Ineighbr(i,:)))+1);
                end
            end
            flagInf = flagInf +1;
        end
        
        itInd(i,1) = Trial(nanmin(Ineighbr(i,:)));
        
        %Identifying Incorrect correction trials
        if any(itInd(i,1) == correction(:,1))
            criRow(i,1) = i;
            cri = cri + 1;
        end
        %setting trial length
        if any(isnan(Ineighbr(i,:)))
            itCount(i,1) = length(frameMap) - Trial(nanmax(Ineighbr(i,:)));
        else
            itCount(i,1) = abs(diff(Trial(Ineighbr(i,:))));
        end
        flagInf = 0;
    end
    %
    [Cneighbr, ic,~] = unique(Cneighbr,'rows','stable');
    [Ineighbr, ii, ~] = unique(Ineighbr,'rows','stable');
    ctCount = ctCount(ic);%unique(ctCount,'rows','stable');
    itCount = itCount(ii);%unique(itCount,'rows','stable');
    ctInd = unique(ctInd,'stable');
    itInd = unique(itInd,'stable');
    
    %Slicing Calcium trace into their respective pieces,correct(Ctrace)/incorrect(Itrace)/correction(CrTrace)
    % and taking out correction trials from Correct/Incorrect recordings
    if ~isempty(crcRow)
        %     crcRow = find(crcRow);
        [~,ia,~] = intersect(ic,crcRow);
        CrneighbrC = Cneighbr(ia,:);
        crCountC = ctCount(ia,:);
        crcInd = ctInd(ia);
        ctInd(ia,:)=  [];
        Cneighbr(ia,:)=  [];
        ctCount(ia,:)=  [];
    end
    if ~isempty(criRow)
        %     criRow = find(criRow);
        [~,ia,~] = intersect(ii,criRow);
        CrneighbrI = Ineighbr(ia,:);
        crCountI = itCount(ia,:);
        criInd = itInd(ia);
        itInd(ia,:)=  [];
        Ineighbr(ia,:)=  [];
        itCount(ia,:)=  [];
    else
        CrneighbrI = zeros(1,2);
        crCountI = [];
    end
    
    %sort delay periods
    sorted = 0;
    for i = 1: length(delay)
        for j = 1 : length(Cneighbr(:,1))
            if delay(i) < Trial(max(Cneighbr(j,:))) && delay(i) > Trial(min(Cneighbr(j,:)))|| (any(find(isnan(Cneighbr(j,:)))) && delay(i) > Trial(min(Cneighbr(j,:))))
                dc(dccount) = delay(i);
                dccount = dccount +1;
                sorted = 1;
                break
            end
        end
        if ~sorted
            for j = 1 : length(Ineighbr(:,1))
                if delay(i) < Trial(max(Ineighbr(j,:))) && delay(i) > Trial(min(Ineighbr(j,:)))|| (any(find(isnan(Ineighbr(j,:)))) && delay(i) > Trial(min(Ineighbr(j,:))))
                    di(dicount) = delay(i);
                    dicount = dicount +1;
                    sorted = 1;
                    break
                end
            end
        end
        if ~sorted
            for j = 1 : length(CrneighbrC(:,1))
                if (delay(i) < Trial(max(CrneighbrC(j,:))) && delay(i) > Trial(min(CrneighbrC(j,:)))) || (any(find(isnan(CrneighbrC(j,:)))) && delay(i) > Trial(min(CrneighbrC(j,:))))
                    dcc(dcccount) = delay(i);
                    dcccount = dcccount +1;
                    sorted = 1;
                    break
                end
            end
        end
        if ~sorted && ~isempty(criRow)
            for j = 1 : length(CrneighbrI(:,1))
                if delay(i) < Trial(max(CrneighbrI(j,:))) && delay(i) > Trial(min(CrneighbrI(j,:))) || (any(find(isnan(CrneighbrI(j,:)))) && delay(i) > Trial(min(CrneighbrI(j,:))))
                    dic(diccount) = delay(i);
                    diccount = diccount +1;
                    break
                end
            end
        end
        sorted = 0;
    end
    
    %Discecting the calcium traces in their respective trials
    if itteration == 1
        trace = transpose(normalize(ms.FiltTraces,'zscore'));
        binarize = Binarize(ms);
    else
        trace = transpose(normalize(msC.FiltTraces,'zscore'));        
        binarize = Binarize(msC);
    end
    btrace = transpose(binarize.binarizedTraces);
    % trace = transpose(normalize(ms.FiltTraces,'range'));                        %total trace
    Ctrace = nan(length(ms.FiltTraces(1,:)),max(ctCount),length(ctInd(:,1)));   %Correct Trace
    Itrace = nan(length(ms.FiltTraces(1,:)),max(itCount),length(itInd(:,1)));   %Incorrect Trace
    CRCtrace = nan(length(ms.FiltTraces(1,:)),max(crCountC),length(crcInd(:,1)));   %Correct correction trace
    if ~isempty(criRow)
        CRItrace = nan(length(ms.FiltTraces(1,:)),max(crCountI),length(criInd(:,1)));   %Incorrect correction trace
    end
    dctrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(dc));
    ditrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(di));
    dcctrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(dcc));
    dicctrace = nan(length(ms.FiltTraces(1,:)),dFrames+1,length(dic));
    
    for i = 1: length(ms.FiltTraces(1,:))
        %Front Anchored
        for j =1: length(Ctrace(1,1,:))
            Ctrace(i,1:ctCount(j),j) = trace(i,ctInd(j):ctInd(j)+ctCount(j)-1);
        end
        for j =1: length(Itrace(1,1,:))
            Itrace(i,1:itCount(j),j) = trace(i,itInd(j):itInd(j)+itCount(j)-1);
        end
        for j = 1: length(CRCtrace(1,1,:))
            CRCtrace(i,1:crCountC(j),j) = trace(i,crcInd(j):crcInd(j)+crCountC(j)-1);
        end
        if ~isempty(criRow)
            for j = 1: length(CRItrace(1,1,:))
                CRItrace(i,1:crCountI(j),j) = trace(i,criInd(j):criInd(j)+crCountI(j)-1);
            end
        else
            CRItrace = [];
        end
        %Delay
        for j =1: length(dc)
            dctrace(i,:,j) = trace(i,dc(j):dc(j)+dFrames); %Changed +60 to +dFrames
        end
        for j =1: length(di)
            ditrace(i,:,j) = trace(i,di(j):di(j)+dFrames);
        end
        for j =1: length(dcc)
            dcctrace(i,:,j) = trace(i,dcc(j):dcc(j)+dFrames);
        end
        for j =1: length(dic)
            dicctrace(i,:,j) = trace(i,dic(j):dic(j)+dFrames);
        end
    end                
    %Split-half reliability with random trial halves
    if itteration == 1
        out.DelaySplithalfcorrect = singleSplitShuffle(dctrace,0,1);
        out.DelaySplithalfincorrect = singleSplitShuffle(ditrace,0,1);
        out.DelaySplithalfccor = singleSplitShuffle(dcctrace,0,1);
        out.DelaySplithalficor = singleSplitShuffle(dicctrace,0,1);
        
%         out.dctrace = sortPeaks(dctrace);
%         out.ditrace = sortPeaks(ditrace);
%         out.dcctrace = sortPeaks(dcctrace);
%         out.dicctrace = sortPeaks(dicctrace);
        
        out.FrontSplithalfcorrect = singleSplitShuffle(Ctrace,ctCount,2);
        out.FrontSplithalfincorrect = singleSplitShuffle(Itrace,itCount,2);
        out.FrontSplithalfccor = singleSplitShuffle(CRCtrace,crCountC,2);
        out.FrontSplithalficor = singleSplitShuffle(CRItrace,crCountI,2);
        
%         out.frontCtrace = sortPeaks(Ctrace);
%         out.frontItrace = sortPeaks(Itrace);
%         out.frontCRCtrace = sortPeaks(CRCtrace);
%         out.frontCRItrace = sortPeaks(CRItrace);
    else
        out.ScorrDelaycorrect(:,itteration-1) = FullShuffleSplit(dctrace,0,1);
        out.ScorrDelayincorrect(:,itteration-1) = FullShuffleSplit(ditrace,0,1);
        out.ScorrDelayccor(:,itteration-1) = FullShuffleSplit(dcctrace,0,1);
        out.ScorrDelayicor(:,itteration-1) = FullShuffleSplit(dicctrace,0,1);
        
        out.ScorrFrontcorrect(:,itteration-1) = FullShuffleSplit(Ctrace,ctCount,2);
        out.ScorrFrontincorrect(:,itteration-1) = FullShuffleSplit(Itrace,itCount,2);
        out.ScorrFrontccor(:,itteration-1) = FullShuffleSplit(CRCtrace,crCountC,2);
        out.ScorrFronticor(:,itteration-1) = FullShuffleSplit(CRItrace,crCountI,2);
    end
    %Back Anchored
    bCtrace = nan(length(ms.FiltTraces(1,:)),max(ctCount),length(ctInd(:,1)));
    bItrace = nan(length(ms.FiltTraces(1,:)),max(itCount),length(itInd(:,1)));
    bCRCtrace = nan(length(ms.FiltTraces(1,:)),max(crCountC),length(crcInd(:,1)));
    if ~isempty(criRow)
        bCRItrace = nan(length(ms.FiltTraces(1,:)),max(crCountI),length(criInd(:,1)));
    end
    for i = 1: length(ms.FiltTraces(1,:))
        for j =1: length(bCtrace(1,1,:))
            bCtrace(i, (end- ctCount(j))+1:end,j) = trace(i,ctInd(j):ctInd(j)+ctCount(j)-1);
        end
        for j =1: length(bItrace(1,1,:))
            bItrace(i, (end- itCount(j))+1:end,j) = trace(i,itInd(j):itInd(j)+itCount(j)-1);
        end
        for j = 1: length(bCRCtrace(1,1,:))
            bCRCtrace(i,(end - crCountC(j))+1 :end,j) = trace(i,crcInd(j):crcInd(j)+crCountC(j)-1);
        end
        if ~isempty(criRow)
            for j = 1: length(bCRItrace(1,1,:))
                bCRItrace(i,(end-crCountI(j))+1 :end ,j) = trace(i,criInd(j):criInd(j)+crCountI(j)-1);
            end
        else
            bCRItrace = [];
        end
    end
    
    if itteration == 1 
        out.BackSplithalfcorrect = singleSplitShuffle(bCtrace,ctCount,3);
        out.BackSplithalfincorrect = singleSplitShuffle(bItrace,itCount,3);
        out.BackSplithalfccor = singleSplitShuffle(bCRCtrace,crCountC,3);
        out.BackSplithalficor = singleSplitShuffle(bCRItrace,crCountI,3);
%         out.backCtrace = sortPeaks(bCtrace);
%         out.backItrace = sortPeaks(bItrace);
%         out.backCRCtrace = sortPeaks(bCRCtrace);
%         out.backCRItrace = sortPeaks(bCRItrace);
    else
        out.ScorrBackcorrect(:,itteration-1) = FullShuffleSplit(bCtrace,ctCount,3);
        out.ScorrBackincorrect(:,itteration-1) = FullShuffleSplit(bItrace,itCount,3);
        out.ScorrBackccor(:,itteration-1) = FullShuffleSplit(bCRCtrace,crCountC,3);
        out.ScorrBackicor(:,itteration-1) = FullShuffleSplit(bCRItrace,crCountI,3);
    end
end
%Delay
out.Dpercentdc = prctile(out.ScorrDelaycorrect,99,2);
out.DcorrectNumberPassed = length(find(out.DelaySplithalfcorrect> out.Dpercentdc));
out.DcorrectCells = find(out.DelaySplithalfcorrect> out.Dpercentdc);
out.DcorrectPassed = dctrace(out.DcorrectCells,:,:);

out.Dpercentdi = prctile(out.ScorrDelayincorrect,99,2);
out.DincorrectNumberPassed = length(find(out.DelaySplithalfincorrect> out.Dpercentdi));
out.DincorrectCells = find(out.DelaySplithalfincorrect> out.Dpercentdi);
out.DincorrectPassed = ditrace(out.DincorrectCells,:,:);

out.Dpercentdcc = prctile(out.ScorrDelayccor,99,2);
out.DccorNumberPassed = length(find(out.DelaySplithalfccor> out.Dpercentdcc));
out.DccorCells = find(out.DelaySplithalfccor> out.Dpercentdcc);
out.DccorrPassed = dcctrace(out.DccorCells,:,:);

out.Dpercentdic = prctile(out.ScorrDelayicor,99,2);
out.DicorNumberPassed = length(find(out.DelaySplithalficor> out.Dpercentdic));
out.DicorCells = find(out.DelaySplithalficor> out.Dpercentdic);
out.DicorPassed = dicctrace(out.DicorCells,:,:);

%Front
out.Fpercentdc = prctile(out.ScorrFrontcorrect,99,2);
out.FcorrectNumberPassed = length(find(out.FrontSplithalfcorrect> out.Fpercentdc));
out.FcorrectCells = find(out.FrontSplithalfcorrect> out.Fpercentdc);
out.FcorrectPassed = Ctrace(out.FcorrectCells,:,:);

out.Fpercentdi = prctile(out.ScorrFrontincorrect,99,2);
out.FincorrectNumberPassed = length(find(out.FrontSplithalfincorrect> out.Fpercentdi));
out.FincorrectCells = find(out.FrontSplithalfincorrect> out.Fpercentdi);
out.FincorrectPassed = Itrace(out.FincorrectCells,:,:);

out.Fpercentdcc = prctile(out.ScorrFrontccor,99,2);
out.FccorNumberPassed = length(find(out.FrontSplithalfccor> out.Fpercentdcc));
out.FccorCells = find(out.FrontSplithalfccor> out.Fpercentdcc);
out.FccorPassed = CRCtrace(out.FccorCells,:,:);

out.Fpercentdic = prctile(out.ScorrFronticor,99,2);
out.FicorNumberPassed = length(find(out.FrontSplithalficor> out.Fpercentdic));
out.FicorCells = find(out.FrontSplithalficor> out.Fpercentdic);
out.FicorPassed = CRItrace(out.FicorCells,:,:);

%Back
out.Bpercentdc = prctile(out.ScorrBackcorrect,99,2);
out.BcorrectNumberPassed = length(find(out.BackSplithalfcorrect> out.Bpercentdc));
out.BcorrectCells = find(out.BackSplithalfcorrect> out.Bpercentdc);
out.BcorrectPassed = bCtrace(out.BcorrectCells,:,:);

out.Bpercentdi = prctile(out.ScorrBackincorrect,99,2);
out.BincorrectNumberPassed = length(find(out.BackSplithalfincorrect> out.Bpercentdi));
out.BincorrectCells = find(out.BackSplithalfincorrect> out.Bpercentdi);
out.BincorrectPassed = bItrace(out.BincorrectCells,:,:);

out.Bpercentdcc = prctile(out.ScorrBackccor,99,2);
out.BccorNumberPassed = length(find(out.BackSplithalfccor> out.Bpercentdcc));
out.BccorCells = find(out.BackSplithalfccor> out.Bpercentdcc);
out.BccorPassed = bCRCtrace(out.BccorCells,:,:);

out.Bpercentdic = prctile(out.ScorrBackicor,99,2);
out.BicorNumberPassed = length(find(out.BackSplithalficor> out.Bpercentdic));
out.BicorCells = find(out.BackSplithalficor> out.Bpercentdic);
out.BicorPassed = bCRItrace(out.BicorCells,:,:);

toc
end

function rast = sortPeaks(rast)
[~,maxind] = max(rast,[],2);
[~, rastsort] = sort(maxind);
rast = rast(rastsort,:);
end

function [out] = FullShuffleSplit(trial,cl,cut)
if ~isempty(trial)
    if ~(length(trial(1,1,:)) == 0)
        out = zeros(length(trial(:,1,1)),1);
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
                h1 = mean(temp((1:round(length(temp(:,1))/2)),:,i),1);
                h2 = mean(temp((round(length(temp(:,1))/2)+1:end),:,i),1);
                %         if ~isempty(find(isnan(h1),1))
                %             h1(find(isnan(h1))) = 0;
                %         elseif ~isempty(find(isnan(h2),1))
                %             h2(find(isnan(h2))) = 0;
                %         end
                out(i) = corr2(h1,h2);
            end
        end
    else
        out = NaN;
    end
else
    out = NaN;
end
end