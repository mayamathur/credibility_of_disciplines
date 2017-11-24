%%
% (C) Denes Szucs; 2016
% Please note that this code is provided as it is without any support.

%% INITIALIZE pdftoolbox
clear
clear java
javaaddpath('D:\matlab\pdfparse\PDFBox-0.7.3\lib\PDFBox-0.7.3.jar')
javaaddpath('D:\matlab\pdfparse\FontBox-0.1.0\lib\FontBox-0.1.0.jar')
%
pdfdoc = org.pdfbox.pdmodel.PDDocument;
reader = org.pdfbox.util.PDFTextStripper;


%%
% START SEARCH and read data from PDF folder

% {'t(' >> these are GENERAL psychology STYLE ; 't' >> Nature style
FindWhat     = { 'F(' 't(' };
FindExtender = { 'F'  't'  };
StrLength = 65; % length of string to be stored/parsed

cd 'D:\DATA\PDFsForAnalysis\pdfs'
Folder  = 'D:\DATA\PDFsForAnalysis\pdfs\';

% folders with pdf files
SubDirs = {'Psychological Science'  'Cognitive Psychology' 'Cognition' 'Acta Psychologica' 'JECP'  ...
           'Nature Neuroscience' 'Neuron' 'Brain' 'The Journal of Neuroscience' 'Cerebral Cortex'  ...
           'Neuroimage' 'Cortex' 'Biological Psychology' 'Neuropsychologia' 'Neuroscience' ...
           'Biological Psychiatry' 'J Psychiatric Research' 'Neurobiology of Aging' };

% initialize
clear Files FIndx TIndx BIndx RIndx ErrFiles FDir TDir BDir RDir
Fstr = repmat('@',[2*1E5 StrLength+1]); 
Fst1 = repmat('@',[2*1E5 StrLength+1]); 
Tstr = repmat('@',[2*1E5 StrLength+1]); 
FDir = repmat(0,[2*1E5 1]); FIndx = repmat(0,[2*1E5 1]);
TDir = repmat(0,[2*1E5 1]); TIndx = repmat(0,[2*1E5 1]);
Dirs = [];
Fnum = 1; Tnum = 1;
Errnum = 0; Encrynum = 1; SuccessFiles = 0;
% d=1 ; i=1; fi = 1; s = 1; % For testing
tic
for d = 1:length(SubDirs)
    DirList = dir(strcat(Folder,char(SubDirs(d)),'\','*.pdf'));        
    for i = 1:length(DirList)
        Dirs(d).FileNames(i) = cellstr(DirList(i).name);
        fprintf('%u%s%u%s%s%s%s\n',i,'/',...
                 length(DirList),' : File: ',char(SubDirs(d)) ,'  :  ', DirList(i).name);
        pdfdoc = pdfdoc.load(strcat(Folder,char(SubDirs(d)),'\',DirList(i).name));
        if pdfdoc.isEncrypted
            pdfdoc.getDocument().close;
            fprintf('%s%s\n','  ++++ > Encripted file  > ', DirList(i).name);
            EncryptedFiles(Encrynum).Name = DirList(i).name;
            Encrynum = Encrynum + 1;
        else
            try
                pdfstr = reader.getText(pdfdoc); % critical point; this gives java error for some files: TRY/CATCH
                pdfdoc.getDocument().close;
                pdfstr = deblank(char(pdfstr));

                SuccessFiles = SuccessFiles + 1;
                SuccessFileSeries(SuccessFiles) = SuccessFiles;
                %Files(SuccessFiles).Name = DirList(SuccessFiles).name;

                for fi = 1:length(FindWhat)
                    Where1 = strfind(pdfstr,char(FindWhat(fi)) ); % protected finds; F( , t(
                    Where2 = strfind(pdfstr,char(FindExtender(fi)) ); % extender: F , t
                    if ~isempty(Where1) | ~isempty(Where2)
                        % get rid of doubles in 2 relative to 1
                        % F( F(  F(
                        Doubles = find(ismember(Where2, Where1));
                        Where2(Doubles) = []; 
                        % get rid of too close ones protecting type 1
                        IdM = [repmat(1,[1 length(Where1)])  repmat(2,[1 length(Where2)])];
                        [WhereAll WInd] = sort([Where1 Where2]);
                        IdM = IdM(WInd);
                        WhereTooClose = 1;
                        while ~isempty(WhereTooClose)
                            WhereTooClose = find(diff(WhereAll)<20)+1;  % diff: y(2) - y(1); etc...
                            Protected = find(IdM(WhereTooClose)==1);
                            % next line is not good cos while N. line is protected, may delete N-1 which may also be protected!!
                            %WhereTooClose(Protected) = WhereTooClose(Protected)-1; % change to the one before                            
                            WhereTooClose(Protected) = [];
                            WhereAll(WhereTooClose)=[];
                        end
                        % check whether we have enough space after the last found items; if not, remove them
                        WhereAll(find((WhereAll+StrLength) > length(pdfstr))) = [];
                        
                        % analyze each string
                        for s = 1:length(WhereAll)
                            FoundStr = pdfstr(WhereAll(s):WhereAll(s)+StrLength);
                            FoundStr1 = FoundStr; % save orignal line for error checking if needed
                            % remove Enters, etc (just for dispay check)
                            FoundStr(FoundStr<14)=32;
                            % change strange occas. STRANGE signs in some journals to = sign : CORTEX
                            % it's good to do all this here to avoid huge array size
                            % BUT: beware of removing important stuff; e.g. p< 0.05
                            % can test if enabling  HelpFlag2 = 1;
                            
                            FoundStr(strfind(FoundStr,'¼')) = '='; % = in Cortex BUT > in Neuropsychologia! * ie. p value errors in Neuropsychologia must be double checked
                            
                            FoundStr = strrep(FoundStr,'H11005','=     '); % The Journal of Neuroscience; this comes at t(df)= & p=
                            FoundStr = strrep(FoundStr,'H11021','<     '); % The Journal of Neuroscience; this comes at p<
                            FoundStr = strrep(FoundStr,'H11022','>     '); % The Journal of Neuroscience; this comest at p > 0.05
                            FoundStr = strrep(FoundStr,'H11002','     -'); % The Journal of Neuroscience: minus
                            FoundStr = strrep(FoundStr,'H11006','     ±'); % The Journal of Neuroscience: plus/minus
                            FoundStr = strrep(FoundStr,'H9257  2'   ,'ETA2    '); % The Journal of Neuroscience
                            FoundStr = strrep(FoundStr,'H9257  p  2','ETAPAr2    '); % The Journal of Neuroscience
                            FoundStr = strrep(FoundStr,'H11003','     *'); % The Journal of Neuroscience: *
                            FoundStr = strrep(FoundStr,'C17','ETA');       % BRAIN
                            FoundStr = strrep(FoundStr,'p40','p=0');       % BRAIN
                            FoundStr = strrep(FoundStr,'p50','p=0');       % BRAIN
                            FoundStr = strrep(FoundStr,'= C0','=  -');     % BRAIN
                            FoundStr = strrep(FoundStr,'=C0','= -');       % BRAIN
                            FoundStr = strrep(FoundStr,'P5','p<');         % BRAIN; p<
                            FoundStr = strrep(FoundStr,'P  5','P  <');     % BRAIN; p<                            
                            FoundStr = strrep(FoundStr,'P4','p>');         % BRAIN; p>
                            
                            % remove spaces for checking 'p='
                            FoundStrNOSPACE = FoundStr;
                            FoundStrNOSPACE(FoundStrNOSPACE==32)=[]; % remove spaces for analysis
                            
                            % only analyze text involving = < > signs and 'p' char; this is likely to have p value                            
                            % this is necessary here otherwise F/Tstr are enormous
                            if  ( ...
                                (length(strfind(FoundStr(1:end),'='))>1)       | ...   % min TWO = signs
                                ( (length(strfind(FoundStr(1:end),'<'))>0) &     ...   % min 1 = AND 1 <
                                  (length(strfind(FoundStr(1:end),'='))>0)   ) | ...
                                ( (length(strfind(FoundStr(1:end),'>'))>0) &     ...   % min 1 = AND 1 > (for n.s.)
                                  (length(strfind(FoundStr(1:end),'='))>0)   )   ...
                                ) ...
                                & ...
                                ( ...
                                (~isempty(strfind(lower(FoundStrNOSPACE),'p='))) | ... % AND there is a p= ; this restricts very well for T TESTS, less relevant for F tests (there 'p' is OK)
                                (~isempty(strfind(lower(FoundStrNOSPACE),'p>'))) | ... % AND there is a p> ; this restricts very well for T TESTS, less relevant for F tests (there 'p' is OK)
                                (~isempty(strfind(lower(FoundStrNOSPACE),'p<'))) | ... % AND there is a p< ; this restricts very well for T TESTS, less relevant for F tests (there 'p' is OK)
                                (~isempty(strfind(lower(FoundStrNOSPACE),'n.s.'))) |  ...
                                (~isempty(strfind(lower(FoundStrNOSPACE),'ns')))      ...
                                 ) ; 
                                HelpFlag2 = 1;
                            else 
                                HelpFlag2 = 0;
                            end

                            %HelpFlag2 = 1; % if this is active; there's no =<>p check; all strings go thru
                            if fi==3 % 
                                HelpFlag2 = 1; % no filtering otherwise nothing will be found!
                            end
                            
                            % analysis runs if any of the above are true
                            if HelpFlag2==1
                                switch fi
                                    case 1 % F test
                                        Fstr(Fnum,:) = FoundStr;
                                        Fst1(Fnum,:) = FoundStr1;
                                        FDir(Fnum,1)   = d; % directory number
                                        FIndx(Fnum,1) = i;
                                        Fnum = Fnum + 1;
                                    case 2 % t test
                                        Tstr(Tnum,:) = FoundStr;
                                        TDir(Tnum,1)   = d; % directory number
                                        TIndx(Tnum,1) = i;
                                        Tnum = Tnum + 1;
                                end % switch
                            end % end of analysis if HeLPFLAG
                        end % for s
                    end % ~isempty
                end % fi
            catch MEsg
                pdfdoc.getDocument().close;
                Errnum = Errnum + 1;                
                fprintf('%s%s\n','  **** > FILE ERROR  > ', DirList(i).name);
                disp(MEsg.message(1:56));
                fprintf('%s%s\n','  **** > Error message truncated here < ', DirList(i).name);
                ErrFiles(Errnum).Name = DirList(i).name;
                ErrFiles(Errnum).Mesg = MEsg;
            end% try
        end % not encrypted
    end
    ParsedInEachFolder(d)  = length(DirList); % Total number of parsed files in each folder
    SuccessInEachFolder(d) = SuccessFiles; % No of successfully parsed files in each folder (consec)
end
toc
% Elapsed time is 3118.371983 seconds / 60 =  51.97 mins for ALL


%% Initial summary statistics

SuccessInEachFolder = [SuccessInEachFolder(1) diff(SuccessInEachFolder)];
DirInfoVAR.Folder = Folder;
DirInfoVAR.SubDirs = SubDirs;
DirInfoVAR.SuccessInEachFolder = SuccessInEachFolder;
DirInfoVAR.ParsedInEachFolder = ParsedInEachFolder;
DirInfoVAR.SuccessFiles = SuccessFiles;
DirInfoVAR.ParsedALL = sum(ParsedInEachFolder)

% delete extra length matrix data
Fstr(Fnum:end,:) = [];
Fst1(Fnum:end,:) = [];

Tstr(Tnum:end,:) = [];
FDir(Fnum:end,:) = [];
TDir(Tnum:end,:) = [];
% added: 03 Dec 2015
FIndx(Fnum:end,:) = [];
TIndx(Tnum:end,:) = [];

% stats
fprintf('\n%s%u\n','Files parsed: ',sum(ParsedInEachFolder));
fprintf('%s%u\n','Successful reads: ',SuccessFiles);
fprintf('%s\n','Successful (top) and total (bottom) reads in each folder and success %:');
disp(SubDirs)
disp(SuccessInEachFolder)
disp(ParsedInEachFolder)
disp(round((SuccessInEachFolder ./ ParsedInEachFolder)*100))
fprintf('%s%u\n%s%u\n%s%u\n','Initial F values (potential) = ',size(Fstr,1) , ...
                       'Initial T values (potential) = ',size(Tstr,1) , ...
                       'Initial r values (potential) = ',size(Rstr,1)  )

% ** SAVE POINT 1
% save 'RAW F,T,B,R data for ALL'  Tstr Fstr Bstr Rstr   FDir TDir BDir RDir   DirInfoVAR Dirs   FIndx TIndx BIndx RIndx
% check?
% uint16(Fst1(2,1:20))


%% Standardize formats for easier checking and making sure results are correct

% Fstr standardization: it is clearer to do here than in analysis part
% bring all to the standard reporting format, easier to check outcome: F (1,23) = xx
for i = 1:size(Fstr,1)
    a = Fstr(i,:);
    % correct strings with no parentheses
    [Fstart,Fend] = regexpi(a,'F\s*\d*,\d*\s*=');  % e.g. >> 'F  2,14   ='
    if ~isempty(Fstart)
        Fnd = a(Fstart(1):Fend(1));
        Fnd = strrep(Fnd,'F ','F(');
        Fnd = strrep(Fnd,' =',')=');
        Fstr(i,Fstart(1):Fend(1)) = Fnd;
    end
end
% The following pass only makes visual check easier but not necessary for analysis
FDelFlag = zeros([size(Fstr,1) 1]);
for i = 1:size(Fstr,1)
    a = Fstr(i,:); 
    a(a==32)=[];  % remove spaces to simplify search
    % strsplit(a,{'',''})
    % regexpi(a,'\d*','match') % all nums
    DF1 = regexpi(a,'F(\d*','match');
    if isempty(DF1)
        FDelFlag(i) = 1;
    end
end
for i = 1:size(Fstr,1)
    fprintf('%u %s %s %s \t%u \n', i,' : ',Fstr(i,:), '  <<>>   Del:',FDelFlag(i))
end
fprintf('%s%u%s%u\n','Initial Fstr size: ',size(Fstr,1),' ; To Delete:',length(find(FDelFlag)) );
% delete non-standard entries
Fstr(find(FDelFlag),:) = [];
FIndx(find(FDelFlag),:) = [];
FDir(find(FDelFlag),:) = [];
FDelFlag(find(FDelFlag),:) = [];
%FDelFlag = 0;


% ***
% Tstr standardization
% bring all to the standard reporting format, easier to check outcome: F (1,23) = xx
for i = 1:size(Tstr,1)
    a = Tstr(i,:);
    % correct strings with no parentheses
    [Tstart,Tend] = regexpi(a,'t\s*\d*\s*=');  % e.g. >> 'F  2,14   ='
    %regexpi(a,'t\s*\d*\s*=','match')
    if ~isempty(Tstart)
        Fnd = a(Tstart(1):Tend(1));
        Fnd = strrep(Fnd,'t ','t(');
        Fnd = strrep(Fnd,' =',')=');
        Tstr(i,Tstart(1):Tend(1)) = Fnd;
    end
end
% The following pass only makes visual checks easier but not necessary for analysis
TDelFlag = zeros([size(Tstr,1) 1]);
for i = 1:size(Tstr,1)
    a = Tstr(i,:); 
    a(a==32)=[];  % remove spaces to simplify search
    % strsplit(a,{'',''})
    % regexpi(a,'\d*','match') % all nums
    DF = regexpi(a,'t[(\d*.d*)]+=','match'); % parentheses are command characters for finding groups/tokens!!
    if ~isempty(DF)
        if isempty(sscanf(char(DF(1))','t(%f')) 
            TDelFlag(i) = 1;
        end
    else
        TDelFlag(i) = 1;
    end
end
for i = 1:size(Tstr,1)
    fprintf('%u %s %s %s %u \n', i,' : ',Tstr(i,:), '  <<>>   Del: ',TDelFlag(i))
end
fprintf('%s%u%s%u\n','Initial Tstr size: ',size(Tstr,1),' ; To Delete:',length(find(TDelFlag)) );
% delete non-standard entries
Tstr(find(TDelFlag),:) = [];
TIndx(find(TDelFlag),:) = [];
TDir(find(TDelFlag),:) = [];
TDelFlag(find(TDelFlag),:) = [];
%TDelFlag = 0;

