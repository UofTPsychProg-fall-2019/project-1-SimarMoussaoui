function TSWM6()
%%%In this version of the program, the positions of the display elements
%%%jitter form trial to trial to help eliminate trunk-based coordination
%%%across trials
%cID = computerID()

% In the case that you dont want to have the eyeTracker ON, please have it
% set to 0. That way, it can run without involving the eyetracker. Debug
% will be set to 1 so that only 5 trials will run.

eyeTrackerON = 0;
debug=1;

% Warning suggesting to include in case there can be display problems under
% the below condition
if ~eyeTrackerON
	disp('Warning: sync tests are skipped - this could produce display problems! Please delete the line below if this is a real experiment')
	Screen('Preference', 'SkipSyncTests', 1)
end

% I want to save each participant under a prefix of SMptcp (my initials and
% participant number), starting at 1.
ptcpPrefix = 'SMptcp';
ptcpNo = '1';

%%%Permitted eye drift from fixation before trial aborts (in visual
%%%degrees)
allowableFixation = 2;
allowableTargetFixation = 3;
fs = 44100; % Hz
t.mid = 0:1/fs:0.25; % seconds

t.short = 0:1/fs:0.075; % seconds
t.long = 0:1/fs:0.5; % seconds
f = 440; % Hz
f2 = 1760;
sounds.good = [sin(2.*pi.*f.*t.short),sin(2.*pi.*f2.*t.mid)];
sounds.noGood = [(sin(2.*pi.*f2.*t.short)).*0.5];


% Set my data directory in my C folder, under folder Simar, under folder
% Data. For checking purposes, please update the path for the data to be
% saved in.
dataDir = 'C:\Simar\Data';

% I want to state the last participant number, and make sure that the
% participant numbers are consistently going up.
ptcpDirInfo = GetPtcpDirInfo(dataDir, ptcpPrefix);
[pMax, pindexMax] = max(ptcpDirInfo.ptcpNo);
if pMax > 0
    disp(['the last participant was ' ptcpDirInfo.ptcpDirName{pindexMax}]);
end
% Below, gives the participant number and data folder
prompt = {'Participant #:', 'Data Folder:'};
defaults = {num2str(pMax+1), dataDir};
answer = inputdlg(prompt, 'Data Selection', 1, defaults);
ptcpNo = answer{1};
dataDir = answer{2};

ptcpID = [ptcpPrefix, ptcpNo];
ptcpPath = [dataDir,filesep,ptcpID];

theRand = rand; %%%Making sure that the random number generator is working
disp(theRand)
theRand = rand;
disp(theRand)

% Below are the prompts that I created for the researcher to indicate if
% the file exists, if some blocks completed, and what block condition it is
% (which may be redundent since for my experiment this time around, only
% one condition is used (Ego; 1). But for my future experiments which may
% implement more conditions, I included other conditions.
if exist(ptcpPath)
    str = ['Data folder ', ptcpPath, ' already exists. Some blocks already complete? Please enter block no. of the last block completed: '];
    startingBlock = input(str); %Displays above prompt, receives user input (y/n)
    str = ['Ok cool. Now please enter the block types that are left to complete: (1 = Ego, 2 = Allo, format = [1 2 2 2] ): '];
    blockPlan = input(str); %Displays above prompt, receives user input (y/n)
else
    str = ['Data folder ', ptcpPath, ' does not exist. Create? (Y/N) '];
    yesNo = input(str, 's'); %Displays above prompt, receives user input (y/n)
    if yesNo == 'y' | yesNo == 'Y'
        % directory doesn't exist; create it.
        mkdir(ptcpPath);
        startingBlock = 0;
        if theRand > .5
            blockPlan = [1 1 1]; %this would be different for the next experiment which would include more condition types
            if debug
                blockPlan = [1];
            end
            disp('The condition for this block is EGOCENTRIC')
            %disp('Screen condition comes first, three blocks')
        else
            blockPlan = [1 1 1];
            if debug
                blockPlan = [1];
            end
            disp('The condition for this block is EGOCENTRIC')
            % disp('Frame condition comes first, three blocks')
        end
    else
        disp('okay, bye!')
        return
        %%%%Exit program
    end
end
% the next thing to do in order to start the experiment is to press any key
% to start after one second of waiting.
disp('Today, you will view an array of squares and, after a saccade,')
disp('try to click on the memorized location of one of them.')
disp('In just a bit, the experimenter will give you more detailed instructions.')
disp('Press a key to start!')
pause(1);
KbWait

% Below will provide me with a log of eye offset from fixationpoints. I
% have set up the computers to allow me to watch the log from another
% computer to make sure offset is within acceptable range.

currentLogID = fopen([dataDir,filesep,'currentLog.txt'], 'w');
sessionLogID = fopen([ptcpPath,filesep,'log.txt'], 'w');

% I want to present the time it took, using code below:

fprintf(currentLogID,[num2str(clock),'\r\n'])
fprintf(sessionLogID,[num2str(clock),'\r\n'])

fprintf(currentLogID,[ptcpID, '\r\n\r\n'])
fprintf(sessionLogID,[ptcpID, '\r\n\r\n'])

warning off MATLAB:DeprecatedLogicalAPI5    %problem with KbCheck
warning off MATLAB:DeprecatedLogicalAPI

RubberBand = 0; %Refer to later lines of code...

% timing
fpT = 1;
fpRangeT = .1;  % FP presented for fpT +/- fpRangeT
barT = .4;  % .2;

if ~eyeTrackerON
    barT = .4;
end

maskT = .2; %  .2;

interFPdelayT = .5;

T1T = 0.05;

rand('state',sum(100*clock));

if eyeTrackerON
    % Initialization of the connection with the Gazetracker.
end

fname = 'SSD.edf';
fnameData = 'TSI.mat'; % this is the file name which consists of all my data that I want to save
%====== some parameters
eye_used = -1; % I am only tracting the left eye movements
el = [];


KbName('UnifyKeyNames')
pauseKey = 'p';

%parameters; how far the person is from the screen and the screen width
eyeScreenDist_m = 0.885;
ScreendWidth_m = 1.15;
if debug
    eyeScreenDist_m = 0.2;
    ScreendWidth_m = .35;
end

% Fixing display on the screen
% [window, screenRect]=Screen('OpenWindow', 0,[0,0,1920,1080]);
[window, screenRect]=Screen('OpenWindow', 0);
Screen('TextSize',window, 30);
% [window, screenRect]=Screen('OpenWindow', 0,[0,0,1920/2,1080/2]);

%Fixing the resolution of the screen
screenResH = screenRect(3)-screenRect(1);
screenResV = screenRect(4)-screenRect(2);

HC = screenResH/2; VC = screenResV/2;
ScreendWidth_deg = 2.* atan(ScreendWidth_m./2 ./eyeScreenDist_m) ./ pi .*180;
pixPERdeg = 2.*HC ./ ScreendWidth_deg;
cueTime = 0.5;
tISI = 0;

% saccade velocity
VelThreshDegPerSec = 30;
VelThreshPixPerSec = VelThreshDegPerSec .* pixPERdeg;
VelThreshPix = VelThreshPixPerSec ./ 500; % assuming that eye sampling = 500 Hz
VelThreshPix_squared = VelThreshPix.^2;

if eyeTrackerON
    
    [el] = calibrateEL(el,window,false);

end

el.simarDebug = true;

% setting up the blocks - for keeping track purposes
for iBlock = 1:length(blockPlan)
    thisBlock = blockPlan(iBlock);
    blockNo = startingBlock + iBlock;
	
    path = [ptcpPath, filesep, 'block',num2str(blockNo)];

    if exist(path, 'file') == 7
        error('the folder I''m trying to save to already exists! Shutting down to avoid overwriting data.')
    else
        mkdir(path);
    end
    
    saveAs = [path,filesep,fnameData]; %saving the data...

    %=== independent variables (some of which will be wrapped in)
    if thisBlock == 1 %%%Block type 1 = horizontal saccade
        FP_degx = [15 0];
        FP_degy = [0 15];
    elseif thisBlock == 2 %%%Block type 2 = vertical saccade
        FP_degx = [0 15];
        FP_degy = [15 0];
    else
        error('Invalid block number!')
    end
    nFp = length(FP_degx);
    % FP_degy = ones(1, nFp).*-5;
    
    nDir = 2; %The number of distractors

    %-- min FP distance for saccade detection
    FPdistDEG = 0.05 .* 15;      %based on pretest 10% yields strong SSD
    FPdistPIX = FPdistDEG .* pixPERdeg;
    FPdistPIX_squared = FPdistPIX.^2;

	
    targDir = [0 180];
    nXings = 8;
    intersectionSquareSize = 10; %%%visual degrees wide
    % intersectionSquareHeight = 3; %%%Fixation box center point (above screen center)
    intersectionSquareHeight = 0; %%%Fixation box center point (above screen center)
    % xing_degX = [linspace(-(intersectionSquareSize/2), (intersectionSquareSize/2), sqrt(nXings)) linspace(-(intersectionSquareSize/2), (intersectionSquareSize/2), sqrt(nXings)) linspace(-(intersectionSquareSize/2), (intersectionSquareSize/2), sqrt(nXings))]';
    % xing_degY = [xing_degX(1).*ones(1,3) xing_degX(2).*ones(1,3) xing_degX(3).*ones(1,3)]' + intersectionSquareHeight;
    xing_degX = [1.91; 4.62; -4.62; -1.91; -1.91; -4.62; 4.62; 1.91]';
    xing_degY = [4.62; 1.91; 1.91; 4.62; -4.62; -1.91; -1.91; -4.62]' + intersectionSquareHeight;

    xing.xmin = min(xing_degX)-((intersectionSquareSize/sqrt(nXings))/2);
    xing.xmax = max(xing_degX)+((intersectionSquareSize/sqrt(nXings))/2);
    xing.ymin = min(xing_degY)-((intersectionSquareSize/sqrt(nXings))/2);
    xing.ymax = max(xing_degY)+((intersectionSquareSize/sqrt(nXings))/2);
    
    %Below is for allocentric condition - not used in this experiment
    alloBox_degX = [min(xing_degX) max(xing_degX)].*2.2;
    alloBox_degY = [max(xing_degY) min(xing_degY)].*2.2;


    nRepeats = 3; %%%Changed from 4 to 3 
    %nRepeats = 4; 
    nTrials = nRepeats.*nXings.*nFp.*nDir; % calculated the number of trials based on conditions and variables...
    % Screen('CloseAll')
    
    %-- if the eye tracker is on, only present 5 trials.
    
    if eyeTrackerON == 0;
        maxTrials = 5;
        %maxTrials = 20;
    end
    
    if eyeTrackerON == 1;
        maxTrials = nTrials;
        %maxTrials = 90;
    end
    
    %--- build variable lists
    totalFPx = [];
    totalFPy = [];
    totalXingsX = [];
    totalXingsY = [];

    % Screen('closeall')
    for iRepeats = 1 : nRepeats
        for iXings = 1 : nXings %Xings means crossings - used from previous experiment but really just means location of the stimuli and distractors
            totalFPx = [totalFPx FP_degx -FP_degx];
            totalFPy = [totalFPy FP_degy -FP_degy];
            totalXingsX = [totalXingsX xing_degX(iXings).*ones(1,nFp.*nDir)];
            totalXingsY = [totalXingsY xing_degY(iXings).*ones(1,nFp.*nDir)];
        end
    end

    %-- convert to pixels
    totalFPx = round(totalFPx .* pixPERdeg) + HC;
    totalFPy = VC - round(totalFPy .* pixPERdeg);

    totalXingsX = round(totalXingsX .* pixPERdeg) + HC;
    totalXingsY = VC - round(totalXingsY .* pixPERdeg);
    
    alloBox_X = round(alloBox_degX .* pixPERdeg) + HC;
    alloBox_Y = VC - round(alloBox_degY .* pixPERdeg);

    %--standard display specs
    nScreen = 10;     % trial needs frames: blank - fp1 - fp2 - cue - fp+T1 - T2 - response

    %-- Open a graphics window on the main Screen

    white=WhiteIndex(window);
    black=BlackIndex(window);
    grey=(white+black)/2;
    inc=white-grey;
    colBack = black;
    colFP = white./1.5;
    % HideCursor;

    % FP specs
    fixL = 5;       % size FP
    targL = 5;
    lnwd = 2;



    % theData= zeros(nTrials, 11); % collects the data
    stopkey=KbName('space');
    goodTrialCount = 0;
    badTrialCount = 0;
    iTrial = 1;
    iFile = 1;
    trialRejected = false;
    HideCursor;
    keyResult.recalibrate = false;
	
    % setting conditions for each block...
	if iBlock >= 2
		prevBlock = blockPlan(iBlock-1)
		if thisBlock == prevBlock
			conditionChange = false;
			blockScreen(window, screenRect, grey, white, conditionChange)
		else
			conditionChange = true;
			blockScreen(window, screenRect, grey, white, conditionChange)
		end
	end
	
	theData = []
	
    while ~isempty(find(isnan(totalFPx)==0)) & (iTrial<=maxTrials) %run until trials used up or max # of trials
        theTimes.trialStart = GetSecs;


        if trialRejected
            %--flush variables and Screens - or memory explodes

                windowPtrs=Screen('Windows');
                for iPtr = 1:length(windowPtrs)
                    isOffscreen=Screen(windowPtrs(iPtr),'IsOffscreen');
                    if isOffscreen
                        Screen(windowPtrs(iPtr),'Close')
                    end
                end

            keyResult.recalibrate = false;

            badTrialCount = badTrialCount+1;
            sound(sounds.noGood,fs);
            fadeDuration = 1;
            gaussScaling = fadeDuration/4;
            fadeStart = GetSecs;
            badTrialCount = badTrialCount+1;
            
            % Here I am setting up when a participant needs a break or I
            % need to recalibrate. The P button will pause and give me the
            % options I can take, to recalibrate, quit... etc.
            
            while GetSecs < fadeStart+fadeDuration
                [pauseCheck] = myKBInput_PausePlease(pauseKey);
                if pauseCheck
                    [keyResult] = pauseScreen(window, screenRect, grey, white);
                    if keyResult.recalibrate
                        windowPtrs=Screen('Windows');
                        for iPtr = 1:length(windowPtrs)
                            isOffscreen=Screen(windowPtrs(iPtr),'IsOffscreen');
                            if isOffscreen
                                Screen(windowPtrs(iPtr),'Close')
                            end
                        end
                        [el] = calibrateEL(el,window,true);
                    end
                    trialRejected = true;
                    continue
                end
                currentLum = (normpdf((GetSecs-fadeStart)/gaussScaling))/normpdf(0)*white;
                Screen(window,'FillRect',[0,0,currentLum]);
                Screen('Flip', window, 0, 0, 2);
            end
            trialRejected = false;


            % Below is updating the log 
            fprintf(currentLogID,'Fixation error:\r\n')
            fprintf(currentLogID,[errorIn, '\r\n'])
            fprintf(currentLogID,['So far there have been ' num2str(goodTrialCount) ' good trials and ' num2str(badTrialCount) ' bad\r\n'])
            fprintf(sessionLogID, 'Fixation error:\r\n')
            fprintf(sessionLogID,[errorIn, '\r\n'])
            fprintf(sessionLogID,['Sof far there have been ' num2str(goodTrialCount) ' good trials and ' num2str(badTrialCount) ' bad\r\n'])
        end

        fclose(currentLogID);
        fclose(sessionLogID);
        currentLogID = fopen([dataDir,filesep,'currentLog.txt'], 'a');
        sessionLogID = fopen([path,filesep,'log.txt'], 'a');
        errorIn = [];
        
        % I want to give the participants a break half way through the
        % blocks, and to do that, I will present the below:
        if iTrial == round(nTrials/2)
            instString = ['Thanks, half-way there - let''s take a short break. Click to continue.'];
            Screen(window,'FillRect',grey);
            Screen(window,'DrawText',instString,50,50,white);
            Screen('Flip', window, 0, 0, 2);

            GetClicks;
        end
        if eyeTrackerON

            Eyelink('openfile', 'SSD.edf');

            %--- some of these commands need to be repeated, I believe
            Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
            Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');
            Eyelink('Command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
            Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,BLINK,SACCADE,BUTTON');
            Eyelink('Command', 'heuristic_filter = 0')
            if ScreendWidth_m > 1.5
                redctFact = 1/4; calR = round([0+redctFact.*screenResH 0+redctFact.*screenResV screenResH-redctFact.*screenResH screenResV-redctFact.*screenResV]);
                eval(['Eyelink(''Command'', ''Screen_pixel_coords = ', num2str(calR), ''')']);
            end

            % STEP 5
            % start recording eye position
            Eyelink('StartRecording');
            % record a few samples before we actually start displaying
            WaitSecs(0.1);
            % mark zero-plot time in data file
            Eyelink('Message', 'SYNCTIME');
        end
        
        % As mentioned in the beginning of the code, I am implementing
        % jitter (codes below)
        stimJitter = (rand-.5 )*150;
        stimJitterV = (rand-.5 )*150;

        % pick trial
        selTrial = ceil(rand.*nTrials); %Randi(nTrials);
        while isnan(totalFPx(selTrial))
            selTrial = ceil(rand.*nTrials); %Randi(nTrials);
        end
        
        % add jitter to the trials
        selFPx = totalFPx(selTrial)+stimJitter;
        selFPy = totalFPy(selTrial)+stimJitterV;
        selXingsX = totalXingsX(selTrial)+stimJitter;
        selXingsY = totalXingsY(selTrial)+ stimJitterV;
        
        thisalloBox_X = alloBox_X +stimJitter;
        thisalloBox_Y = alloBox_Y + stimJitterV;

        %-- pick distractor intersection (location) 1
        iPickDistr = ceil(rand.*length(totalXingsX)); %Randi(length(totalXingsX));
        selXingsXd1 = totalXingsX(iPickDistr)+stimJitter;
        selXingsYd1 = totalXingsY(iPickDistr) + stimJitterV;
        while (selXingsXd1 == selXingsX) & (selXingsYd1 == selXingsY)
            iPickDistr = ceil(rand.*length(totalXingsX)); %Randi(length(totalXingsX));
            selXingsXd1 = totalXingsX(iPickDistr)+stimJitter;
            selXingsYd1 = totalXingsY(iPickDistr)+ stimJitterV;
        end

        %-- pick distractor intersection (location) 2
        iPickDistr = ceil(rand.*length(totalXingsX)); %Randi(length(totalXingsX));
        selXingsXd2 = totalXingsX(iPickDistr)+stimJitter;
        selXingsYd2 = totalXingsY(iPickDistr) + stimJitterV;
        while ((selXingsXd2 == selXingsX) & (selXingsYd2 == selXingsY)) | ((selXingsXd2 == selXingsXd1) & (selXingsYd2 == selXingsYd1))
            iPickDistr = ceil(rand.*length(totalXingsX)); %Randi(length(totalXingsX));
            selXingsXd2 = totalXingsX(iPickDistr)+stimJitter;
            selXingsYd2 = totalXingsY(iPickDistr) + stimJitterV;
		end
		
		%-- pick distractor intersection 3 - NOT Needed for my experiment
		%but including it just in case!
        iPickDistr = ceil(rand.*length(totalXingsX)); %Randi(length(totalXingsX));
        selXingsXd3 = totalXingsX(iPickDistr)+stimJitter;
        selXingsYd3 = totalXingsY(iPickDistr) + stimJitterV;
        while ((selXingsXd3 == selXingsX) & (selXingsYd3 == selXingsY)) | ((selXingsXd3 == selXingsXd1) & (selXingsYd3 == selXingsYd1)) | ((selXingsXd3 == selXingsXd2) & (selXingsYd3 == selXingsYd2))
            iPickDistr = ceil(rand.*length(totalXingsX)); %Randi(length(totalXingsX));
            selXingsXd3 = totalXingsX(iPickDistr)+stimJitter;
            selXingsYd3 = totalXingsY(iPickDistr) + stimJitterV;
		end
		
        % create position squares for stimuli
        fpx_1 = HC + stimJitter;
        fpy_1 = VC + stimJitterV;
        fpx_2 = selFPx;
        fpy_2 = selFPy;
        theRectFP = [fpx_1 fpy_1 fpx_1 fpy_1] + [-1 -1 1 1].*(fixL./2);
        theRectFP2 = [fpx_2 fpy_2 fpx_2 fpy_2] + [-1 -1 1 1].*(fixL./2);
        targDims.x = [-10,10];
        targDims.y = [-10,10];


        %============= prepare off-Screen windows: blank - fp1 - fp2 - fp+T1 - T2 - response
        for iScreen = 1:nScreen
            offwin(iScreen)=Screen(window,'OpenOffScreenWindow');
        end
        iScreen = 0;

        %--blank frame
        iScreen = iScreen + 1; iBlank = iScreen;
        Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);

        %--fp1 (fixation point 1) frame
        iScreen = iScreen + 1; iFPDim = iScreen;
        Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
        Screen(offwin(iScreen),'FillOval',colFP.*0.75,theRectFP);

        %--fp1 bright
        iScreen = iScreen + 1; iFP = iScreen;
        Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
        Screen(offwin(iScreen),'FillOval',colFP,theRectFP);

        %--fp2 frame
        iScreen = iScreen + 1; iFP2 = iScreen;
        Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
        Screen(offwin(iScreen),'FillOval',colFP,theRectFP2);

        %--fp1 + bar1 frame NOT FOR MY EXPERIMENT
        iScreen = iScreen + 1; iFPbar = iScreen;
        Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
        Screen(offwin(iScreen),'FillOval',colFP,theRectFP);
        fromH = 0;
        fromV = screenResV;
        toH = screenResH/2;
        toV = 0;

        %pickColr = Randi(3);
		
		colourOptions{1} = [0 white 0]./1.5;
		colourOptions{2} = [0 0 white]./1.2;
		colourOptions{3} = [white 0 0]./1.5;
		colourOptions{4} = [white.*2./3 white./3 0];     %orange
        colourOptions{5} = [white white 0]./1.6;          %yellow
        colourOptions{6} = [white 0 white]./1.6;          %purple
        colourOptions{7} = [0 white white]./1.7;          %cyan
        
%         colourOptions{1} = white.*[1 .3 .3];
% 		colourOptions{2} = white.*[.1 .9 .1];
% 		colourOptions{3} = white.*[.3 .7 1];
% 		colourOptions{4} = white.*[1 1 0];
        
        % colours for Jiali's EEG experiment
        % [1 .3 .3; .1 .9 .1; .3 .7 1; 1 1 0]
		
		data = colourOptions;
		numColr = length(data);
		rndIDX = randperm(numColr);
		
		pickColr = rndIDX(1);
		pickD1Colr = rndIDX(2);
		pickD2Colr = rndIDX(3);
		pickD3Colr = rndIDX(4);
		
		colBar = colourOptions{pickColr};
		colBard1 = colourOptions{pickD1Colr};
		colBard2 = colourOptions{pickD2Colr};
		colBard3 = colourOptions{pickD3Colr};
		
		
%         if pickColr == 1
%             colBar =   [white 0 0]./1.5;
%             pickColrDist = Randi(2);
%             if pickColrDist == 1
%                 colBard2 = [0 white 0]./1.5;
%                 colBard1 = [0 0 white]./1.5;
%             else
%                 colBard1 = [0 white 0]./1.5;
%                 colBard2 = [0 0 white]./1.5;
%             end
%         elseif pickColr == 2
%             colBar =   [0 white 0]./1.5;
%             pickColrDist = Randi(2);
%             if pickColrDist == 1
%                 colBard1 = [white 0 0]./1.5;
%                 colBard2 = [0 0 white]./1.5;
%             else
%                 colBard2 = [white 0 0]./1.5;
%                 colBard1 = [0 0 white]./1.5;
%             end
%         else
%             colBar =   [0 0 white]./1.5;
%             pickColrDist = Randi(2);
%             if pickColrDist == 1
%                 colBard1 = [white 0 0]./1.5;
%                 colBard2 = [0 white 0]./1.5;
%             else
%                 colBard2 = [white 0 0]./1.5;
%                 colBard1 = [0 white 0]./1.5;
%             end
%         end

        %-- distractor bar - presenting my distractors
        nDistr = 2;
        if nDistr > 0
            % 		Screen(offwin(iScreen),'DrawLine', colBard1, flashbarX1Leftd1,flashbarY1Leftd1,flashbarX1Rightd1,flashbarY1Rightd1)
            Screen(offwin(iScreen),'FillRect', colBard1, [selXingsXd1+targDims.x(1),selXingsYd1+targDims.y(1),selXingsXd1+targDims.x(2),selXingsYd1+targDims.y(2)])
            if nDistr > 1
                % 			Screen(offwin(iScreen),'DrawLine', colBard2, flashbarX1Leftd2,flashbarY1Leftd2,flashbarX1Rightd2,flashbarY1Rightd2)
                Screen(offwin(iScreen),'FillRect', colBard2, [selXingsXd2+targDims.x(1),selXingsYd2+targDims.y(1),selXingsXd2+targDims.x(2),selXingsYd2+targDims.y(2)])
            else
                selXingsXd2 = NaN; % if I didnt want more than 1 distractor - NOT THE CASE IN THIS EXPERIMENT
                selXingsYd2 = NaN;
			end
			
			if nDistr > 2
				Screen(offwin(iScreen),'FillRect', colBard3, [selXingsXd3+targDims.x(1),selXingsYd3+targDims.y(1),selXingsXd3+targDims.x(2),selXingsYd3+targDims.y(2)])
			else
				selXingsXd3 = NaN; % if I didnt want more than 2 distractors - NOT THE CASE IN THIS EXPERIMENT
				selXingsYd3 = NaN;
			end
        else
            selXingsXd1 = NaN; % if I didnt want any distractors - NOT THE CASE IN THIS EXPERIMENT
            selXingsYd1 = NaN;
            selXingsYd2 = NaN;
            selXingsXd2 = NaN;
			selXingsYd3 = NaN;
            selXingsXd3 = NaN;
        end
        %-- experimental bar
        Screen(offwin(iScreen),'FillRect', colBar, [selXingsX+targDims.x(1),selXingsY+targDims.y(1),selXingsX+targDims.x(2),selXingsY+targDims.y(2)])


        %--fp2 + bar2 frame
        iScreen = iScreen + 1; iFPRecall = iScreen;
        Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
        fromH = screenResH; fromV = screenResV; toH = screenResH/2; toV = 0;
		
		if thisBlock == 1 %%%Block type 1 is without mem array relocation
			shakeBox.width = 0;
			shakeBox.height = 0;
		elseif thisBlock == 2 %%%Block type 1 has mem array relocation (allocentric condition NOT IN MY CURRENT EXPERIMENT - FOR FUTURE)
			shakeBox.width = 0;
			shakeBox.height = 0;
		else
			error('Invalid block type!')
		end
% 		ShakeBox.topLeftOffsetX = screenResH/2 - shakeBox.width
% 		ShakeBox.topLeftOffsetY = (screenResV - shakeBox.height)/2
% 		
		shakeBox.offsetX = (rand*(shakeBox.width))-(shakeBox.width/2)
		shakeBox.offsetY = (rand*(shakeBox.height))-(shakeBox.height/2)
		
		

        if nDistr > 0
            % 		Screen(offwin(iScreen),'DrawLine', colBard1, flashbarX1Leftd1,flashbarY1Leftd1,flashbarX1Rightd1,flashbarY1Rightd1)
			if thisBlock == 1
				%%%Don't draw the recall anchor if it's not the allocentric
				%%%condition
				memAnchorDims.x1 = NaN;
				memAnchorDims.y1 = NaN;
				memAnchorDims.x2 = NaN;
				memAnchorDims.y2 = NaN;
				
				memAnchorCentersX = NaN;
				memAnchorCentersY = NaN;
			elseif thisBlock == 2
				%Screen(offwin(iScreen),'FillRect', colBard1, [shakeBox.offsetX+selXingsXd1+targDims.x(1),shakeBox.offsetY+selXingsYd1+targDims.y(1),shakeBox.offsetX+selXingsXd1+targDims.x(2),shakeBox.offsetY+selXingsYd1+targDims.y(2)])
                %Screen(offwin(iScreen),'FillRect', colFP, [shakeBox.offsetX+thisalloBox_X(1),shakeBox.offsetY+thisalloBox_Y(1),shakeBox.offsetX+thisalloBox_X(2),shakeBox.offsetY+thisalloBox_Y(2)])
                %Screen(offwin(iScreen),'FillRect', colBack, [shakeBox.offsetX+thisalloBox_X(1)+lnwd,shakeBox.offsetY+thisalloBox_Y(1)+lnwd,shakeBox.offsetX+thisalloBox_X(2)-lnwd,shakeBox.offsetY+thisalloBox_Y(2)-lnwd])

                memAnchorDims.x1 = NaN;
				memAnchorDims.y1 = NaN;
				memAnchorDims.x2 = NaN;
				memAnchorDims.y2 = NaN;
				
				memAnchorCentersX = NaN;
				memAnchorCentersY = NaN;
            end
            
            %Screen(offwin(iScreen),'FillOval',colFP,theRectFP2);

%             if nDistr > 1
%                 % 			Screen(offwin(iScreen),'DrawLine', colBard2, flashbarX1Leftd2,flashbarY1Leftd2,flashbarX1Rightd2,flashbarY1Rightd2)
% %                 Screen(offwin(iScreen),'FillRect', colBard2, [shakeBox.offsetX+selXingsXd2+targDims.x(1),shakeBox.offsetY+selXingsYd2+targDims.y(1),shakeBox.offsetX+selXingsXd2+targDims.x(2),shakeBox.offsetY+selXingsYd2+targDims.y(2)])
% 
%             else
%                 selXingsXd2 = NaN;
%             end
		else
            selXingsXd1 = NaN;
            selXingsYd1 = NaN;
            selXingsYd2 = NaN;
            selXingsXd2 = NaN;
			selXingsYd3 = NaN;
            selXingsXd3 = NaN;
        end
		
% 		Screen(offwin(iScreen),'FillRect', colBar, [shakeBox.offsetX+selXingsX+targDims.x(1),shakeBox.offsetY+selXingsY+targDims.y(1),shakeBox.offsetX+selXingsX+targDims.x(2),shakeBox.offsetY+selXingsY+targDims.y(2)])
        
        % BELOW CAN BE IGNORED - NOT USING THIS FOR CURRENT EXPERIMENT
		memTargetDims.x1 = shakeBox.offsetX+selXingsX+targDims.x(1);
		memTargetDims.y1 = shakeBox.offsetY+selXingsY+targDims.y(1);
		memTargetDims.x2 = shakeBox.offsetX+selXingsX+targDims.x(2);
		memTargetDims.y2 = shakeBox.offsetY+selXingsY+targDims.y(2);
		
		memTargetCentersX = shakeBox.offsetX+selXingsX;
		memTargetCentersY = shakeBox.offsetY+selXingsY;

        %--mask 1 frame
        iScreen = iScreen + 1; iMask1 = iScreen;
        Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
        Screen(offwin(iScreen),'FillOval',colFP,theRectFP);
		%nLines = 40;
		nLines = 1500;
		
		for iLine = 1 : nLines
			pickColr = ceil(rand.*numColr); %ceil(rand.*numColr)
            colMask = colourOptions{pickColr};
			fromH = ceil(rand.*screenResH); %Randi(screenResH);
            fromV = ceil(rand.*screenResV); %randi(screenResV);
            Screen(offwin(iScreen),'FillRect', colMask, [fromH+targDims.x(1),fromV+targDims.y(1),fromH+targDims.x(2),fromV+targDims.y(2)])
        end

		%--mask 2 frame
		iScreen = iScreen + 1;
		iMask2 = iScreen;
		Screen(offwin(iScreen),'FillRect',colBack, [0 0 screenResH screenResV]);
		Screen(offwin(iScreen),'FillOval',colFP,theRectFP2);
		for iLine = 1 : nLines
			pickColr = ceil(rand.*numColr); %ceil(rand.*numColr)
            colMask = colourOptions{pickColr};
			fromH = ceil(rand.*screenResH); %Randi(screenResH);
            fromV = ceil(rand.*screenResV); %randi(screenResV);
			Screen(offwin(iScreen),'FillRect', colMask, [fromH+targDims.x(1),fromV+targDims.y(1),fromH+targDims.x(2),fromV+targDims.y(2)])
        end

        %================ bring off-Screens to front
        if eyeTrackerON
            % Check recording status, stop display if error
            errorEL=Eyelink('CheckRecording');
            if (errorEL~=0), Screen('closeall'); Eyelink('Shutdown'); 'eyelink connection not working'; return; end;

            %- first sample eye position

            % STEP 6 - monitor eye
            if eye_used == -1 % if we don't know which eye to use, first find eye that's being tracked
                eye_used = Eyelink( 'EyeAvailable'); % get eye that's tracked
                if eye_used == el.BINOCULAR; % if both eyes are tracked
                    eye_used = el.LEFT_EYE; % use left eye
                end
            end

        end

        %-- frame 1: FP1 ------------------------------------------------
        [x,y,keyIsDown1] = GetMouse;
        keyIsDown1 = [0 0 0];
        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iFPDim),window);
        Screen('Flip', window, 0, 0, 2);
        if eyeTrackerON
            Eyelink('Message', 'FP1');
        end

        wasZero = false;
        while 1 && ~trialRejected

            % 		[pauseCheck] = myKBInput_PausePlease(pauseKey);

            pauseCheck = false;
            [kbKeyDown, responseClock, keyCode] = KbCheck;
            if kbKeyDown
                KBkeys = KbName(keyCode);
                if iscell(KBkeys)
                    KBkey = KBkeys{1};
                else
                    KBkey = KBkeys;
                end

                if strcmp(KBkey,pauseKey)
                    pauseCheck = true;
                else
                end
            else
            end

            if pauseCheck
                [keyResult] = pauseScreen(window, screenRect, grey, white);
                if keyResult.recalibrate
                    windowPtrs=Screen('Windows');
                    for iPtr = 1:length(windowPtrs)
                        isOffscreen=Screen(windowPtrs(iPtr),'IsOffscreen');
                        if isOffscreen
                            Screen(windowPtrs(iPtr),'Close')
                        end
                    end
                    [el] = calibrateEL(el,window,true);
                end
                trialRejected = true;
                continue
            end

            [notUsed,notUsed,keyIsDown1] = GetMouse;
            if keyIsDown1(1) == 0
                wasZero = true;
            end
            if wasZero == true && keyIsDown1(1) == 1
                break
            end
        end
        if trialRejected == true;
            continue
        end
        bufferSize = 10;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = 0.1;
        allowableOffset = 20;
        fixationOffset.x = 0;
        fixationOffset.y = 0;
        realFixation.x = 0;
        realFixation.y = 0;

        [notUsed, realFixation, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation,allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);

        fprintf(currentLogID,notes);
        fprintf(sessionLogID,notes);

        fixationOffset.x = fixLocation.x - realFixation.x;
        fixationOffset.y = fixLocation.y - realFixation.y;

        fprintf(currentLogID,['Fixation offset = ' num2str(sqrt((fixationOffset.x^2)+(fixationOffset.y^2))/pixPERdeg) '\r\n']);
        fprintf(sessionLogID,['Fixation offset = ' num2str(sqrt((fixationOffset.x^2)+(fixationOffset.y^2))/pixPERdeg) '\r\n']);


        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iFP),window);
        Screen('Flip', window, 0, 0, 2);

        bufferSize = 3;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = (fpT+ 2.*(rand-.5).*fpRangeT);
        allowableOffset = allowableTargetFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        fprintf(currentLogID,notes);
        fprintf(sessionLogID,notes);

        if trialRejected
            errorIn = 'Initial fixation';
            continue
        end


        %-- frame 2: FP1 + mamory array ------------------------------------------------
        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iFPbar),window);
        Screen('Flip', window, 0, 0, 2);
        if eyeTrackerON
            Eyelink('Message', 'bar');
        end

        bufferSize = 3;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = barT;
        allowableOffset = allowableFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        fprintf(currentLogID,notes);
        fprintf(sessionLogID,notes);


        if trialRejected
            errorIn = 'Memory encoding';
            continue
        end


        %-- frame 3: mask 1 ------------------------------------------------
        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iMask1),window);
        Screen('Flip', window, 0, 0, 2);
        if eyeTrackerON
            Eyelink('Message', 'mask1');
        end

        bufferSize = 3;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = maskT;
        allowableOffset = allowableFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        fprintf(currentLogID,notes);
        fprintf(sessionLogID,notes);
        if trialRejected
            errorIn = 'First mask';
            continue
        end



        %-- again frame 1: FP1 ------------------------------------------------
        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iFP),window);
        Screen('Flip', window, 0, 0, 2);
        if eyeTrackerON
            Eyelink('Message', 'FP1');
        end

        bufferSize = 3;
        fixLocation.x = fpx_1;
        fixLocation.y = fpy_1;
        holdTime = interFPdelayT;
        allowableOffset = allowableFixation;
        [trialRejected, notUsed, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation, allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg);
        fprintf(currentLogID,notes);
        fprintf(sessionLogID,notes);
        if trialRejected
            errorIn = 'Just before saccade';
            continue
        end


        %-- frame 4: FP2 ------------------------------------------------
        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iFP2),window);
        Screen('Flip', window, 0, 0, 2);
        if eyeTrackerON
            Eyelink('Message', 'FP2');
        end
        WaitSecs(fpT+ 2.*(rand-.5) .* fpRangeT)

          
        %-- frame 4=6: FP2 ------------------------------------------------
        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iFP2),window);
        Screen('Flip', window, 0, 0, 2);
        if eyeTrackerON
            Eyelink('Message', 'FP2');
        end

        WaitSecs(fpT+ 2.*(rand-.5) .* fpRangeT)

		
		%-- frame 4=6: FP2 ------------------------------------------------
        %Screen(window,'WaitBlanking');
		Screen('CopyWindow',offwin(iFPRecall),window);
        Screen('Flip', window, 0, 0, 2);
        if eyeTrackerON
            Eyelink('Message', 'RecallStim');
        end
		
%         mouseLoc.x = round((((rand*(xing.xmax-xing.xmin))+xing.xmin)*pixPERdeg)+HC);
%         mouseLoc.y = VC - round((((rand*(xing.ymax-xing.ymin))+xing.ymin)*pixPERdeg));
        
        mouseLoc.x = round(fpx_2);
        mouseLoc.y = round(fpy_2);
        
        SetMouse(mouseLoc.x, mouseLoc.y);
        x = HC; y = VC;
        while ~trialRejected %~keyIsDown(1)
            [pauseCheck] = myKBInput_PausePlease(pauseKey);
            if pauseCheck
                [keyResult] = pauseScreen(window, screenRect, grey, white);
                if keyResult.recalibrate
                    windowPtrs=Screen('Windows');
                    for iPtr = 1:length(windowPtrs)
                        isOffscreen=Screen(windowPtrs(iPtr),'IsOffscreen');
                        if isOffscreen
                            Screen(windowPtrs(iPtr),'Close')
                        end
                    end
                    [el] = calibrateEL(el,window,true);
                end
                trialRejected = true;
                continue
            end
            %Screen(window,'WaitBlanking');
            %     Screen(window,'FillOval',colBack,[x y x y] + [-1 -1 1 1].*(targL./2));
            % 	colFP
            [x,y,keyIsDown] = GetMouse; % update cursor position

            if RubberBand == 0
                Screen('CopyWindow',offwin(iFP2),window);
				Screen('CopyWindow',offwin(iFPRecall),window);
                Screen(window,'FillRect',colBar,[x y x y] + [targDims.x(1) targDims.y(1) targDims.x(2) targDims.y(2)]);
                Screen(window,'FillRect',colBack,[x y x y] + [targDims.x(1)+lnwd targDims.y(1)+lnwd targDims.x(2)-lnwd targDims.y(2)-lnwd]);
				mouseDims.X1 = x+targDims.x(1); 
				mouseDims.Y1 = y+targDims.y(1);
				mouseDims.X2 = x+targDims.x(2);
				mouseDims.Y2 = y+targDims.y(2);
            else
                if x < HC
                    x1 = x;
                    x2 = HC;
                else
                    x1 = HC;
                    x2 = x;
                end
                Screen('CopyWindow',offwin(iFP2),window);
				Screen(window,'FillRect',colBar,[x y x y] + [targDims.x(1) targDims.y(1) targDims.x(2) targDims.y(2)].*(fixL./2));
				mouseDims.X1 = x+targDims.x(1);
				mouseDims.Y1 = y+targDims.y(1);
				mouseDims.X2 = x+targDims.x(2);
				mouseDims.Y2 = y+targDims.y(2);
            end
            Screen('Flip', window, 0, 0, 2);
            if keyIsDown(1)
                responseX = x;
                responseY = y;
                break
            end
        end
        
        if trialRejected == true;
            continue
        end

        if eyeTrackerON
            Eyelink('Message', 'response');
        end

        HideCursor;
        WaitSecs(tISI);

        %--store data
        %   theData(iFile,:) = [iTrial, selFPx, selFPy, selXingsX, selXingsY, responseX, responseY, selXingsXd1, selXingsYd1, selXingsXd2, selXingsYd2];
        theData.iTrial(iFile) = iTrial;
        theData.selFPx(iFile) = selFPx;
        theData.selFPy(iFile) = selFPy;
        theData.FPx(iFile) = fpx_1;
        theData.FPy(iFile) = fpy_1;
        theData.FPx2(iFile) = fpx_2;
        theData.FPy2(iFile) = fpy_2;
        theData.selXingsX(iFile) = selXingsX;
        theData.selXingsY(iFile) = selXingsY;
        theData.responseX(iFile) = responseX;
        theData.responseY(iFile) = responseY;
        theData.selXingsXd1(iFile) = selXingsXd1;
        theData.selXingsYd1(iFile) = selXingsYd1;
        theData.selXingsXd2(iFile) = selXingsXd2;
        theData.selXingsYd2(iFile) = selXingsYd2;
		theData.selXingsXd3(iFile) = selXingsXd3;
        theData.selXingsYd3(iFile) = selXingsYd3;
        theData.mousePosInit.x(iFile) = mouseLoc.x;
        theData.mousePosInit.y(iFile) = mouseLoc.y;
        theData.stimJitterX(iFile) = stimJitter;
        theData.stimJitterY(iFile) = stimJitterV;
		theData.blockNo(iFile) = blockNo;
		theData.blockType(iFile) = thisBlock;
		
		theData.memTargetDims.X1(iFile) = memTargetDims.x1;
		theData.memTargetDims.Y1(iFile) = memTargetDims.y1;
		theData.memTargetDims.X2(iFile) = memTargetDims.x2;
		theData.memTargetDims.Y2(iFile) = memTargetDims.y2;
	
		theData.memAnchorDims.X1(iFile) = memAnchorDims.x1;
		theData.memAnchorDims.Y1(iFile) = memAnchorDims.y1;
		theData.memAnchorDims.X2(iFile) = memAnchorDims.x2;
		theData.memAnchorDims.Y2(iFile) = memAnchorDims.y2;
		
		theData.memTargetCentersX(iFile) = memTargetCentersX;
		theData.memTargetCentersY(iFile) = memTargetCentersY;

		theData.memAnchorCentersX(iFile) = memAnchorCentersX;
		theData.memAnchorCentersY(iFile) = memAnchorCentersY;

		
		theData.mouseDims.X1(iFile) = mouseDims.X1;
		theData.mouseDims.Y1(iFile) = mouseDims.Y1;
		theData.mouseDims.X2(iFile) = mouseDims.X2;
		theData.mouseDims.Y2(iFile) = mouseDims.Y2;
		
		theData.colBar{iFile} = colBar;
		theData.colBard1{iFile} = colBard1;
		theData.colBard2{iFile} = colBard2;
		theData.colBard3{iFile} = colBard3;
		
		
		save(saveAs, '-mat', 'theData');

        %-- frame 8: blank again ------------------------------------------------------------------
        %Screen(window,'WaitBlanking');
        Screen('CopyWindow',offwin(iBlank),window);
        Screen('Flip', window, 0, 0, 2);

        %   trialOK = 1;
        if eyeTrackerON==1
            zerostr = '000';
            targfname = ['data', zerostr(1:3-length(num2str(iFile))), num2str(iFile),'.edf'];
            try
                transfer_EDF_mosaic(path,fname,targfname,1);
                translate_EDF(iFile, path, targfname);
            catch
                try
                    transfer_EDF_mosaic(path,fname,targfname,1);
                    translate_EDF(iFile, path, targfname);
                catch
                    try
                        transfer_EDF_mosaic(path,fname,targfname,1);
                        translate_EDF(iFile, path, targfname);
                    catch
                        clear all;
                        return;
                    end
                end
            end
        end


        %--flush variables and Screens - or memory explodes
        windowPtrs=Screen('Windows');
        for iPtr = 1:length(windowPtrs)
            isOffscreen=Screen(windowPtrs(iPtr),'IsOffscreen');
            if isOffscreen
                Screen(windowPtrs(iPtr),'Close')
            end
        end
        % avoid trial repeats
        totalFPx(selTrial) = NaN;
        fprintf(currentLogID,['Trial ', num2str(iTrial), ' took ' , num2str(GetSecs-theTimes.trialStart) 's and looks good!\r\n']);
        fprintf(currentLogID,[num2str(sum(~isnan(totalFPx))),' trials remain\r\n\r\n']);
        fprintf(sessionLogID,['Trial ', num2str(iTrial), ' took ' , num2str(GetSecs-theTimes.trialStart) 's and looks good!\r\n']);
        fprintf(sessionLogID,[num2str(sum(~isnan(totalFPx))),' trials remain\r\n\r\n']);

        iTrial = iTrial + 1;

        iFile = iFile + 1;

        goodTrialCount = goodTrialCount+1;
    end
end
if eyeTrackerON==1
    Eyelink('closefile');
    Eyelink('Shutdown');
end
Screen('closeall');
ShowCursor;
end



function cID = computerID()
cID = '';
ni = java.net.NetworkInterface.getNetworkInterfaces;
while ni.hasMoreElements
    addr = ni.nextElement.getHardwareAddress;
    if ~isempty(addr)
        addrStr = dec2hex(int16(addr)+128);
        cID = [cID, '.', reshape(addrStr,1,2*length(addr))];
    end
end

end


function ptcpDirInfo = GetPtcpDirInfo(dataDir, ptcpPrefix)
ptcpDirs = dir([dataDir,filesep,ptcpPrefix,'*']);
if length(ptcpDirs) >= 1
    for iPtcp = 1:length(ptcpDirs)
        thisPtcp = ptcpDirs(iPtcp).name;
        thisPtcpNo = str2num(thisPtcp(length(ptcpPrefix)+1:end));
        ptcpDirInfo.ptcpDirName{iPtcp} = thisPtcp;
        ptcpDirInfo.ptcpNo(iPtcp) = thisPtcpNo;
    end
else
    ptcpDirInfo.ptcpDirName = [];
    ptcpDirInfo.ptcpNo = [];
end
end

function [trialRejected, eye, notes] = fixationCheck(holdTime,bufferSize,fixationOffset,fixLocation,allowableOffset,eyeTrackerON,eye_used,el,pixPERdeg)
startTime = GetSecs;
buffered = false;
bufferPos = 1;
eye.x = 0;
eye.y = 0;
trialRejected = false;
notes = '';
while GetSecs < startTime+holdTime
    if eyeTrackerON==1 %
        error=Eyelink('CheckRecording');
        if(error~=0)
            break;
        end
        if Eyelink( 'NewFloatSampleAvailable') > 0
            % get the sample in the form of an event structure
            evt = Eyelink( 'NewestFloatSample');
            % if we do, get current gaze position from sample
            liveEye.x = evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
            liveEye.y = evt.gy(eye_used+1);
            % do we have valid data and is the pupil visible?
            if liveEye.x~=el.MISSING_DATA && liveEye.y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
                if buffered || bufferPos > bufferSize
                    eye.x = mean(eyeBuffer.x);
                    eye.y = mean(eyeBuffer.y);
                    %%%Check to see if the eye position is within
                    %%%range:
                    fixDiff.x = fixLocation.x - eye.x - fixationOffset.x;
                    fixDiff.y = fixLocation.y - eye.y - fixationOffset.y;
                    fixDiff.rho = sqrt((fixDiff.x^2)+(fixDiff.y^2));
                    notes = ['fixation offset was ' num2str(fixDiff.rho/pixPERdeg) ' degrees\r\n'];
                    if (fixDiff.rho/pixPERdeg) > allowableOffset
                        trialRejected = true;
                        break
                    end
                    buffered = true;
                    bufferPos = 1;
                end
                eyeBuffer.x(bufferPos)=liveEye.x;
                eyeBuffer.y(bufferPos)=liveEye.y;
                bufferPos = bufferPos+1;
            end
        end
    else
        % Query current mouse cursor position (our "pseudo-eyetracker") -
        [eye.x, eye.y, notUsed]=GetMouse;
    end
end
end

function [pauseCheck] = myKBInput_PausePlease(pauseKey)
pauseCheck = false;
[kbKeyDown, responseClock, keyCode] = KbCheck;
if kbKeyDown
    KBkeys = KbName(keyCode);
    if iscell(KBkeys)
        KBkey = KBkeys{1};
    else
        KBkey = KBkeys;
    end

    if strcmp(KBkey,pauseKey)
        pauseCheck = true;
    else
    end
end
end

function [keyResult] = pauseScreen(window, screenRect, grey, white)

myKeys.unPauseKey = 'u';
myKeys.recalibrate = 'r';
myKeys.quit = 'q';
Screen(window,'FillRect',grey);

Screen('TextSize',window, 36);

myPrompt = ['Paused!'];
Screen(window,'DrawText',myPrompt,50,50,white);

myPrompt = [myKeys.unPauseKey ' = unpause'];
Screen(window,'DrawText',myPrompt,50,150,white);

myPrompt = [myKeys.recalibrate ' = recalibrate'];
Screen(window,'DrawText',myPrompt,50,250,white);

myPrompt = [myKeys.quit ' = quit'];
Screen(window,'DrawText',myPrompt,50,350,white);

Screen('Flip', window, 0, 0, 2);

validResp = false;
while ~validResp
    [keyResult, validResp] = myKBInput_PausePrompt(myKeys);
end


if keyResult.quit
    quitProg
end

end

function [keyResult, validResp] = myKBInput_PausePrompt(myKeys)

validResp = false;
keyResult.unpause = false;
keyResult.recalibrate = false;
keyResult.quit = false;

pauseKey = myKeys.unPauseKey;
recalKey = myKeys.recalibrate;
quitKey = myKeys.quit;

[kbKeyDown, responseClock, keyCode] = KbCheck;
if kbKeyDown
    KBkeys = KbName(keyCode);
    if iscell(KBkeys)
        KBkey = KBkeys{1};
    else
        KBkey = KBkeys;
    end

    if strcmp(KBkey,pauseKey)
        keyResult.unpause = true;
        validResp = true;
    elseif strcmp(KBkey,recalKey)
        keyResult.recalibrate = true;
        validResp = true;
    elseif strcmp(KBkey,quitKey)
        keyResult.quit = true;
        validResp = true;
    end
end

end

function [el] = calibrateEL(el,window,recal)
if recal
    Eyelink('CloseFile')
    Eyelink('Stoprecording')
end
Eyelink('Shutdown');
if Eyelink('Initialize');
    return;
end;
EyelinkInit(0);
% STEP 3
% Provide eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and eyelink key  values).
el=EyelinkInitDefaults(window);
% make sure that we get gaze data from the eyelink
Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');
Eyelink('Command', 'link_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY');
Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,BLINK,SACCADE,BUTTON');
Eyelink('Command', 'heuristic_filter = 0')
Eyelink('Command', 'set_calibration_colors([100 0 0], [0 0 0])');

el.calibrationtargetsize = 1;
el.calibrationtargetwidth = .5;
el.displayCalResults = 1;
el.eyeimgsize = 50;

% STEP 4
% Calibrate the eye tracker using the standard calibration routines
%eyelink('trackersetup');
%     % you must send this command with value NO for custom calibration
%     % you must also reset it to YES for subsequent experiments


Eyelink('command', 'calibration_type = HV9');
% you must send this command with value NO for custom calibration
% you must also reset it to YES for subsequent experiments
Eyelink('command', 'generate_default_targets = NO');

% STEP 5.1 modify calibration and validation target locations
Eyelink('command','calibration_targets = 1360,540 960,140 560,940 560,540 560,140 960,540 960,940 1360,140 1360,940');
Eyelink('command','validation_targets = 1360,540 960,140 560,940 560,540 560,140 960,540 960,940 1360,140 1360,940');


%Eyelink('StartSetup');
errorID=EyelinkDoTrackerSetup(el)
% do a final check of calibration using driftcorrection
success=EyelinkDoDriftCorrection(el,[],[],1,1)
if isfield(el,'simarDebug')
    sca
    x = 1
end

if success ~= 1
    return;
end
end

function [] = quitProg()
Priority(0);
ShowCursor
Screen('CloseAll');
error('Experiment terminated!')
end

function [] = blockScreen(window, screenRect, grey, white, conditionChange)
Screen(window,'FillRect',grey);

% Below are prompts that show up for certain situations; very
% self-explanatory.
if ~conditionChange
myPrompt = ['Block complete!'];
Screen(window,'DrawText',myPrompt,50,50,white);

myPrompt = ['Please take a short break if you like'];
Screen(window,'DrawText',myPrompt,50,150,white);

myPrompt = ['If there''s anything you need, please call for the experimenter'];
Screen(window,'DrawText',myPrompt,50,250,white);

myPrompt = ['When you''re ready, click mouse to continue'];
Screen(window,'DrawText',myPrompt,50,350,white);

else
	
myPrompt = ['We''re going to make a change now'];
Screen(window,'DrawText',myPrompt,50,50,white);

myPrompt = ['Please call for the experimenter'];
Screen(window,'DrawText',myPrompt,50,150,white);	
end

Screen('Flip', window, 0, 0, 2);
GetClicks;

end