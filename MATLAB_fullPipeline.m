% 2.13.2026
% Batch run model scaling, inverse kinematics, inverse dynamics, joint
% reactions, and body kinematics (analysis tool)
% Note: this code assumes that you have generated proper .trc and .mot
% files, from the NORAXON SYSTEM ONLY. May need to first use the
% MATLAB_fullPipeline_CSVtoTRCMOT.m code to convert Noraxon exported CSVs to usable TRC and MOT files.

% Data Structure:
    % Parent folder
        % One copy of the base, unscaled model for analysis
        % Template Scaling, IK, ID, Analysis XMLs
        % Subject subfolders
            % All subject trial .trc's
            % All subject trials GRF .mots
            % subjectData.mat <-- has subjectMass, subjectName

% Code will add to each subfolder:
    % Kinematics .mot file
    % Kinetics .sto file
    % JointReactions .sto file
    % BodyKinematics .sto files (3: pos, vel, acc)

% 1. Load in template files
% 2. Select parent directory to work out of (contains subject subfolders)
% Within each subject subfolder:
    % 3. Create scaled model
    % For each subject's trial:
        % 4. Set IK parameters and run IK
        % 5. Set ID parameters and run ID
        % 6. Set Analysis parameters and run for JointReactions and BodyKinematics

%% SETUP
clear
clc
close all

import org.opensim.modeling.* % Set up proper OpenSim stuff

%% SELECT WORKING DIRECTORY
parentDirectory = uigetdir('', 'Select the parent folder that contains all subject subfolders'); % Select parent directory
cd(parentDirectory);
a = dir(parentDirectory); b = [a.isdir]; c = a(b); d = {c.name}; 
subjectSubfoldersList = d(:, 3:end); % Get list of all subject subfolders in the parent directory
clear("a", "b", "c", "d");

%% LOAD TEMPLATE FILES

kinematicsTool = InverseKinematicsTool("KinematicsTemplate.xml");
kineticsTool = InverseDynamicsTool("KineticsTemplate.xml");
analyzeTool = AnalyzeTool("AnalyzeTemplate.xml");
baseModel = Model("LFB_model.osim");
baseModel.initSystem();

%% SCALING PERFORMS ONCE FOR EACH SUBJECT
% For each trial in each subject subfolder in the parent directory
for i = 1:length(subjectSubfoldersList)
    subjectFolderAddress = fullfile(parentDirectory, subjectSubfoldersList{i}); % Set address of current subject subfolder
    cd(subjectFolderAddress); % Set current directory to this subject subfolder
    load('subjectData.mat'); % Load the subject data (name and mass)
    
    % Get list of all .trc trials
    a = dir('*.trc');
    subjectTRCTrials = {a.name}; clear("a");

    %% SCALING
    % Couldn't get the MATLAB API calls to work properly, so this edits the
    % XML directly, saves it as a new one, and then executes it to scale
    % the model for each subject

    xmlPath = fullfile(parentDirectory, "NoraxonScaleTemplate.xml"); % Read in template XML
    xmlFile = xmlread(xmlPath);
    
    lfbAddress = fullfile(parentDirectory, "LFB_model.osim"); % Set the base model
    xmlFile.getElementsByTagName('model_file').item(0).setTextContent(lfbAddress);

    noraxonMarkersAddress = fullfile(parentDirectory, "LFB_Noraxon_Markers.xml"); % Assign marker set
    xmlFile.getElementsByTagName('marker_set_file').item(0).setTextContent(noraxonMarkersAddress);

    trialForScale = subjectTRCTrials{1}; % Use whatever first trial as pseudostatic to scale
    % Note: for CATT participants, may want to instead use the walking
    % calibration as pseudostatic
    trialScaleAddress = fullfile(subjectFolderAddress, trialForScale);
    xmlFile.getElementsByTagName('marker_file').item(0).setTextContent(trialScaleAddress);
    xmlFile.getElementsByTagName('marker_file').item(1).setTextContent(trialScaleAddress);

    xmlFile.getElementsByTagName('mass').item(0).setTextContent(num2str(subjectMass)); % Set subject mass

    outputModelName = append(subjectName, '_scaledModel.osim'); % Set the scaled model output
    xmlFile.getElementsByTagName('output_model_file').item(0).setTextContent(outputModelName);

    scaleXMLName = [subjectName '_scaleTool.xml'];
    xmlwrite(scaleXMLName, xmlFile); % Save the edited XML
    
    scaleTool = ScaleTool(scaleXMLName);
    scaleTool.run(); % Execute the new scaling XML
    
    scaledModel = Model(outputModelName);
    scaledModel.initSystem(); % Load and initialize the new scaled model

    for j = 1:length(subjectTRCTrials) % Loop through the trials

        %% IK
    
        kinematicsTool.setModel(scaledModel);

        % Get the name of the file for this trial
        markerFile = subjectTRCTrials{j};

        % Create name of trial from .trc file name
        trialName = regexprep(markerFile,'.trc','');

        % Get trc data to determine time range
        markerData = MarkerData(markerFile);

        % Get initial and intial time 
        initial_time = markerData.getStartFrameTime();
        final_time = markerData.getLastFrameTime();

        % Setup the kinematicsTool for this trial
        kinematicsTool.setName(trialName);
        kinematicsTool.setMarkerDataFileName(markerFile);
        kinematicsTool.setStartTime(initial_time);
        kinematicsTool.setEndTime(final_time);
        outputIK = [trialName '_ik.mot'];
        kinematicsTool.setOutputMotionFileName(outputIK);

        fprintf(['Performing IK on ' trialName '\n']);
        % Run IK
        kinematicsTool.run();


        %% ID
        kineticsTool.setModel(scaledModel);

        % Edit ExternalLoads xml
        xmlPath = fullfile(parentDirectory, "GRFTemplate.xml"); % Read in template XML
        xmlFile = xmlread(xmlPath);
        grfName = [trialName '_GRF.mot'];
        xmlFile.getElementsByTagName('datafile').item(0).setTextContent(grfName);
        grfXMLName = [trialName '_GRF.xml'];
        xmlwrite(grfXMLName, xmlFile); % Save the edited XML

        kineticsTool.setExternalLoadsFileName(grfXMLName);
        kineticsTool.setCoordinatesFileName(outputIK);

        % 4. Set Time Range and Output
        kineticsTool.setStartTime(initial_time);
        kineticsTool.setEndTime(final_time);
        outputID = [trialName '_ID.sto'];
        kineticsTool.setResultsDir(cd);
        outputIDpath = fullfile(cd, outputID);
        kineticsTool.setOutputGenForceFileName(outputIDpath);

        kineticsTool.setLowpassCutoffFrequency(6.0); % Filter data
        allForces = scaledModel.getForceSet();
        excluded = ArrayStr();
        for p = 0:allForces.getSize()-1
            excluded.append(allForces.get(p).getName());
        end
        kineticsTool.setExcludedForces(excluded);        
        % 5. Run the Tool
        fprintf(['Performing ID on ' trialName '\n']);
        kineticsTool.run();

        %% ANALYZE - JOINTREACTIONS AND BODYKINEMATICS

        % Change subjectID
        analyzeTool.setName(trialName); 

        % Change .mot kinematics file to use
        analyzeTool.setCoordinatesFileName(outputIK);

        % Change .xml GRF file to use=
        analyzeTool.setExternalLoadsFileName(grfXMLName);

        % Change model
        analyzeTool.setModelFilename(outputModelName);

        fprintf(['Performing Analysis on ' trialName '\n']);
        % Run!
        analyzeTool.run();
    end
end