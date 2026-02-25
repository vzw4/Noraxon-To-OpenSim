% 3.15.2025
% Vicki Wang
% Human Engineering Research Laboratories

% This code converts Noraxon files into:
% 1 - marker data for the entire trial (.trc)
% 2 - GRF and CoP data for the trial (.mot)
% It will run through all the .csv files in a directory

% Data should be organized by separate subject subfolders in a parent
% directory, each subject subfolder should contain every trial CSV you'd
% like to analyze. Parent directory needs:
    % "NoraxonTrajectoriesAndForcesAll.mat" file
    % 

%%
clear
clc
close all

%% Functions

% Read data from the Noraxon .csv files
function trial_data = readNoraxon(trial_read, trajectoriesFixed, trajectoriesOriginal, forces)
    % Read time data
    trial_data.time = table2array(trial_read(:, 'time'));
    for i = 1:length(trajectoriesOriginal) % Read virtual marker data
        trial_data.(trajectoriesFixed{i}) = table2array(trial_read(:, trajectoriesOriginal{i}));
    end
        
    for i = 1:length(forces) % Read insole data
        trial_data.(forces{i}) = table2array(trial_read(:, forces(i)));
    end
end

% Noraxon insoles offset on the x axis by 1.8" - correct for this
% Also, change NaNs to zeros, and make corresponding GRF also zero
function [fixCoPx, fixCoPy, fixGRF] = fixCoPNaNs(rawCoPx, rawCoPy, rawGRF)
    fixCoPx = rawCoPx;
    fixCoPy = rawCoPy;
    fixGRF = rawGRF;
    for t = 1:length(rawCoPx)
        if isnan(rawCoPx(t))
            fixCoPx(t) = 0;
            fixCoPy(t) = 0;
            fixGRF(t) = 0;
        else
            fixCoPx(t) = rawCoPx(t) - 50;
            fixCoPy(t) = rawCoPy(t) - 50;
        end
    end
end

% Switch y and z, then negate the new z
function rotated_data = rotateAroundX(trial_data, trajectories)
        rotated_data = trial_data;
        for i = 1:length(trajectories)
            if contains(trajectories{i}, '_y')
                markerName = erase(trajectories{i}, '_y');
                YZmarkers = {};
                YZmarkers{1} = append(markerName, '_y');
                YZmarkers{2} = append(markerName, '_z');
                tempY = trial_data.(YZmarkers{1});
                tempZ = trial_data.(YZmarkers{2});
                rotated_data.(YZmarkers{1}) = tempZ;
                rotated_data.(YZmarkers{2}) = -tempY;
            end
        end
    end

%% Import data

% Read directory of current subject
parentDirectory = uigetdir('', 'Select the parent folder that contains all subject subfolders'); % Select parent directory
cd(parentDirectory);
a = dir(parentDirectory); b = [a.isdir]; c = a(b); d = {c.name}; 
subjectSubfoldersList = d(:, 3:end); % Get list of all subject subfolders in the parent directory
clear("a", "b", "c", "d");

%%

% Load template info
load("NoraxonTrajectoriesAndForcesAll.mat") % Loads trajectory and force names with original and fixed names
templateTRCname = fullfile(parentDirectory, "NoraxonMarkerBlank.txt");
templateMOTname = fullfile(parentDirectory, "grf_motBLANK.txt");

for i = 1:length(subjectSubfoldersList)
    subjectFolderAddress = fullfile(parentDirectory, subjectSubfoldersList{i}); % Set address of current subject subfolder
    cd(subjectFolderAddress); % Set current directory to this subject subfolder
    load('subjectData.mat'); % Load the subject data (name and mass)
    
    % Get list of all .trc trials
    a = dir('*.csv');
    subjectCSVTrials = {a.name}; clear("a");

        for j = 1:length(subjectCSVTrials) % Loop through the trials

    
            % Change these as needed to pull files from your system
            csvName = subjectCSVTrials{j};
            trialID = csvName(1:end-4);
            data_read = readtable(subjectCSVTrials{j});
        
            currentMarkerSet = {};
        
            trial_trc = readcell(templateTRCname);
            trial_mot = readcell(templateMOTname);
            
            trial_data = readNoraxon(data_read, trajectoriesFixed, trajectoriesOriginal, forces);    
            rotated_data = rotateAroundX(trial_data, trajectoriesFixed);
            
            %% ROTATE AND TRANSLATE CENTER OF PRESSURE
            
            % Get rid of NaNs in CoP data
            [rotated_data.noNaN_CoP_LT_x, rotated_data.noNaN_CoP_LT_y, rotated_data.noNaN_GRF_LT] = fixCoPNaNs(rotated_data.Insole_LTInsole_CenterOfPressure_x_mm_, rotated_data.Insole_LTInsole_CenterOfPressure_y_mm_, rotated_data.LTInsole_Total___);
            [rotated_data.noNaN_CoP_RT_x, rotated_data.noNaN_CoP_RT_y, rotated_data.noNaN_GRF_RT] = fixCoPNaNs(rotated_data.Insole_RTInsole_CenterOfPressure_x_mm_, rotated_data.Insole_RTInsole_CenterOfPressure_y_mm_, rotated_data.RTInsole_Total___);
            
            % Midpoint - x
            footVector.frontLT_x = (rotated_data.MedMetatarsophalangealJointLT_x + rotated_data.LatMetatarsophalangealJointLT_x)/2;
            footVector.frontLT_y = (rotated_data.MedMetatarsophalangealJointLT_y + rotated_data.LatMetatarsophalangealJointLT_y)/2;
            footVector.frontLT_z = (rotated_data.MedMetatarsophalangealJointLT_z + rotated_data.LatMetatarsophalangealJointLT_z)/2;
            footVector.frontRT_x = (rotated_data.MedMetatarsophalangealJointRT_x + rotated_data.LatMetatarsophalangealJointRT_x)/2;
            footVector.frontRT_y = (rotated_data.MedMetatarsophalangealJointRT_y + rotated_data.LatMetatarsophalangealJointRT_y)/2;
            footVector.frontRT_z = (rotated_data.MedMetatarsophalangealJointRT_z + rotated_data.LatMetatarsophalangealJointRT_z)/2;
            
            % Center trajectories with the heel as the origin
            footVector.frontLT_x_centered = footVector.frontLT_x - rotated_data.HeelBackLT_x;
            footVector.frontLT_y_centered = footVector.frontLT_y - rotated_data.HeelBackLT_y;
            footVector.frontLT_z_centered = footVector.frontLT_z - rotated_data.HeelBackLT_z;
            footVector.frontRT_x_centered = footVector.frontRT_x - rotated_data.HeelBackRT_x;
            footVector.frontRT_y_centered = footVector.frontRT_y - rotated_data.HeelBackRT_y;
            footVector.frontRT_z_centered = footVector.frontRT_z - rotated_data.HeelBackRT_z;
            
            % Find rotation matrices through each time point
            for i = 1:length(rotated_data.time)
                footVecLT = [footVector.frontLT_x_centered(i); footVector.frontLT_y_centered(i); footVector.frontLT_z_centered(i)];
                footVecRT = [footVector.frontRT_x_centered(i); footVector.frontRT_y_centered(i); footVector.frontRT_z_centered(i)];
        
                CopAdjust = 30;
        
                CoPVecLT = [rotated_data.noNaN_CoP_LT_y(i)+CopAdjust; 0; rotated_data.noNaN_CoP_LT_x(i)];
                CoPVecRT = [rotated_data.noNaN_CoP_RT_y(i)+CopAdjust; 0; rotated_data.noNaN_CoP_RT_x(i)];
                normFootLT = norm(footVecLT);
                normFootRT = norm(footVecRT);
                footUnitVecLT = footVecLT/normFootLT;
                footUnitVecRT = footVecRT/normFootRT;
                insoleUnitVec = [1; 0; 0];
        
                % LEFT
                % p0 = footUnitVecLT;
                % p1 = insoleUnitVec;
                p0 = insoleUnitVec;
                p1 = footUnitVecLT;
                % calculate cross and dot products
                C = cross(p0, p1) ; 
                D = dot(p0, p1) ;
                NP0 = norm(p0) ; % used for scaling
                if ~all(C==0) % check for colinearity    
                    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
                    rotMatrixLT = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
                else
                    rotMatrixLT = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
                end
            
                % RIGHT
                % p0 = footUnitVecRT;
                % p1 = insoleUnitVec;
                p0 = insoleUnitVec;
                p1 = footUnitVecRT;
                % calculate cross and dot products
                C = cross(p0, p1) ; 
                D = dot(p0, p1) ;
                NP0 = norm(p0) ; % used for scaling
                if ~all(C==0) % check for colinearity    
                    Z = [0 -C(3) C(2); C(3) 0 -C(1); -C(2) C(1) 0] ; 
                    rotMatrixRT = (eye(3) + Z + Z^2 * (1-D)/(norm(C)^2)) / NP0^2 ; % rotation matrix
                else
                    rotMatrixRT = sign(D) * (norm(p1) / NP0) ; % orientation and scaling
                end 
            
                rotatedCoPLT = rotMatrixLT*CoPVecLT;
                rotatedCoPRT = rotMatrixRT*CoPVecRT;
        
                rotated_data.CoP_LT_x_centered(i) = rotatedCoPLT(1,1);
                rotated_data.CoP_LT_y_centered(i) = rotatedCoPLT(2,1);
                rotated_data.CoP_LT_z_centered(i) = rotatedCoPLT(3,1);
                rotated_data.CoP_RT_x_centered(i) = rotatedCoPRT(1,1);
                rotated_data.CoP_RT_y_centered(i) = rotatedCoPRT(2,1);
                rotated_data.CoP_RT_z_centered(i) = rotatedCoPRT(3,1);
            end
           
            % Translate CoPs to move with the foot through space
            for t = 1:length(rotated_data.time)
                rotated_data.CoP_LT_x(t) = rotated_data.CoP_LT_x_centered(t) + rotated_data.HeelBackLT_x(t);
                rotated_data.CoP_LT_y(t) = rotated_data.CoP_LT_y_centered(t) + rotated_data.HeelBackLT_y(t);
                rotated_data.CoP_LT_z(t) = rotated_data.CoP_LT_z_centered(t) + rotated_data.HeelBackLT_z(t);
                rotated_data.CoP_RT_x(t) = rotated_data.CoP_RT_x_centered(t) + rotated_data.HeelBackRT_x(t);
                rotated_data.CoP_RT_y(t) = rotated_data.CoP_RT_y_centered(t) + rotated_data.HeelBackRT_y(t);
                rotated_data.CoP_RT_z(t) = rotated_data.CoP_RT_z_centered(t) + rotated_data.HeelBackRT_z(t);
            end
        
            % Convert from mm to m
            rotated_data.CoP_LT_x = rotated_data.CoP_LT_x / 1000;
            rotated_data.CoP_LT_y = rotated_data.CoP_LT_y / 1000;
            rotated_data.CoP_LT_z = rotated_data.CoP_LT_z / 1000; 
            rotated_data.CoP_RT_x = rotated_data.CoP_RT_x / 1000;
            rotated_data.CoP_RT_y = rotated_data.CoP_RT_y / 1000;
            rotated_data.CoP_RT_z = rotated_data.CoP_RT_z / 1000; 
            
            %% CONVERT GRFS TO NEWTONS
            rotated_data.GRF_LT = (rotated_data.noNaN_GRF_LT ./ 100) .*subjectMass.*4.448;
            rotated_data.GRF_RT = (rotated_data.noNaN_GRF_RT ./ 100) .*subjectMass.*4.448;
            
            %% CALCULATE GROUND REACTION MOMENTS
            rotated_data.M_LT_x = rotated_data.GRF_LT.*rotated_data.CoP_LT_z;
            % rotated_data.M_LT_y = 0;
            rotated_data.M_LT_z = rotated_data.GRF_LT.*rotated_data.CoP_LT_x;
            rotated_data.M_RT_x = rotated_data.GRF_RT.*rotated_data.CoP_RT_z;
            % rotated_data.M_RT_y = 0;
            rotated_data.M_RT_z = rotated_data.GRF_RT.*rotated_data.CoP_RT_x;
        
            
            %% SET UP .TRC FILES
            
            % Set up header info
            trial_trcName = append(subjectName, '_', trialID, '.trc');
            dataCount = length(rotated_data.time);
            
            % Add specific header info
            trial_trc{2,4} = trial_trcName;
            trial_trc{4,3} = dataCount;
            trial_trc{4,8} = dataCount;
            
            % Add frame, time, and marker coordinate data
            for i = 1:dataCount
                trial_trc{i+6,1} = i; % Frame#
                trial_trc{i+6,2} = i*0.002-0.002; % Time
                for k = 1:length(trajectoriesFixed) % Marker coordinate data
                    trial_trc{i+6,k+2} = rotated_data.(trajectoriesFixed{k})(i);
                end
            end
            
            % Remove <missing> in cells
            mask = cellfun(@(x) any(isa(x,'missing')), trial_trc);
            trial_trc(mask) = {[]}; % or whatever value you want to use
            
        
            %% WRITE MARKER TRC FILES
        
            % Trial
            trial_fileName = trial_trcName(1:end-4);
            fullfilePath = fullfile(trial_fileName);
            writecell(trial_trc(2:end,1:end), fullfilePath, 'Delimiter', '\t');
        
            % Change files from .txt to .trc
            file1 = append(fullfilePath, '.txt');
            file2 = strrep(file1,'.txt','.trc');
            copyfile(file1,file2);
            delete(file1);
        
        
            %% SET UP .MOT FILE
        
            % Set up header info
            trial_motName = append(subjectName, '_', trialID, '_GRF', '.mot');
            dataCount = length(trial_data.time);
        
            % Add specific header info
            trial_mot{1,1} = trial_motName;
            trial_mot{3,1} = append('nRows=',num2str(dataCount));
        
            % Add frame, time, and marker coordinate data
            for i = 1:dataCount
                trial_mot{i+7,1} = i*0.002-0.002; % Time
                trial_mot{i+7,2} = 0; % Left GRF x
                trial_mot{i+7,3} = rotated_data.GRF_LT(i); % Left GRF y
                trial_mot{i+7,4} = 0; % Left GRF z
                trial_mot{i+7,5} = rotated_data.CoP_LT_x(i); % Left CoP x
                trial_mot{i+7,6} = rotated_data.CoP_LT_y(i); % Left CoP y
                trial_mot{i+7,7} = rotated_data.CoP_LT_z(i); % Left CoP z
                %
                trial_mot{i+7,8} = 0; % Right GRF x
                trial_mot{i+7,9} = rotated_data.GRF_RT(i); % Right GRF y
                trial_mot{i+7,10} = 0; % Right GRF z
                trial_mot{i+7,11} = rotated_data.CoP_RT_x(i); % Right CoP x
                trial_mot{i+7,12} = rotated_data.CoP_RT_y(i); % Right CoP y
                trial_mot{i+7,13} = rotated_data.CoP_RT_z(i); % Right CoP z
                %
                trial_mot{i+7,14} = rotated_data.M_LT_x(i); % Left Moment x
                trial_mot{i+7,15} = 0; % Left Moment y
                trial_mot{i+7,16} = rotated_data.M_LT_z(i); % Left Moment z
                %
                trial_mot{i+7,17} = rotated_data.M_RT_x(i); % Right Moment x
                trial_mot{i+7,18} = 0; % Right Moment y
                trial_mot{i+7,19} = rotated_data.M_RT_z(i); % Right Moment z
            end
        
            % Remove <missing> in cells
            mask = cellfun(@(x) any(isa(x,'missing')), trial_mot);
            trial_mot(mask) = {[]};
        
            %% WRITE GRF MOT FILE
        
            trial_fileName = trial_motName(1:end-4);
            fullfilePath = fullfile(trial_fileName);
            writecell(trial_mot, fullfilePath, 'Delimiter', '\t');
        
            % Change files from .txt to .mot
            file1 = append(fullfilePath, '.txt');
            file2 = strrep(file1,'.txt','.mot');
            copyfile(file1,file2);
            delete(file1);
        end
end