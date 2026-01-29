%%
clear
clc
close all

%%
import org.opensim.modeling.* % Set up proper OpenSim stuff

currentFolder = pwd();
AnalyzeSetup = fullfile(currentFolder, "AnalyzeSetup_Template.xml"); % Assign the template XML that sets up the Analyze tool
analyzeTool = AnalyzeTool(AnalyzeSetup);

% Pick the scaled .osim model
[modelName, modelAddress] = uigetfile('*.osim', 'Select the .osim model');
% Pick .mot file (kinematics)
[motName, motAddress] = uigetfile('*.mot', 'Select the kinematics .mot file');
% Pick GRF.xml file
[xmlName, xmlAddress] = uigetfile('*.xml', 'Select the .xml GRF file');


%%

grfFile = fullfile(xmlAddress, xmlName);
kinFile = fullfile(motAddress, motName);
modelFile = fullfile(modelAddress, modelName);

% Change subjectID
trialID = input('Enter trial ID here: ', "s");
analyzeTool.setName(trialID); 

% Change .mot kinematics file to use
analyzeTool.setCoordinatesFileName(kinFile);

% Change .xml GRF file to use
analyzeTool.setExternalLoadsFileName(grfFile);

% Change model
model = fullfile(modelFile);
analyzeTool.setModelFilename(modelName);

% Change results directory
analyzeTool.setResultsDir(modelAddress);

% Run!
analyzeTool.run();


%% EXTERNAL MOMENTS



%% READING IN DATA
% Read in segment center of mass data from BodyKinematics output from OpenSim
% Pick .sto file that contains "pos"
    pickCOMFile = dir(fullfile(currentFolder, '*pos*.sto')).name;
    comFile = readtable(pickCOMFile, FileType="text"); 

% Read in joint center data from JointReactions output from OpenSim
    pickJRFile = dir(fullfile(currentFolder, '*ReactionLoads*.sto')).name;
    JRFile = readtable(pickJRFile, FileType="text");

% Read in GRF and CoP data from .mot file
    pickGRFFile = dir(fullfile(currentFolder, '*GRF*.mot')).name;
    grfData = readtable(pickGRFFile, FileType="text");
    grfColumns = ["time", ...
        "ground_force_vx", "ground_force_vy", "ground_force_vz", ...
        "ground_force_px", "ground_force_py", "ground_force_pz", ...
        "L_ground_force_vx", "L_ground_force_vy", "L_ground_force_vz", ...
        "L_ground_force_px", "L_ground_force_py", "L_ground_force_pz", ...
        "ground_torque_x", "ground_torque_y", "ground_torque_z", ...
        "L_ground_torque_x", "L_ground_torque_y", "L_ground_torque_z"]; % Assign the table variable names
    grfData.Properties.VariableNames = grfColumns; % Rename table columns
    grfData = grfData(7:end, :); % Remove first six rows of NaN values

    F_R = [grfData.ground_force_vx, ...
        grfData.ground_force_vy, ...
        grfData.ground_force_vz];

    COP_R = [grfData.ground_force_px, ...
        grfData.ground_force_py, ...
        grfData.ground_force_pz];

    F_L = [grfData.L_ground_force_vx, ...
        grfData.L_ground_force_vy, ...
        grfData.L_ground_force_vz];

    COP_L = [grfData.L_ground_force_px, ...     
        grfData.L_ground_force_py, ...
        grfData.L_ground_force_pz];

% Read in segment mass info from scaled OpenSim model
    osimFile = modelName;
    xDoc = xmlread(osimFile); % Read XML

    bodyList = xDoc.getElementsByTagName('Body'); % Get all Body elements
    nBodies = bodyList.getLength;
    
    % Preallocate
    bodyNames  = strings(nBodies,1);
    bodyMasses = zeros(nBodies,1);
    
    for i = 0:nBodies-1
        bodyNode = bodyList.item(i);
        bodyNames(i+1) = string(bodyNode.getAttribute('name')); % Body name
        massNodes = bodyNode.getElementsByTagName('mass'); % Find mass
    
        if massNodes.getLength > 0
            massValue = massNodes.item(0).getFirstChild.getData;
            bodyMasses(i+1) = str2double(massValue);
        else
            bodyMasses(i+1) = NaN;  % safety fallback
        end
    end
    
    bodyMassTable = table(bodyNames, bodyMasses, 'VariableNames', {'Body','Mass_kg'}); % Put into a table
   
%% COMPUTE A GENERAL UPPER BODY CENTER OF MASS
% Set which bodies are considered 'upper body' and will be included
    upperBodies = {'lumbar4', 'lumbar3', 'lumbar2', 'lumbar1', 'torso'}; % Note that arms are being excluded

% Compute a conglomerated upper body center of mass
    nFrames = height(comFile);

    rTrunk = zeros(nFrames, 3); % Preallocate
    totalMass = 0; % Initialize

    for i = 1:length(upperBodies)
        segment = upperBodies{i}; % Identify the trunk segment
        idx = ismember(bodyMassTable{:, 1}, segment); % Find the index of that segment in the body mass table
        mass = bodyMassTable{idx, 2}; % Find the mass of that segment

        rx = comFile.([segment '_X']); % Find the COM coordinates
        ry = comFile.([segment '_Y']);
        rz = comFile.([segment '_Z']);

        rTrunk = rTrunk + (mass * [rx ry rz]); % Multiply each mass by its position, and add to sum
        totalMass = totalMass + mass; % Add mass to total trunk mass
    end

    trunkCOM = rTrunk/totalMass; % Divide the summed mass*position values by the total trunk mass to obtain the weighted average COM for the trunk

% Create array for L5 COM
L5S1COM = [JRFile.('L5_S1_IVDjnt_on_lumbar5_in_ground_px'), JRFile.('L5_S1_IVDjnt_on_lumbar5_in_ground_py'), JRFile.('L5_S1_IVDjnt_on_lumbar5_in_ground_pz')];

% L5 orientation for vector projection correction
L5Orientation = deg2rad([comFile.('lumbar5_Ox'), comFile.('lumbar5_Oy'), comFile.('lumbar5_Oz')]);

%% ROTATE COP, COM SUCH THAT L5 MATCHES ASSUMED ORIENTATION FOR MOMENT CALCULATION
% +X is to the front, +Y is up, +Z is to the right
% X = Roll, Y = Yaw, Z = Pitch

% L5Orientation is the 'local frame'. All OpenSim data is given in 

%% Inputs
% L5Euler: 3xN matrix of [Roll; Yaw; Pitch] over time (radians)
% COP: 3xN matrix of GRF COP positions in global frame
% GRF: 3xN matrix of GRFs in global frame
% COM: 3xN matrix of segment COM positions in global frame (optional)
% mass: scalar or 1xN vector of segment masses
% g: gravity vector in global frame, e.g., [0; -9.81; 0]
% L5S1_joint_center assumed at origin [0;0;0]

L5Euler = L5Orientation';

numFrames = size(L5Euler,2);
M_ext_L5_total = zeros(3,numFrames);

%% Construct rotation matrices for all frames
% Preallocate
R_global_to_L5 = zeros(3,3,numFrames);


for i = 1:numFrames
    roll  = L5Euler(1,i);
    yaw   = L5Euler(2,i);
    pitch = L5Euler(3,i);

    % Rotation matrices around each axis
    Rx = [1 0 0;
          0 cos(roll) -sin(roll);
          0 sin(roll)  cos(roll)];

    Ry = [cos(yaw) 0 sin(yaw);
          0 1 0;
         -sin(yaw) 0 cos(yaw)];

    Rz = [cos(pitch) -sin(pitch) 0;
          sin(pitch)  cos(pitch) 0;
          0 0 1];

    % Combined rotation: L5 local -> global
    R_L5 = Rx * Ry * Rz;

    % Global -> L5 frame
    R_global_to_L5(:,:,i) = R_L5.';
end

%% Transform COP and GRF to L5 frame and compute moments
% Vectorized cross product using column-wise operation
L_COP_L5 = zeros(3,numFrames);
L_GRF_L5   = zeros(3,numFrames);
R_COP_L5 = zeros(3,numFrames);
R_GRF_L5   = zeros(3,numFrames);
Rot_COM_L5   = zeros(3,numFrames);
M_ext_L5_leftGRF   = zeros(3,numFrames);
M_ext_L5_rightGRF   = zeros(3,numFrames);
M_ext_trunkGravity   = zeros(3,numFrames);

LCOP = COP_L';
RCOP = COP_R';
COM = trunkCOM';
LGRF = F_L';
RGRF = F_R';

for i = 1:numFrames
    L_COP_L5(:,i) = R_global_to_L5(:,:,i) * LCOP(:,i);
    L_GRF_L5(:,i)   = R_global_to_L5(:,:,i) * LGRF(:,i);

    R_COP_L5(:,i) = R_global_to_L5(:,:,i) * RCOP(:,i);
    R_GRF_L5(:,i)   = R_global_to_L5(:,:,i) * RGRF(:,i);

    Rot_COM_L5(:,i) = R_global_to_L5(:,:,i) * COM(:,i);    

        forceGrav = [0 -totalMass*9.81 0];
        gravVec = repmat(forceGrav, height(trunkCOM), 1)';
    Rot_grav(:,i) = R_global_to_L5(:,:,i) * gravVec(:,i);    

    % External moment
    M_ext_L5_leftGRF(:,i) = cross(L_COP_L5(:,i), L_GRF_L5(:,i));
    M_ext_L5_rightGRF(:,i) = cross(R_COP_L5(:,i), R_GRF_L5(:,i));
    M_ext_trunkGravity(:, i) = cross(Rot_COM_L5(:,i), Rot_grav(:,i)); % moment from gravity

    M_ext_L5_total(:,i) = M_ext_L5_leftGRF(:,i) + M_ext_L5_rightGRF(:,i) + M_ext_trunkGravity(:, i);
end

%% Output
% M_ext_L5_total is 3xN matrix: [Mx; My; Mz] in L5 local frame over time

exMomentsFileName = append(trialID, "_externalMoments");

L5S1_External_Moments = M_ext_L5_total';

% Normalize to subject by mass
% mass = subjectMass * 0.453592;
% L5S1_External_Moments = L5S1_External_Moments/mass;

% 4th order low pass Butterworth filter
fc=10;
fs=500;
[g,h]=butter(4,fc/(fs/2),'low');
% freqz(g,h,[],fs);

L5S1_External_Moments = filtfilt(g, h, L5S1_External_Moments);

save(exMomentsFileName, 'L5S1_External_Moments');