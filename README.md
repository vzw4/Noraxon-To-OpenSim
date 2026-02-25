# Noraxon-To-OpenSim
Batch converts Noraxon IMU/Insole data into .trc and .mot file for OpenSim biomechanical analysis

This code utilizes trajectory outputs from Noraxon IMUs as virtual bony markers to create a .trc file of coordinates for each marker. It also converts Noraxon insole data into a vGRFs and CoP .mot to use in OpenSim. It then uses a MATLAB wrapper to run kinematics, kinetics, JointReactions, and BodyKinematics analysis in OpenSim. 

Notes: This code assumes that each trial begins with the subject standing in a neutral pose. It uses the first 0.2s of each trial to create a static model for scaling in OpenSim.

For this code to work, you need to organize your data like this:
1. Parent Directory
   
   a. grf_motBLANK.txt (blank template file to load GRF data into)

   b. NoraxonMarkerBlank.txt (blank template file to load virtual marker data into)
   
   c. NoraxonTrajectoriesAndForces.mat (This helps rename the columns from the Noraxon data to be more manageable)
     
   d. LFB_Noraxon_Markers.xml (Marker data for this model and from this equipment)
   
   e. NoraxonScaleTemplate.xml
   
   f. KinematicsTemplate.xml
   
   g. KineticsTemplate.xml

   h. GRFTemplate.xml

   i. AnalyzeTemplate.xml

   j. LFB_model.osim (the base OpenSim model)

   k. Subject subfolders (a separate folder for each subject)

      1. all trial CSV files exported by Noraxon ("Export as Single CSVs")
  
      2. subjectData.mat (create this manually in MATLAB - contains subjectMass and subjectName. Mass should be in kg.)

Then, follow the steps below to generate the OpenSim outputs for each trial.

3. Run MATLAB_fullPipeline_CSVtoTRCMOT.m

4. Run MATLAB_fullPipeline.m
