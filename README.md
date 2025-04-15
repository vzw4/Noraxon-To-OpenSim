# Noraxon-To-OpenSim
Converts Noraxon IMU/Insole data into .trc and .mot file for OpenSim biomechanical analysis

This code utilizes trajectory outputs from Noraxon IMUs as virtual bony markers to create a .trc file of coordinates for each marker. It also converts Noraxon insole data into a vGRFs and CoP .mot to use in OpenSim.

Notes: This code assumes that each trial begins with the subject standing in a neutral pose. It uses the first 0.2s of each trial to create a static model for scaling in OpenSim.

Steps:
1. Ensure you have the following files:
   a. grf_motBLANK.txt
   b. marker_trcBLANK.txt
   c. NoraxonToOpenSim.m
   d. NoraxonTrajectoriesAndForces.mat
   e. Data from whatever you measured with the Noraxon IMUs, in .csv form. Example data is included
   f. NoraxonMarkersCLOSE.xml
   g. Scale_Setup.xml
3. Edit lines 26-29 to read files from your system 
4. Change the parameters in lines 17-19 (Subject name, bodyweight in lbs, Noraxon data file name)
5. Run NoraxonToOpenSim. Output should be two .trc files (one for the trial and one in a static pose) and one .mot file (contains insole data)
6. Edit line 50 in grf_template.xml to read the GRF .mot file you just created
7. In OpenSim now - load gait2354_simbody.osim model
8. Scale model, load Scale_Setup.xml and change fields on the 'Settings' tab
9. Run Inverse Kinematics, load the trial .trc
10. Run Inverse Dynamics, load the GRF .xml
