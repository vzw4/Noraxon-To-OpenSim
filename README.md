# Noraxon-To-OpenSim
Converts Noraxon IMU/Insole data into .trc and .mot file for OpenSim biomechanical analysis

This code utilizes trajectory outputs from Noraxon IMUs as virtual bony markers to create a .trc file of coordinates for each marker. It also converts Noraxon insole data into a vGRFs and CoP .mot to use in OpenSim.

Notes: This code assumes that each trial begins with the subject standing in a neutral pose. It uses the first 0.2s of each trial to create a static model for scaling in OpenSim.

Steps:
1. Ensure you have the following files, all in one directory:
   
   a. grf_motBLANK.txt (blank template file to load GRF data into)

   b. marker_trcBLANK.txt (blank template file to load virtual marker data into)
   
   c. NoraxonToOpenSimSINGLE.m (MATLAB code that will take the Noraxon .csv files and convert them to all the OpenSim files you need to run inverse dynamics. Note that this is for a single .csv, feel free to edit to parse through all files in a directory, etc)
   
   d. NoraxonTrajectoriesAndForces.mat (This helps rename the columns from the Noraxon data to be more manageable)
   
   e. Data from whatever you measured with the Noraxon IMUs, in .csv form. Example data is included as 'ExampleNoraxonData.csv'
   
   f. LFB_Noraxon_Markers.xml (Marker data for this model and from this equipment)

   g. GRFXMLTemplate.xml (blank template file that will be utilized in the OpenSim inverse dynamics tool) 
   
   h. The following .xml files, all of which are templates that will be converted into OpenSim setup files to quickly population fields:

         i. ScaleXMLTemplate.xml
   
         ii. IK_Template.xml
   
         iii. ID_Template.xml

   i. LFB_model.osim (the base OpenSim model)
   
3. Run NoraxonToOpenSimSINGLE.m. Output should be a .trc file, a _GRF.mot file, Scale_Setup.xml, IK_Setup.xml, ID_Setup.xml. The other .txt files created can be deleted.
4. In OpenSim now - load LFB_model.osim model
5. Scale model: load Scale_Setup.xml and run. Once done, remove LFB_model.osim
6. Inverse Kinematics: load IK_Setup.xml and run
7. Inverse Dynamics: load ID_Setup.xml and run
8. Right-click IDResults and save to a .sto file
