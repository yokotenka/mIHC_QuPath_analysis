# mIHC_QuPath_analysis
QuPath pipeline for Claire


## For Importing image into QuPath:

1. convert_fused_image.ijm.    
Requirements: Bio-Formats plugin for imageJ.   
Description: Convert the timepoints in the fused image to colour channels.       
Method:   
-make a folder with all the fused images in it.  
-run the script and when prompted, select the folder.  

2. Get metadata xml file.   
Requirements: Need Bio-Formats command line tools.  
Description: Use the bio-formats tools to extract the metadata.   
Method:  
-download the zip file from here:	https://docs.openmicroscopy.org/bio-formats/6.6.0/users/comlinetools/   
-unzip    
-cd into the unzipped folder.   
-run the following in terminal:
      $ ./tiffcomment PATH_TO_ORIGINAL_FUSED_IMAGE_BEFORE_CONVERSION | ./xmlindent >> /output/folder/{name}.xml   

3. rename_channels.groovy            
Method:       
-open the image you saved from 1 in QuPath.    
-open rename_channels.groovy in script         
-change the variable filename to the location of the xml file         
-run.      


## For phenotyping:

1. Complete above steps

2. Select a region of interest and run stardist as shown here:
      https://github.com/ninatubau/QuPath_scripts

3. Run the SPIAT_threshold.groovy file