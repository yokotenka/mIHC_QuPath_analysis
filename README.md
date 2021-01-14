# mIHC_QuPath_analysis
QuPath pipeline for Claire


Steps: 

1.convert_fused_image.ijm.    
Requirements: Bio-Formats plugin for imageJ.   
Description: Convert the timepoints in the fused image to colour channels.
Method:   
-Make a folder with all the fused images in it.  
-run the script and when prompted, select the folder.  

2. Get metadata xml file
Requirements: Need Bio-Formats command line tools.  
Description: Use the bio-formats tools to extract the metadata.   
Method:  
-Download the zip file from here:	https://docs.openmicroscopy.org/bio-formats/6.6.0/users/comlinetools/   
-unzip    
-cd into the unzipped folder.   
Run the following in terminal   
$ ./tiffcomment /path/to/image | ./xmlindent >> /output/folder/{name}.xml   

3. rename_channels.groovy        
Method:    
-open the image you saved from 1 in QuPath. 
-open rename_channels.groovy in script     
-Change the variable filename to the location of the xml file     
-run.  

