# mIHC_QuPath_analysis
QuPath pipeline for Claire


Steps: 

1. convert_fused_image.ijm\n
Requirements: Bio-Formats plugin for imageJ\n
Convert the timepoints in the fused image to colour channels.\n
-Make a folder with all the fused images in it\n
-run the script and when prompted, select the folder\n 

2. Need Bio-Formats command line tools\n
-Download the zip file from here:	https://docs.openmicroscopy.org/bio-formats/6.6.0/users/comlinetools/ \n
-unzip \n
-cd into the unzipped folder. \n
Run the following in terminal \n
$ ./tiffcomment /path/to/image | ./xmlindent >> /output/folder/{name}.xml \n

3. 
-open the image you saved from 1 in QuPath. \n
-open rename_channels.groovy in script \n
-Change the variable filename to the location of the xml file \n
-run the file \n

