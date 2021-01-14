/*
 * Author: Kenta
 * Description: Converts images with CZT dimension 1x1x15 to 15x1x1  
 * 				ie. changes timepoints into colour channels
 */
 
macro "Nuclei Binary Mask" {
	setBatchMode(true);

	// Input is the folder containing all the images
	input = getArgument(); 
	if(input =="");{
		input = getDirectory("Select the folder of image sets");
	}
	
	// Each image
	pdir = getFileList(input); 

	// Create output folder 
	output = input+"Processed/";
	if(!File.exists(output)) {
		File.makeDirectory(output);
	}

	// Process each image
	for (m=0; m<pdir.length; m++){
		
		if(endsWith(pdir[m], ".tif")){
			// image filename
			imageName = output 
						+ "ColourChannel_"
						+ pdir[m];
			processFolder(input+pdir[m], imageName);
		}
	}
	print("Done!");
}
	

/*
 * Create a hyperstack from the stack using Bio-formats plugin
 * @params input The location of the image to be converted
 * @params output The location of the converted image
 */
function processFolder(input, output) {

	// Get image dimensions
	run("Bio-Formats Macro Extensions");
	Ext.setId(input);
	Ext.getSizeX(width);
	Ext.getSizeY(height);

	// Open image using BF importer
	run("Bio-Formats Importer", 
		"open="+input+" "+
		"color_mode=Default crop rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT "+
		"x_coordinate_1=0 y_coordinate_1=0 "+
		"width_1="+width+" "+
		"height_1="+height);
	
	// Get Dimensions for the image
	getDimensions(width, height, channels, slices, frames);

	// If image 
	if (frames <= 1){
		print(input+" did not have enough frames to process");
		close()
		return;
		}
	
	// Convert to hyperstack
	run("Stack to Hyperstack...",
		"order=xyczt(default) channels="+frames+" slices=1 frames=1 display=Color");

	// Save
	saveAs("Tiff", output);
	close();
	}