// This script will:
// 1. take czis from input folder
// 2. create ROIs with stardist plugin 
// 3. filter out ROIs too small to be cell
// 4. measure values in ch2 and ch3
// 5. ask user to create background ROIs, and save measurements to sep file
// 6. put results into one csv per image and save to output folder

close("*");

print("Choose the folder containing images you would like to process.");
input_folder = getDirectory("Choose the folder containing images you would like to process.");
print("Choose where you would like to save your measurements.");
output_folder = getDirectory("Choose where you would like to save your measurements.");
list_of_images = getFileList(input_folder);

// Change this to the number of channels you have 
expected_num_of_channels = 4

//Change this to the channel with NEUN stain
neun_channel = 1

// Minimum area threshold
minSize = 50;  

run("Set Measurements...", "area mean integrated stack limit add decimal=3");


for (i = 0; i < lengthOf(list_of_images); i++) { 
	
	raw_image = list_of_images[i];
	image_base_name = replace(raw_image, ".czi", "");

	
	// This will ensure the program only opens .czi files
	if (endsWith(raw_image, ".czi")){
		
		run("Bio-Formats Windowless Importer", "open=" + input_folder + raw_image);
		
		Stack.getDimensions(width, height, channels, slices, frames);
		
		// If your image has more than one slice, will create Max Intensity z-projection
		if (slices > 1){
			run("Z Project...", "projection=[Max Intensity]");
		}
		
		// Make sure your image has the expected number of channels
		if (channels == expected_num_of_channels){
			
			
			
			//Create split image with only NeuN channel for StarDist
			run("Duplicate...", "title=NEUN_Channel duplicate channels=1");
			selectImage("NEUN_Channel");
			
			
			//CREATE ROIs with StarDist//
			imageTitle = getTitle();
			print("Running StarDist on: " + imageTitle);
			run("StarDist 2D", "model=Versatile (fluorescent nuclei) " +
    			"input=" + imageTitle + " normalize percentileBottom=1.0 percentileTop=100.0 " +
    			"probThresh=0.479071 nmsThresh=0.3 output=Both excludeBoundary=2");
			
			// Go back to image with all 4 channels
			selectImage("MAX_" + raw_image);
			roiManager("Show All");
			
			//Remove ROIs with area smaller than threshold set at top of code
			n = roiManager("Count");
			roiIndicesToDelete = newArray(); // will hod indices of ROIs to be deleted
			for (j = 0; j < n; j++) {  
			    roiManager("Select", j);
			    roiManager("measure");
			    area = getResult("Area", j);
			    if (area < minSize) {
			        roiIndicesToDelete = Array.concat(roiIndicesToDelete, j);
			        print(area);

			    }
			}

			
			// delete small ROIs in reverse order to avoid changing indices
			Array.reverse(roiIndicesToDelete);
			for (l = 0;  l < lengthOf(roiIndicesToDelete); l++) {
			    roiManager("Select", roiIndicesToDelete[l]);
			    roiManager("Delete");
			}
			
			//save ROIs
			saveAs("ROI Manager", output_folder+"ROIs_"+image_base_name);

			close("Results");
			
			//go back to unsplit .czi
			selectImage("MAX_" + raw_image);
			n = roiManager("count"); // get new updated value
			
			roiManager("Show All"); 
			for (j = 0; j < n; j++) {
		  	 	roiManager("Select", j);
		  	 	Stack.setChannel(2);
				run("Measure");
				Stack.setChannel(3);
				run("Measure");
				
			}
			saveAs("Results", output_folder+"Results_"+image_base_name+".csv");
			
		} 
		close("Results");
		close("ROI Manager");
		
		
		//get background measurements
		
		selectImage("MAX_" + raw_image);
		Stack.setChannel(1);
		waitForUser("Make 5 circles about the same size as your cells in the background of the NEUN channel. Save each one as an ROI.");
		
		
		roiManager("Show All");
		n = roiManager("Count");
		for (j = 0; j < n; j++) {
			roiManager("Select", j);
		  	Stack.setChannel(2);
			run("Measure");
			Stack.setChannel(3);
			run("Measure");
		}
		saveAs("Results", output_folder+"BG_"+image_base_name+".csv");
		close("Results");
		close("ROI Manager");
		close("*");
	}