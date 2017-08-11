// written by Moritz.Kirschmann@zmb.uzh.ch 2015
// Recursively converts files to TIFF using Bio-Formats reader.
// FlattenFolders is not tested yet.

requires("1.39u");

ext = "lif"; // this variable controls the extension of source files
filter_by_name =true;
filtered_by_name = "erging"; // filtered by this string
split_ch=0;
rgb_out=0;
saveAsSequence= true;
removeLastChars=0;

run("Close All");

inDir = getDirectory("Choose Directory Containing " + ext + " Files ");
outDir = getDirectory("Choose Directory for TIFF Output ");
//outDir = inDir;
setBatchMode(false);
processFiles(inDir, outDir, "");
print("-- Done --");

function processFiles(inBase, outBase, sub) {
	flattenFolders = false; // this flag controls output directory structure
	print("-- Processing folder: " + sub + " --");
	list = getFileList(inBase + sub);
	if (!flattenFolders) File.makeDirectory(outBase + sub);
	for (i=0; i<list.length; i++) {
		path = sub + list[i];
    	upath = toUpperCase(path);
    	if (endsWith(path, "/")) {
      		// recurse into subdirectories
      		processFiles(inBase, outBase, path);
    	} else if (endsWith(upath, "." + toUpperCase(ext))) {
      		print("-- Processing file: " + path + " --");
      		
			// All subimages withing the file list[i] are loaded
     		run("Bio-Formats Importer", "open=[" + inBase + path + "] color_mode=Composite open_all_series view=Hyperstack stack_order=XYCZT");
			//run("Bio-Formats Importer", "open=[" + inBase + path + "] color_mode=Default view=[Standard ImageJ] stack_order=Default");
			//waitForUser("did it open?");
			//run("Bio-Formats Importer", "open='" + inBase + path + "' color_mode=Default view=[Standard ImageJ] stack_order=Default use_virtual_stack");
      		getDimensions(width, height, channels, slices, frames);
     		print("slices "+slices+" channels "+channels +" frames " + frames);

			// run over all opened images which are substacks/images of the LIF and converts them on by one


			//while(stillImagesOpen) {
			numberOfImagesPerImportFile=nImages;
			for (j=numberOfImagesPerImportFile; j>=1; j--) {
    		    selectImage(j);
        
       			print("Image "+j+" of "+nImages +" is "+getTitle);
        		title=getTitle;
				//selectWindow(list[i]+" - "+imageNumber);

				// if the number of channels does not match reorder the stack into a hyperstack
				/*if (channels>1) {
					//slices=slices/n_channels;
					run("Stack to Hyperstack...", "order=xyzct channels=" + channels + " slices=" + slices + " frames=1 display=Color");
					print("hyperstacked!");
				}*/

				//if (indexOf(title, "erg")>0) {
				title=replace(title, "/","_");
				if (removeLastChars>0){
					titleLen=lengthOf(title);
					title= substring(title, 0, titleLen-removeLastChars);
				}
				
      			outFile = outBase + replace(title, "."+ ext, "") + ".tif";
				
      			
				print("outFile: "+outFile);
				
				if (indexOf(outFile, filtered_by_name)<0) {
					if (filter_by_name) {
						print("not processing "+ title + " because it does not contain "+filtered_by_name);
						close();
				
					} 
				}else{
				
				
				if (flattenFolders) outFile = replace(outFile, "/", "_");
      			//run("Bio-Formats Exporter", "save=[" + outFile + "] compression=Uncompressed");
     			//saveAs("Tiff", outFile);
      			title=getTitle;

      			//selectWindow(title);
     		 	if ((split_ch==1) &&(channels>1)) {
      				print("case: 1");
      				run("Duplicate...", "duplicate");
      				rename("split");
					run("Split Channels");
					for (k=1; k<=channels; k++) {
						selectWindow("C"+k+"-"+"split");
						l=0;
						while (File.exists(replace(outFile, ".tif", "_CH_"+k+"_"+l+".tif")) ) {
							l++;
						}
						saveAs("Tiff", replace(outFile, ".tif", "_CH_"+k+"_"+l+".tif"));
						title=getTitle();
						//waitForUser("Channel"+k);
						close();
					}
					selectWindow(title);
					close();
				} else if (rgb_out==1) {
					print("case: 2");
					selectWindow(title);
					run("Duplicate...", "duplicate");
					rename("rgbPrep");
					if (channels>1) {
						//Stack.setDisplayMode("color"); 
			
						for (k=1; k<=channels; k++) {
							//Stack.setDimensions(k, 1, 1)
							Stack.setChannel(k);
					
							setMinAndMax(0, 4095);
							call("ij.ImagePlus.setDefault16bitRange", 12);
							//Stack.setDisplayMode("color");
						}
					} else {
						setMinAndMax(0, 4095);
						call("ij.ImagePlus.setDefault16bitRange", 12);
					}
			
					//call("ij.ImagePlus.setDefault16bitRange", 12);
					run("RGB Color");

					l=0;
					while (File.exists(replace(outFile, ".tif", "_"+l+"_RGB.tif")) ) {
						l++;
					}
					saveAs("Tiff", replace(outFile, ".tif", "_"+l+"_RGB.tif"));
					title=getTitle();
			
					//waitForUser("RGB");
					close();
				} else {
					print("case: 3");
					if  (!saveAsSequence) {
						l=0;
						while (File.exists(replace(outFile, ".tif", "_"+l+".tif")) ) {
							l++;
						}
						if (l==0) {
							print("outFile wo gespeichert wird:"+ outFile);
							saveAs("Tiff", outFile);
							close();
						} 
					} else  {
						newFolderName=outBase+File.separator+replace(replace(title, ext,"")," ","_");
						print("Making a folder at "+ newFolderName);
						File.makeDirectory(newFolderName);
						if (indexOf(newFolderName, " ")>0) print("Spaces in path"); 
						if (File.isDirectory(newFolderName)) print("Dir was created successfully.");
						
						outFile = newFolderName +File.separator+ replace(replace(title, "."+ ext, "")," ", "_")   ;
						print("outFile wo gespeichert wird:"+ outFile);
						//waitForUser("Check Now!");
						run("Bio-Formats Exporter", "save="+outFile+j+".ome.tif write_each_z_section write_each_timepoint write_each_channel compression=Uncompressed"); 
						close();
					}
	
					//title=getTitle();
		
				}



				/*	if (isOpen(list[i]+" - "+imageNumber)) {
				print("saving next sub image "+  imageNumber);
				} else {
				stillImagesOpen=false;
				}
				//run("Close All");
				*/
   		  	 	run("Collect Garbage");
				}
				//close();
    	//run("Close All");
		
	}
	//run("Close All");
}
//waitForUser("alles zu?");
}
}