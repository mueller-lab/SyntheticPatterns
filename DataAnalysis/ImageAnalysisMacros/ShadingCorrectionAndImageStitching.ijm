//headless mode macro - 20230723

//args = getArgument();
//
//args= split(args, ',');
//
//indir= args[0];
//
//outdir= args[1];

//Choosing flatField folder
flatDir= getDirectory("flatDir");

//Choosing input folder
indir=getDirectory("indir");

//Choosing output folder
outdir= getDirectory("_outdir");

setBatchMode(true);

imageList=getFileList(indir);

//scle = 3.77442 // scale is 3.77442 um per pixel or 3.77442 mm per 1000 px
//scle = 5.66163 // um/px for 4X objective with 640x480 px image resolution
//scle = 7.54884 // um/px for 4X objective with 480x360 px image resolution
scle = 15.09767 // um/px for 2X objective with 480x360 px image resolution 

//Array.print(imageList);

tList= newArray(5);
for (i=0; i<tList.length; i++){
	tList[i] = i+1;	
}

psList=newArray(1);
for (i=0; i<psList.length; i++){
	psList[i] = i+17;	
}

//psList=newArray("01");

chnlList=  newArray("1","2","4"); // newArray("1","2","4");

//make shading correction

gPS_num = 4;
correctShading(indir, outdir, flatDir, psList, chnlList, tList, gPS_num);


// stitch images to grid
for(j=0; j<psList.length; j++){ //for each position
	ps= String.pad(psList[j], 2);
	shadeDir = indir+"shadeCorrected_XY"+ps;

	for(k=0; k<chnlList.length; k++){ //For each channel
		ch=chnlList[k];
	
		for(i=0; i<tList.length; i++){ // for each timepoint
			t= String.pad(i, 4);
			
			stitch_megatron(t, ps, ch, "30", shadeDir, outdir, scle);	
		}
//		makeMontageOfStitched(t, ch, outdir, scle);
	
	}	
}



//
//timeNum = tList.length
//makeGIF_fromStitched(outdir, "CH1", timeNum);
//makeGIF_fromStitched(outdir, "CH2", timeNum);
//makeGIF_fromStitched(outdir, "CH4", timeNum);

print("All done!");


// functions - 

function stitch_megatron(t, ps, ch, olap, indir1, outdir1, scle) { 
// function description - stitch images for given timepoint, position, and channel
	//Stitch the four images for given Time (T00**), pos (XY**), and Channel (CH*)
	// image name - XY01_00025_CH1_T0010.tif
		
	run("Grid/Collection stitching", "type=[Grid: snake by rows] order=[Right & Down                ] "+
	"grid_size_x=2 grid_size_y=2 tile_overlap="+olap+" first_file_index_i=1 "+
	"directory="+indir1+" file_names=XY"+ps+"_000{ii}_CH"+ch+"_T"+t+".tif output_textfile_name=TileConfiguration.txt "+ // 
	"fusion_method=[Max. Intensity] regression_threshold=0.10 max/avg_displacement_threshold=2.50 "+ // Average
	"absolute_displacement_threshold=3.50 computation_parameters=[Save computation time (but use more RAM)] "+ //compute_overlap
	"image_output=[Fuse and display]");
	
	print("Stitched images at time "+t+", pos "+ps+" chnl "+ch);
	outName= "stiched_"+"Image_T"+t+"_CH"+ch+"_XY"+ps+".tif";
	
	run("Set Scale...", "distance=1000 known="+scle+" pixel=1 unit=mm");
//	run("Scale Bar...", "width=2 height=2 thickness=4 font=40 color=White background=None location=[Lower Right] horizontal bold");
	
	outdir2= outdir1+"/"+"XY"+ps+"_CH"+ch+"/";
	
	if (File.exists(outdir2)==0){
		File.makeDirectory(outdir2);
	}
	
	outName1 = outdir2+outName;
	saveAs("Tiff", outName1);
	close("*");
}

function correctShading(indir, outdir, flatDir, psList, chnlList, tList, gPS_num){
	// this functions opens timelapse stack for each channel and each position 
			
	for(j=0; j<psList.length; j++){ //for each position
			ps= String.pad(psList[j], 2);
			
			shadeDir = indir+"shadeCorrected_XY"+ps;
//			File.makeDirectory(shadeDir);

			if (File.exists(shadeDir)==0){
				File.makeDirectory(shadeDir);
			}
	
		for(k=0; k<chnlList.length; k++){ //For each channel
			//open the flatField image
			ch=chnlList[k];
			
			for(l=0; l<gPS_num; l++){ // for each gridPos
				gPS=  String.pad(l+1, 5);
				
				for(i=0; i<tList.length; i++){ // for each timepoint
					t= String.pad(tList[i], 4);
					
					imjN= indir+"XY"+ps+"/"+"T"+t+"/Image_T"+t+"_XY"+ps+"_"+gPS+"_CH"+ch+".tif"; 
					
					print(imjN);
					open(imjN);
					}
				
				run("Images to Stack", "name=Stack title=XY"+ps+"_"+gPS+"_CH"+ch+" use");
				
				outName= "proStack";

				if (ch!="4") {
					run("Split Channels");
			
					
				if (ch=="1") {
					selectWindow("Stack (green)");
					rename(outName);
					close("Stack (red)");
					close("Stack (blue)");
				}
				else {
					selectWindow("Stack (red)");
					rename(outName);
					close("Stack (green)");
					close("Stack (blue)");
				}
				}
				
				if (ch=="4") {
					selectWindow("Stack");
					rename(outName);
				}

				
				if (ch=="1"){
				flatF = flatDir+"Image_CH1.tif"; // The flat field image for CH1 is not correct so using the one for CH2
				open(flatF);
				run("Split Channels");
				
				selectWindow("Image_CH1.tif (green)");
				rename("flatField");
				close("Image_CH1.tif (red)");
				close("Image_CH1.tif (blue)");			
				}
					
				if (ch=="2"){
				flatF = flatDir+"Image_CH2.tif";
				open(flatF);
				run("Split Channels");
				
				selectWindow("Image_CH2.tif (red)");
				rename("flatField");
				close("Image_CH2.tif (green)");
				close("Image_CH2.tif (blue)");									
				}
				
				if (ch=="4"){
				flatF = flatDir+"Image_CH4.tif";
				open(flatF);					
				rename("flatField");		
				}
				
				
				run("BaSiC ", "processing_stack=[proStack] flat-field=[flatField] dark-field=None shading_estimation=[Skip estimation and use predefined shading profiles] "+
				"shading_model=[Estimate flat-field only (ignore dark-field)] setting_regularisationparametes=Automatic temporal_drift=Ignore "+
				"correction_options=[Compute shading and correct images] lambda_flat=2.0 lambda_dark=2.0"); // lambda_flat=0.50 lambda_dark=0.50
				
//				run("BaSiC ", "processing_stack=[Stack (green)] flat-field=[flatField_Green.tif (green)] dark-field=None shading_estimation=[Skip estimation and use predefined shading profiles] "+
//				"shading_model=[Estimate flat-field only (ignore dark-field)] setting_regularisationparametes=Automatic temporal_drift=Ignore "+
//				"correction_options=[Compute shading and correct images] lambda_flat=0.50 lambda_dark=0.50");
				
				selectWindow("Corrected:proStack");
				rename("XY"+ps+"_"+gPS+"_CH"+ch);
				
				run("Image Sequence... ", "dir="+shadeDir+" format=TIFF name=XY"+ps+"_"+gPS+"_CH"+ch+"_T");
				
				close("*");
				
								
	}}}}
