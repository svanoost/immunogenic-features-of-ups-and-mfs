
// Cell Segmentation
// Set the stain vectors for DAB (CD68) and AP (CD163) manually based on single stains
setColorDeconvolutionStains('{"Name" : "doubleIHC", "Stain 1" : "Hematoxylin", "Values 1" : "0.78559 0.5817 0.21089", "Stain 2" : "DAB", "Values 2" : "0.24506 0.31408 0.91722", "Stain 3" : "AP", "Values 3" : "0.038 0.81091 0.58394", "Background" : " 255 255 255"}');

// Select all ROIs
selectAnnotations();

// Perform cell segmentation using the optical density sum
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImageBrightfield": "Optical density sum",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 3.0,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 0.3,  "maxBackground": 2.0,  "watershedPostProcess": true,  "excludeDAB": false,  "cellExpansionMicrons": 0.0,  "includeNuclei": false,  "smoothBoundaries": true,  "makeMeasurements": true}');

// Cell classification
// Select all ROIs
selectAnnotations();

// Run an object classifer trained to classify cells as CD68+, CD163+, CD68+CD163+ or negative
runObjectClassifier("Myeloid");

// Save the detection measurements per ROI, per sample
saveDetectionMeasurements('../Myeloid_counts/')

// Tissue area measurement
// Select all ROIs
selectAnnotations();

// Run a pixel classifier trained to detect tissue and background, and measure the tissue area per ROI
addPixelClassifierMeasurements("tissue_area", "tissue_area")

// Save the annotation measurements per ROI, per sample
saveAnnotationMeasurements('../tissue_area/')
