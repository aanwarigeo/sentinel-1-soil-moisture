/* 
--------------------------
This script outlines a step-by-step workflow for step 1 to estimate field-scale
soil moisture using sentinel-1 GRD SAR data.

1. Calculated Dual-polarimetric Radar Vegetation Index (DpRVIc)
2. Processed Sentinel-1 GRD SAR Backscatter intenstiy in vertical-vertical polarization mode (VVdB). 
3. Stack both DpRVIc and VVdB for the complete timeseries in GeoTIFF format.

For a detailed explanation of each step of the methodology, please refer to the 
research paper by Bhogapurapu et al. (2022).

Author: aanwari
Date: Feb 28, 2025
--------------------------
*/


/* 
--------------------------
1. Define experimental area and map
--------------------------
*/ 

var studyArea = ee.Geometry.Polygon(
  [[[-8.657760938819553,31.426809284220685],
    [-8.655540069755222,31.4196864314692],
    [-8.65012200754514,31.42110555197101],
    [-8.6520961133801,31.427926186352536],
    [-8.657760938819553,31.426809284220685]]], null, false);
  
    
Map.centerObject(studyArea, 16);
Map.addLayer(studyArea, {color: 'red'}, 'extent', false);



/*
-----------------------------------------
2) Cloud Filtering
-----------------------------------------
*/

// Date range
var startDate = '2018-11-01'; 
var endDate   = '2019-06-06';

// Filter ImageCollection to match study's specifications
var s1Collection = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filterDate(startDate, endDate)
        .filter(ee.Filter.eq('instrumentMode', 'IW'))  
        .filter(ee.Filter.eq('resolution', 'H'))  
        .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))  
        .filter(ee.Filter.eq('relativeOrbitNumber_start', 118))
        .select('VV', 'VH', 'angle')  
        .sort('system:time_start', false) 
        .filterBounds(studyArea);

/*
-----------------------------------------
3) Preprocessing, maksing, stacking DpRVIc & S1 Backscatter intensity (VVdb) 
-----------------------------------------
*/

var processingParams = {
  windowSize: 1.5,    // Kernel radius in pixels
  scale: 10           // Spatial resolution in meters
};

var processImage = function(image) {
  // Convert dB to power scale
  var c11 = image.expression(
    '10 ** (VV/10)',
    {'VV': image.select('VV')}
  ).rename('C11');
  
  var c22 = image.expression(
    '10 ** (VH/10)',
    {'VH': image.select('VH')}
  ).rename('C22');

  // Apply mean filter
  var meanFilter = {
    reducer: ee.Reducer.mean(),
    kernel: ee.Kernel.square(processingParams.windowSize, 'pixels'),
    optimization: 'boxcar'
  };
  
  var c11Mean = c11.reduceNeighborhood(meanFilter);
  var c22Mean = c22.reduceNeighborhood(meanFilter);

  // Calculate vegetation index
  var ratio = c22Mean.divide(c11Mean);
  var dpRVIc = ratio.expression(
    '(q * (q + 3)) / ((q + 1) ** 2)',
    {'q': ratio}
  ).rename('DpRVIc');

  // Convert back to dB
  var vvDb = c11Mean.log10().multiply(10).rename('VV');

  // Create masks
  var validMask = c11Mean.gt(c22Mean); // remove all pixel value where c22Mean > c11Mean
  
  var dynamicWorld = ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
    .filterDate(startDate, endDate)
    .filterBounds(studyArea);

  var waterMask = dynamicWorld.select('label')
    .mode().neq(0);  // Exclude water (class 0)
    
  var urbanMask = dynamicWorld.select('label')
    .mode().neq(6); // Exclude urban areas (class 6)
    
  var combinedImage = vvDb.addBands(dpRVIc);

  var maskedImage = combinedImage
    .updateMask(validMask)
    .updateMask(waterMask)
    .updateMask(urbanMask)
    .clip(studyArea);

  // Preserve copy the original metadata
  return maskedImage.copyProperties(image, ['system:id', 'system:time_start']);
};

// Apply processImage function to all s1Collection
var processedCollection = s1Collection.map(processImage);


/*
-----------------------------------------
4) Processed ImageCollection visualization
-----------------------------------------
*/ 

var visParams = {
  bands: ['VV'],
  min: -20,
  max:-5,
  palette: ['654321', 'c2b280', 'ffffff', '94cfff', '0066cc']
};

Map.addLayer(processedCollection.first(), visParams, 'Processed VV');


/*
-----------------------------------------
5) Export Processed ImageCollection to Google drive
-----------------------------------------
*/

var exportData = function(collection, folderName) {
  var count = collection.size().getInfo();
  var imageList = collection.toList(count);
  
  for (var i = 0; i < count; i++) {
    var image = ee.Image(imageList.get(i));
    var imageId = image.get('system:id').getInfo().split('/').pop();
    
    Export.image.toDrive({
      image: image.select(['VV', 'DpRVIc']).toDouble(),
      description: imageId,
      folder: folderName,
      region: studyArea,
      scale: processingParams.scale,
      maxPixels: 1e12,
      fileFormat: 'GeoTIFF',
      formatOptions: {
        cloudOptimized: true
      }
    });
  }
};

// Uncomment the following line to download all images
// exportData(processedCollection, 'Sentinel-1SoilMoistureMorocco');


/* 
--------------------------
End of script
--------------------------
*/ 