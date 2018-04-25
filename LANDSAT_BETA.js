// BAMS for Landsat implementation // Dhemerson Conciani (dh.conciani@gmail.com)
// Originally developed by Aitor Bastarika for Sentinel 2

// Define time interval to composite scenes
var date_start="1985-01-01";
var date_end="2011-12-30";

// If you wanna run the algorithm to a specific area, create var named "Region"
// By default, the screen extent are taked
var Region=ee.Geometry.Rectangle(Map.getBounds());

// Define the Dataset to use
//var image_collection_type="LANDSAT/LC08/C01/T1_SR"; // OLI - Landsat 8
var image_collection_type="LANDSAT/LT05/C01/T1_SR"; // TM  - Landsat 5

// Define bands
if (image_collection_type=="LANDSAT/LC08/C01/T1_SR"){
  var Band_BLUE='B2'
  var Band_RED='B4'
  var Band_NIR='B5'
  var Band_SSWIR='B6'
  var Band_LSWIR='B7'
  var resolution=30
} 

if (image_collection_type=="LANDSAT/LT05/C01/T1_SR"){
  var Band_BLUE='B1'
  var Band_RED='B3'
  var Band_NIR='B4'
  var Band_SSWIR='B5'
  var Band_LSWIR='B7'
  var resolution=30
}

// Calc spectral index
var add_NBR2 = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_LSWIR + '") - b("' + Band_SSWIR + '")) / (b("' + Band_LSWIR + '") + b("' + Band_SSWIR+ '"))').rename('NBR2'))
};
var add_NBR = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_LSWIR + '") - b("' + Band_NIR + '")) / (b("' + Band_LSWIR + '") + b("' + Band_NIR+ '"))').rename('NBR'));
};
var add_MIRBI = function(image) {
  return image.addBands(image.expression('10* float(b("' + Band_LSWIR + '"))/10000 - 9.8 * float( b("' + Band_SSWIR + '")) /10000 + 2').rename('MIRBI'));
};

// Add date variable
var add_date = function(image) {
  //print (image.metadata('system:time_start'))
  return image.addBands(image.metadata('system:time_start'))
}

// Burned bands
var addBurnedQualityBand = function(image) {
    return image.addBands(image.normalizedDifference([Band_LSWIR, Band_NIR]).rename('BurnedQ'))
};

var addNotBurnedQualityBand = function(image) {
    return image.addBands(image.normalizedDifference([Band_NIR, Band_RED]).rename('NotBurnedQ'))
};

// CLOUD MASK ALGORITHM
// Empirical using Blue Band (for Clouds) and SWIR BAnd (for shadows)
var maskNotValid = function(image) {
  var cloudMask = image.select(Band_BLUE).lt(1500).and(image.select(Band_LSWIR).gt(600))
  return image.updateMask(cloudMask);
}

// Define Maximum NBR
function get_post_composite (date_start, date_end){
  var image_collection = ee.ImageCollection(image_collection_type)
    .filterDate(date_start, date_end)
    .map(addBurnedQualityBand)
    .map(maskNotValid)
    .map(add_MIRBI)
    .map(add_NBR2)
    .map(add_NBR)
    .map(add_date)
    .filterBounds(Region);
  var composite = image_collection.qualityMosaic('BurnedQ');
  
  return composite
}

// Define Maximum NDVI
function get_pre_composite (date_start, date_end,image_collection_type){
  var image_collection = ee.ImageCollection(image_collection_type) 
    .filterDate(date_start, date_end)
    .map(addNotBurnedQualityBand)
    .map(maskNotValid)
    .map(add_MIRBI)
    .map(add_NBR2)
    .map(add_NBR)
    .map(add_date)
    .filterBounds(Region);
  var composite = image_collection.qualityMosaic('NotBurnedQ');
  return composite
}

var post_composite=get_post_composite (date_start, date_end,image_collection_type)
var pre_composite=get_pre_composite (date_start, date_end,image_collection_type)
var subtract_post_pre_composite=post_composite.subtract(pre_composite)

Map.addLayer( pre_composite , {bands: [Band_LSWIR, Band_NIR, Band_RED], min: [0, 0, 0], max: [2500, 2500, 2500]}, 'NDVI maximum composite');
Map.addLayer( post_composite , {bands: [Band_LSWIR, Band_NIR, Band_RED], min: [0, 0, 0], max: [2500, 2500, 2500]}, 'NBR maximum composite');

//Seeds
var thr_AMIRBI = ee.Number(subtract_post_pre_composite.select('MIRBI').reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:Burned,scale:resolution}).get('MIRBI'));
var thr_ANBR2 = ee.Number(subtract_post_pre_composite.select('NBR2').reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:Burned,scale:resolution}).get('NBR2'));
var thr_ANBR = ee.Number(subtract_post_pre_composite.select('NBR').reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:Burned,scale:resolution}).get('NBR'));

var thr_MIRBI = ee.Number(post_composite.select('MIRBI').reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:Burned,scale:resolution}).get('MIRBI'));
var thr_NBR2 = ee.Number(post_composite.select('NBR2').reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:Burned,scale:resolution}).get('NBR2'));
var thr_NBR = ee.Number(post_composite.select('NBR').reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:Burned,scale:resolution}).get('NBR'));

var BA_seed= subtract_post_pre_composite.select('MIRBI').gt(thr_AMIRBI).and((subtract_post_pre_composite.select('NBR2').gt(thr_ANBR2)).and(subtract_post_pre_composite.select('NBR').gt(thr_ANBR)).and(post_composite.select('MIRBI').gt(thr_MIRBI)).and(post_composite.select('NBR2').gt(thr_NBR2)).and(post_composite.select('NBR').gt(thr_NBR)))
BA_seed=BA_seed.updateMask(BA_seed.eq(1))

//Second
var thr_AMIRBI = ee.Number(subtract_post_pre_composite.select('MIRBI').reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:Burned,scale:resolution }).get('MIRBI'));
var thr_ANBR2 = ee.Number(subtract_post_pre_composite.select('NBR2').reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:Burned,scale:resolution}).get('NBR2'));
var thr_ANBR = ee.Number(subtract_post_pre_composite.select('NBR').reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:Burned,scale:resolution}).get('NBR'));

var thr_MIRBI = ee.Number(post_composite.select('MIRBI').reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:Burned,scale:resolution}).get('MIRBI'));
var thr_NBR2 = ee.Number(post_composite.select('NBR2').reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:Burned,scale:resolution}).get('NBR2'));
var thr_NBR = ee.Number(post_composite.select('NBR').reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:Burned,scale:resolution}).get('NBR'));


var BA_second= subtract_post_pre_composite.select('MIRBI').gt(thr_AMIRBI).and((subtract_post_pre_composite.select('NBR2').gt(thr_ANBR2)).and(subtract_post_pre_composite.select('NBR').gt(thr_ANBR)).and(post_composite.select('MIRBI').gt(thr_MIRBI)).and(post_composite.select('NBR2').gt(thr_NBR2)).and(post_composite.select('NBR').gt(thr_NBR)))
var BA_second= post_composite.select('NBR2').gt(thr_NBR2).and(post_composite.select('NBR').gt(thr_NBR))
BA_second=BA_second.updateMask(BA_second.eq(1))

var BA_vectors = BA_second.addBands([BA_seed,post_composite.select('system:time_start')]).reduceToVectors({
  geometry: Region,
  maxPixels: 10e12,
  scale: resolution,
  geometryType: 'polygon',
  eightConnected: true,
  labelProperty: 'zone',
  reducer: ee.Reducer.count().combine(ee.Reducer.median().unweighted())
});


BA_vectors=BA_vectors.filter(ee.Filter.gt('count', 0));
BA_vectors=BA_vectors.map(function(f) { return f.set("date",ee.Date(f.get('median')).format('YYYY-MM-dd'))})

// Show Results
Map.addLayer( BA_vectors, {}, 'Result');

// Export to GDrive
Export.table.toDrive({
  collection: BA_vectors,
  description:'BA_Vector',
  fileFormat: 'KML'
});
