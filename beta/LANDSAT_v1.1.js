// BAMS for Landsat implementation // Dhemerson Conciani (dh.conciani@gmail.com)
// Originally developed by Aitor Bastarika for Sentinel 2
// var Burned= Training polygons

// Define time interval to composite scenes
var date_start="1985-01-01";
var date_end="1985-12-30";

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
var add_NDVI = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_NIR + '") - b("' + Band_RED + '")) / (b("' +  Band_NIR + '") + b("' + Band_RED+ '"))').rename('NDVI'))
};

//Add date variable
var add_date = function(image) {
    var timestamp = image.metadata('system:time_start');
    return image.addBands(timestamp);
};

var get_date=function(d){
    return d.set("date",ee.Date(d.get('median')).format('YYYY-MM-dd'))
};

// Burned bands
var addBurnedQualityBand = function(image) {
    return image.addBands(image.normalizedDifference([Band_LSWIR, Band_NIR]).rename('BurnedQ'))
};

var addNotBurnedQualityBand = function(image) {
    return image.addBands(image.normalizedDifference([Band_NIR, Band_RED]).rename('NotBurnedQ'))
};

// CLOUD MASK ALGORITHM
// Empirical using Blue Band (for Clouds) and SWIR Band (for shadows)
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
    .map(add_NDVI)
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
    .map(add_NDVI)
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


// Defining Predictors
var NBR= post_composite.select('NBR').rename('NBR')
var NBR2= post_composite.select('NBR2').rename('NBR2')
var NDVI= post_composite.select('NDVI').rename('NDVI')
var ANBR= post_composite.select('NBR').subtract (pre_composite.select('NBR')).rename('ANBR')
var ANBR2= post_composite.select('NBR2').subtract (pre_composite.select('NBR2')).rename('ANBR2')
var ANDVI= post_composite.select('NDVI').subtract(pre_composite.select('NDVI')).rename('ANDVI') 

//Seed
var thr_ANBR2 = ee.Number(ANBR2.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('ANBR2'));
var thr_ANDVI = ee.Number(ANDVI.reduceRegion({reducer:ee.Reducer.percentile([75]),geometry:burned,tileScale: 16,scale:100}).get('ANDVI'));
var thr_ANBR = ee.Number(ANBR.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('ANBR'));
var thr_NBR2 = ee.Number(NBR2.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('NBR2'));
var thr_NDVI = ee.Number(NDVI.reduceRegion({reducer:ee.Reducer.percentile([75]),geometry:burned,tileScale: 16,scale:100}).get('NDVI'));
var thr_NBR = ee.Number(NBR.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('NBR'));

var BA_seed= ANBR2.gt(thr_ANBR2).and((ANDVI.lt(thr_ANDVI)).and(ANBR.gt(thr_ANBR)).and(NBR2.gt(thr_NBR2)).and(NDVI.lt(thr_NDVI)).and(NBR.gt(thr_NBR)))

//Second
var thr_ANBR2 = ee.Number(ANBR2.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('ANBR2'));
var thr_ANIR = ee.Number(ANDVI.reduceRegion({reducer:ee.Reducer.percentile([99]),geometry:burned,tileScale: 16,scale:100}).get('ANDVI'));
var thr_ANBR = ee.Number(ANBR.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('ANBR'));
var thr_NBR2 = ee.Number(NBR2.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('NBR2'));
var thr_NIR = ee.Number(NDVI.reduceRegion({reducer:ee.Reducer.percentile([99]),geometry:burned,tileScale: 16,scale:100}).get('NDVI'));
var thr_NBR = ee.Number(NBR.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('NBR'));

var BA= ANBR2.gt(thr_ANBR2).and((ANDVI.lt(thr_ANDVI)).and(ANBR.gt(thr_ANBR)).and(NBR2.gt(thr_NBR2)).and(NDVI.lt(thr_NDVI)).and(NBR.gt(thr_NBR)))

BA_seed=BA_seed.updateMask(BA_seed.eq(1))
BA=BA.updateMask(BA.eq(1))

var BA_vectors = BA.addBands([BA_seed,post_composite.select('system:time_start')]).reduceToVectors({
  geometry: Region,
  scale: 30,
  geometryType: 'polygon',
  eightConnected: true,
  maxPixels: 10e12,
  reducer:  ee.Reducer.count().combine(ee.Reducer.median().unweighted())
});


BA_vectors=BA_vectors.filter(ee.Filter.gt('count', 0));
BA_vectors=BA_vectors.map(get_date)

Map.addLayer(BA_vectors)

Export.table.toDrive({collection: BA_vectors,  description: 'BAMS_perimeter',  folder: 'BAMS',  fileFormat: 'KML'});
