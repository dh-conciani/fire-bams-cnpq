// BAMS for Landsat collections (TM/ETM/OLI)
// Dhemerson Conciani; Daniel Borini; Swanni Alvarado 
// Originally developed in GEE by Aitor Bastarika for Sentinel 2
// var Burned= Training polygons

// Define time interval to composite scenes
var date_start="2000-02-01";
var date_end="2000-12-30";

// define mask functions
// landsat 5 and landsat 7
var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2);
};

// landsat 8
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

// If you wanna run the algorithm to a specific area, create var named "Region"
// By default, the screen extent are taked

// import collections
// landsat 5 - TM
var L5coll = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR') 
.select(['B1', 'B3', 'B4', 'B5', 'B7', 'pixel_qa']) 
.filterBounds(Region)
.filterDate(date_start,date_end)
.map (cloudMaskL457);

// landsat 7 - ETM
var L7coll = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
.select(['B1', 'B3', 'B4', 'B5', 'B7', 'pixel_qa'])
.filterBounds(Region)
.filterDate(date_start,date_end)
.map (cloudMaskL457);

// landsat 8 - OLI
var L8coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
.filterBounds(Region)
.filterDate(date_start,date_end)
.map(maskL8sr)
.map(function(image){
  return image.rename(['blueblue', 'B1', 'green', 'B3', 'B4', 'B5', 'B7', 'brigtemp', 'B8', 'B9', 'pixel_qa', 'B11']);//padroniza a B7 
})
.select(['B1', 'B3', 'B4', 'B5', 'B7', 'pixel_qa']);

// merge collections 
var image_collection_type = ee.ImageCollection(L5coll.merge(L7coll).merge(L8coll));
//var image_collection_type = ee.ImageCollection(L5coll.merge(L7coll));
print (image_collection_type)

// define band names
if (image_collection_type){
  var Band_BLUE='B1'
  var Band_RED='B3'
  var Band_NIR='B4'
  var Band_SSWIR='B5'
  var Band_LSWIR='B7'
  var resolution=30
}

// Calc spectral index
var add_NBR2 = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_LSWIR + '") - b("' + Band_SSWIR + '")) / (b("' + Band_LSWIR + '") + b("' + Band_SSWIR + '"))').rename('NBR2'))
};
var add_NBR = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_LSWIR + '") - b("' + Band_NIR + '")) / (b("' + Band_LSWIR + '") + b("' + Band_NIR + '"))').rename('NBR'));
};
var add_MIRBI = function(image) {
  return image.addBands(image.expression('10* float(b("' + Band_LSWIR + '"))/10000 - 9.8 * float( b("' + Band_SSWIR + '")) /10000 + 2').rename('MIRBI'));
};
var add_NDVI = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_NIR + '") - b("' + Band_RED + '")) / (b("' +  Band_NIR + '") + b("' + Band_RED + '"))').rename('NDVI'))
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

// empiric mask
var maskNotValid = function(image) {
  var cloudMask = image.select(Band_BLUE).lt(1500).and(image.select(Band_LSWIR).gt(600))
  return image.updateMask(cloudMask);
}


// Define Maximum NBR
function get_post_composite (date_start, date_end){
  var image_collection = ee.ImageCollection(image_collection_type)
    .filterDate(date_start, date_end)
    .map(maskNotValid)
    .map(addBurnedQualityBand)
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
    .map(maskNotValid)
    .map(addNotBurnedQualityBand)
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

//clip
var post_composite = post_composite.clipToCollection(Region)
var pre_composite = pre_composite.clipToCollection(Region)

//plot
Map.addLayer( pre_composite , {bands: [Band_LSWIR, Band_NIR, Band_RED], min: [0, 0, 0], max: [4000, 4000, 4000]}, 'NDVI maximum composite');
Map.addLayer( post_composite , {bands: [Band_LSWIR, Band_NIR, Band_RED], min: [0, 0, 0], max: [4000, 4000, 4000]}, 'NBR maximum composite');
//Map.addLayer (Region, {color: 'EBD7D3'}, 'UCs');
Map.addLayer (burned, {color: '#F0340B'}, 'Training BA');

// Defining Predictors
var MIRBI= post_composite.select('MIRBI').rename('MIRBI')
var NBR= post_composite.select('NBR').rename('NBR')
var NBR2= post_composite.select('NBR2').rename('NBR2')
var NDVI= post_composite.select('NDVI').rename('NDVI')
var AMIRBI= post_composite.select('MIRBI').subtract(pre_composite.select('MIRBI')).rename('AMIRBI')
var ANBR= post_composite.select('NBR').subtract (pre_composite.select('NBR')).rename('ANBR')
var ANBR2= post_composite.select('NBR2').subtract (pre_composite.select('NBR2')).rename('ANBR2')
var ANDVI= post_composite.select('NDVI').subtract(pre_composite.select('NDVI')).rename('ANDVI') 

//Seed
var thr_AMIRBI = ee.Number(AMIRBI.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('AMIRBI'));
var thr_ANBR2 = ee.Number(ANBR2.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('ANBR2'));
var thr_ANDVI = ee.Number(ANDVI.reduceRegion({reducer:ee.Reducer.percentile([75]),geometry:burned,tileScale: 16,scale:100}).get('ANDVI'));
var thr_ANBR = ee.Number(ANBR.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('ANBR'));
var thr_MIRBI = ee.Number(MIRBI.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('MIRBI'));
var thr_NBR2 = ee.Number(NBR2.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('NBR2'));
var thr_NDVI = ee.Number(NDVI.reduceRegion({reducer:ee.Reducer.percentile([75]),geometry:burned,tileScale: 16,scale:100}).get('NDVI'));
var thr_NBR = ee.Number(NBR.reduceRegion({reducer:ee.Reducer.percentile([25]),geometry:burned,tileScale: 16,scale:100}).get('NBR'));

var BA_seed= AMIRBI.gt(thr_AMIRBI).and(ANBR2.gt(thr_ANBR2)).and((ANDVI.lt(thr_ANDVI)).and(ANBR.gt(thr_ANBR)).and(MIRBI.gt(thr_MIRBI)).and(NBR2.gt(thr_NBR2)).and(NDVI.lt(thr_NDVI)).and(NBR.gt(thr_NBR)))


//Second
var thr_AMIRBI = ee.Number(AMIRBI.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('AMIRBI'));
var thr_ANBR2 = ee.Number(ANBR2.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('ANBR2'));
var thr_ANDVI = ee.Number(ANDVI.reduceRegion({reducer:ee.Reducer.percentile([99]),geometry:burned,tileScale: 16,scale:100}).get('ANDVI'));
var thr_ANBR = ee.Number(ANBR.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('ANBR'));
var thr_MIRBI = ee.Number(MIRBI.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('MIRBI'));
var thr_NBR2 = ee.Number(NBR2.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('NBR2'));
var thr_NDVI = ee.Number(NDVI.reduceRegion({reducer:ee.Reducer.percentile([99]),geometry:burned,tileScale: 16,scale:100}).get('NDVI'));
var thr_NBR = ee.Number(NBR.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:burned,tileScale: 16,scale:100}).get('NBR'));

var BA= AMIRBI.gt(thr_AMIRBI).and(ANBR2.gt(thr_ANBR2)).and(ANDVI.lt(thr_ANDVI).and(ANBR.gt(thr_ANBR)).and(MIRBI.gt(thr_MIRBI)).and(NBR2.gt(thr_NBR2)).and(NDVI.lt(thr_NDVI)).and(NBR.gt(thr_NBR)))

BA_seed=BA_seed.updateMask(BA_seed.eq(1))
BA=BA.updateMask(BA.eq(1))

var BA_vectors = BA.addBands([BA_seed,post_composite.select('system:time_start')]).reduceToVectors({
  geometry: Region,
  geometryType: 'polygon',
  eightConnected: false,
  scale: 30,
  maxPixels: 10e13,
  reducer:  ee.Reducer.count().combine(ee.Reducer.median().unweighted())
});

BA_vectors=BA_vectors.filter(ee.Filter.gt('count', 0));
BA_vectors=BA_vectors.map(get_date)

// plot result
Map.addLayer (BA_vectors)

// export vector
Export.table.toDrive({
  collection: BA_vectors, 
  description: '2000_BA_L7L5_ESEC', 
  folder: 'BA_ESEC'});
