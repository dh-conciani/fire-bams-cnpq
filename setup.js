//BAMS methodology implemented
//Aitor Bastarrika aitor.bastarrika@ehu.eus
//L5 and 7 Implemented by Dhemerson Conciani (dh.conciani@gmail.com)

//******************* Start parameters
//Composite start and end date
var date_start='2016-01-01' 
var date_end='2017-01-01'
// If you need a specific region, define a geometry named Region
// By default, the Map extent is taken
var Region=ee.Geometry.Rectangle(Map.getBounds())
//
var image_collection_type="COPERNICUS/S2"
//var image_collection_type="LANDSAT/LC08/C01/T1_SR"


//************ end parameters
if (image_collection_type=="LANDSAT/LC08/C01/T1_SR"){
  var Band_RED='B4'
  var Band_NIR='B5'
  var Band_SSWIR='B6'
  var Band_LSWIR='B7'
  var resolution=30
} 
if (image_collection_type=="COPERNICUS/S2"){
  var Band_RED='B4'
  var Band_NIR='B8'
  var Band_SSWIR='B11'
  var Band_LSWIR='B12'
  var resolution=20
} 

var add_NBR2 = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_LSWIR + '") - b("' + Band_SSWIR + '")) / (b("' + Band_LSWIR + '") + b("' + Band_SSWIR+ '"))').rename('NBR2'))
};
var add_NBR = function(image) {
  return image.addBands(image.expression('float(b("'+ Band_LSWIR + '") - b("' + Band_NIR + '")) / (b("' + Band_LSWIR + '") + b("' + Band_NIR+ '"))').rename('NBR'));
};
var add_MIRBI = function(image) {
  return image.addBands(image.expression('10* float(b("' + Band_LSWIR + '"))/10000 - 9.8 * float( b("' + Band_SSWIR + '")) /10000 + 2').rename('MIRBI'));
};

var add_date = function(image) {
  //print (image.metadata('system:time_start'))
  return image.addBands(image.metadata('system:time_start'))
}

var addBurnedQualityBand = function(image) {
    return image.addBands(image.normalizedDifference([Band_LSWIR, Band_NIR]).rename('BurnedQ'))
};

var addNotBurnedQualityBand = function(image) {
    return image.addBands(image.normalizedDifference([Band_NIR, Band_RED]).rename('NotBurnedQ'))
};

function maskNotValid(image) {
  //print (image_collection_type)
  if (image_collection_type=="COPERNICUS/S2") {
    var qa = image.select('QA60');
    var blue=image.select('B1')
    // Bits 10 and 11 are clouds and cirrus, respectively.
    var cloudBitMask = Math.pow(2, 10);
    var cirrusBitMask = Math.pow(2, 11);
    // Both flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).or(blue.gt(1500)).eq(0));
    
  // Return the masked and scaled data.
  }
 if (image_collection_type=="LANDSAT/LC08/C01/T1_SR") {
  // binary mask, no clouds and cloud shadows
  var mask = image.select('cfmask').neq([4, 2]).reduce(ee.Reducer.allNonZero());
 }
  return image.updateMask(mask);
}

// Define Maximum NBR
function get_post_composite (date_start, date_end){
  var image_collection = ee.ImageCollection(image_collection_type)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 90))
    .filterDate(date_start, date_end)
    .map(addBurnedQualityBand)
    .map(maskNotValid)
    .map(add_MIRBI)
    .map(add_NBR2)
    .map(add_NBR)
    .map(add_date)
    .filterBounds(Region);
  //print (image_collection)
  var composite = image_collection.qualityMosaic('BurnedQ');
  
  return composite
}

// Define Maximum NDVI
function get_pre_composite (date_start, date_end,image_collection_type){
  var image_collection = ee.ImageCollection(image_collection_type) 
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 90))
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


//var BA_second= subtract_post_pre_composite.select('MIRBI').gt(thr_AMIRBI).and((subtract_post_pre_composite.select('NBR2').gt(thr_ANBR2)).and(subtract_post_pre_composite.select('NBR').gt(thr_ANBR)).and(post_composite.select('MIRBI').gt(thr_MIRBI)).and(post_composite.select('NBR2').gt(thr_NBR2)).and(post_composite.select('NBR').gt(thr_NBR)))
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

Map.addLayer( BA_vectors, {}, 'Result');


Export.table.toDrive({
  collection: BA_vectors,
  description:'BA_Vector',
  fileFormat: 'KML'
});





