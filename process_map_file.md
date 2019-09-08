# Guide to processing maps from GADM to make cartograms
![Cartogram Prep Flowchart](sample_images/flowchart.png?raw=true "Cartogram Preparation Flowchart")
This is the process for preparing map data to generate cartograms.

## Table Of Contents
- [Step 1: Download map from GADM](#step-1-download-map-from-gadm)
- [Step 2: Process map in Mapshaper](#step-2-process-map-in-mapshaper)
  * [Step 2.1: Projection](#step-21-projection)
  * [Step 2.2: Simplification](#step-22-simplification)
  * [Step 2.3: Export](#step-23-export)

## Step 1: Download map from GADM
First, visit https://gadm.org/download_country_v3.html and select the country. Click on "Shapefile", which will download a zip file.
You can also obtain the map from another website. Shapefiles, GeoJSON and TopoJSON are supported by Mapshaper.

## Step 2: Process map in Mapshaper
Then, visit https://mapshaper.org/ and import the zip file.
Click on the correct layer, which should be `gadm36_[country]_1` for a map downloaded from GADM, as seen in the image below:

![Mapshaper Layer Selection](sample_images/select_layer.png?raw=true "Mapshaper Layer selection")

### Step 2.1: Projection
You will need to know the desired projection for this map. Then, click on console (upper right corner) and type the following command to project the map, replacing `[spatial_reference]` with your chosen spatial reference system:
```
$ -proj +init=EPSG:[spatial_reference]
```
Check to make sure that the map is projected correctly.

![Mapshaper Projection](sample_images/projection.png?raw=true "Mapshaper Projection")

### Step 2.2: Simplification
Type the following command to find out the total number of vertices.
```
$ -simplify 100% stats
```

![Mapshaper Simplify Stats](sample_images/simplify1.png?raw=true "Mapshaper Simplify Stats")

The *unique coordinate locations* tells you the number of vertices. Use the following equation to calculate the percentage required.
> [percentage] = 100 * [desired no. of vertices] / [unique coordinate locations]

So, if you would like 50,000 verticies, then the percentage will be 100 * 50,000 / 573,076 ≈ **9**
We recommend 50,000 vertices for an optimal balance between time, space and resolution.

Then, type the following command to execute the simplification, replacing `[percentage]` with the required percentage:
```
$ -simplify [percentage]% stats
```

![Mapshaper Simplify](sample_images/simplify2.png?raw=true "Mapshaper Simplify")

### Step 2.3: Export
In Mapshaper, export the map using the following command, replacing `[country_name]` accordingly.
```
$ -o format=geojson [country_name].json
```

Now, continue following the steps here to generate the cartogram: https://github.com/Flow-Based-Cartograms/go_cart#running-the-cartogram-generator
