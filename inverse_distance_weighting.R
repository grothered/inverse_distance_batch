## Code to perform inverse distance weighting on some observations, and write
## out some contours (and some rasters) with the output

## INPUT DATA: A table of: [name, x, y, v1, v2, v3, v4,.....]

# The code will make the output for every column in the table following [ name,x,y ].
# It will name these outputs based on the names of the columns

# To run it, open R in the same directory as the script (and the data), and type:
# source('inverse_distance_weighting.R')

# INPUT VARIABLES
datafile='station_rainfall.csv' # The input data file name. 
output_directory='shapefiles'  # Output data will be written to this folder. The folder will be created if it doesn't exist

## MORE INPUT VARIABLES THAT MIGHT LESS OFTEN CHANGE
spatial_ref=4326 ## EPSG Code for wgs84 lat/long -- change for other projections
n_cols=200 # Number of cells along the long direction in the output raster
n_rows=200 # Number of cells along the lat direction in the output raster
range_extension=0.1 # Gridded data has x-range and y-range = data range +- (range_extension)*data_range
inverse_distance_power=2.0 # Weights according to 1/d^(inverse_distance_power). 2.0 is a common choice. Higher values lead to a less 'peaky' interpolation. Very large values approach nearest-neighbour interpolation, but will have numerical problems.
nearest_neighbour=FALSE # If TRUE, do nearest-neighbour interpolation, ignoring the inverse_distance_power
contour_levels=NULL # Vector of contour levels, e.g. c(10, 20, 25, 50, 100). If NULL, values will be automatically selected based on the range of the data
roundoff_threshold=1.0e-08 # If the range in the raster values is < roundoff threshold, then it is treated as a constant. In this case, no contours are created. This is needed because trying to make contours from a constant raster causes an error. 

## MAIN CODE
dir.create(output_directory)

# Libraries to deal with spatial data
library(sp)
library(rgdal)
library(raster)
library(gstat)

station_dat=read.csv(datafile,header=T)

# Remove rows with only missing data
keep=station_dat[,1]!=""
station_dat=station_dat[keep,]

# Record dimension of data
l=dim(station_dat)

# Get proj4string for this projection
EPSG_data=make_EPSG()
p4string=EPSG_data[which(EPSG_data[,1]==spatial_ref),3]

# Put into R's spatialPoints data structures
station_dat_pts=SpatialPoints(station_dat[,3:2], proj4string=CRS(p4string))

# Ranges for plot = range of data + 10%
dx=diff(range(station_dat[,3]))*range_extension
dy=diff(range(station_dat[,2]))*range_extension
long_MIN=min(station_dat[,3]) -dx
long_MAX=max(station_dat[,3]) +dx
lat_MIN=min(station_dat[,2]) -dy
lat_MAX=max(station_dat[,2]) +dy

# Treat the case of nearest_neighbour interpolation
if(nearest_neighbour){
    n_max=1
}else{
    n_max=Inf
}

# Points to predict on a grid
stations_grid=expand.grid(seq(long_MIN, long_MAX, len=n_cols), seq(lat_MIN, lat_MAX, len=n_rows))
names(stations_grid)=names(station_dat)[3:2]

stations_grid_sp=SpatialPoints(stations_grid, proj4string=CRS(p4string))
stations_grid_px=SpatialPixels(stations_grid_sp)

# Get ready to plot
myrast=raster(stations_grid_px)
pdf('plot_series.pdf',width=5,height=5,onefile=T)

# Loop over columns 4:n in the table
for(i in 4:(length(station_dat[1,]))){
#for(i in 4:4){
    print(paste('Column ,' ,i))
    
    varname=names(station_dat)[i]

    # Fit model
    temp_data=station_dat[,c(3,2,i)]
    names(temp_data)[3]='X1'
    stations.gstat<-gstat(id='X1', formula = X1 ~ 1, 
                          locations = ~ LON + LAT, 
                          data=temp_data, set = list(idp = inverse_distance_power), nmax=n_max)
    predvals=predict(stations.gstat, stations_grid)[,3]
   
    output_store=t(matrix(predvals, nrow=n_rows, ncol=n_cols,byrow=F))
    
    # Make raster + contour
    values(myrast)=output_store
    myrast=flip(myrast,direction='y')


    do_contour=abs(myrast@data@min - myrast@data@max)>roundoff_threshold # Skip the contour for a constantraster. We need a round-off criterion
    if(do_contour){
        if(is.null(contour_levels)){
            mycontour=rasterToContour(myrast)
        }else{
            mycontour=rasterToContour(myrast,levels=contour_levels) 
        }
    }
    # Plot 
    plot(myrast, main=names(station_dat)[i])
    if(do_contour) plot(mycontour,add=T)
    plot(station_dat_pts,add=T,cex=0.2,pch=19)

    # Write shapefile
    shplayer=names(station_dat)[i]
    shp_dsn=paste(output_directory,'/',shplayer,sep="")

    if(do_contour) writeOGR(mycontour, dsn=shp_dsn, layer=shplayer, 
                            driver='ESRI Shapefile', overwrite=T, delete_dsn=T)
    writeRaster(myrast, 
                filename=paste(output_directory,'/',names(station_dat)[i],'.tif',sep=""),
                driver='GTiff', overwrite=T)
}

dev.off()
