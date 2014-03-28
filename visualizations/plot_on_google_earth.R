library(dismo)
library(rhdf5)


plot.interaction.locs<-function(latlon.file,interactions.file,plot.name,event.type,strength.thresh,disp.thresh){
	#INPUTS:
	#latlon.file - [string] full path to the file containing lat / lon data (hdf5 file)
	#interactions.file - [string] full path to the file containing dyadic interactions (csv file)
	
	#read in lat/lon data
	print('reading in lat/lon data...') 
	lats<-h5read(latlon.file,'/lats')
	lons<-h5read(latlon.file,'/lons')
	
	#interpolate latitude and longitude data for better viewing
	#print('interpolating NaNs in lat/lon data...')
	#for(t in 2:dim(lats)[2]){
	#	for(i in 1:dim(lats)[1]){
	#		if( is.na(lats[i,t]) | is.nan(lats[i,t])){
	#			lats[i,t] <- (lats[i,t+1] + lats[i,t-1]) / 2
	#			lons[i,t] <- (lons[i,t+1] + lons[i,t-1]) / 2
	#		}
	#	}
	#}
	
	#read in interaction data
	print('reading in dyadic interactions data...')
	ints <- read.table(interactions.file,header=TRUE,sep=',')
	
	#get map boundaries
	xmin = min(lons,na.rm=T)-.0005
	xmax = max(lons,na.rm=T)+.0005
	ymin = min(lats,na.rm=T)-.0005
	ymax = max(lats,na.rm=T)+.0005
	
	#set map extent
	e<-extent(xmin,xmax,ymin,ymax)
	
	#get map
	g<-gmap(e,type='satellite',lonlat=TRUE)
	
	#create plot
	png(file=plot.name,width=10,height=8,units='in',res=300)
	plot(g,interpolate=TRUE)
	
	#extract only interactions of the specified type
	ints <- ints[which(ints$type==event.type),]
	
	#plot locations in which interactions occurred
	print('plotting interactions...')
	print(nrow(ints))
	for(i in 1:nrow(ints)){
		int<-ints[i,]
		if(int$strength.mult > strength.thresh){
			if(int$disp.mult > disp.thresh){
				t<-int$t2
				i<-int$leader
				j<-int$follower
				xi<-lons[i,t]
				yi<-lats[i,t]
				xj<-lons[j,t]
				yj<-lons[j,t]
				points(xi,yi,pch=20,cex=0.01,col='red')
				points(xj,yj,pch=20,cex=0.01,col='red')
				
			}
		}
		
	}
	
	dev.off()
		
	
	
	
}

