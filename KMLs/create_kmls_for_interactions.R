dir <- '/Users/arianasp/Desktop/Baboons/google_earth/slinky_interactions/noise_thresh_5'

ids <- read.csv("/Users/arianasp/Desktop/Baboons/data/csv_raw/IDs.csv")

ints <- as.data.frame(read.csv(paste(dir, '/interactions_list.csv',sep=''),header=TRUE))

system(paste('mkdir ', dir ,'/kmls',sep=''))

for(int_num in 1:dim(ints)[1]){

# READ IN DATA
if(file.exists(paste(dir,'/latlons/',ints$FILE[int_num],'.csv',sep=''))){
mydata <- read.csv(paste(dir,'/latlons/',ints$FILE[int_num],'.csv',sep=''))
colours <- c('baboon_red.png','baboon_blue.png','baboon_black.png')

# START WRITING
sink(paste(dir,'/kmls/',ints$FILE[int_num],'.kml',sep=''))

# start output
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
cat("<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\">\n")
cat("<Document>\n")


for (i in c(1:nrow(ids))[-c(ints$LEADER[int_num],ints$FOLLOWER[int_num])]) {
	cat(paste("<Style id=\"track-",ids$Collar[i],"\"><IconStyle><scale>0.5</scale><Icon><href>baboon_black.png</href></Icon></IconStyle><LineStyle><color>ff000000</color><colorMode>normal</colorMode></LineStyle></Style>\n",sep=""))
}
cat(paste("<Style id=\"track-",ids$Collar[ints$FOLLOWER[int_num]],"\"><IconStyle><scale>0.5</scale><Icon><href>baboon_red.png</href></Icon></IconStyle><LineStyle><color>ff0000ff</color><colorMode>normal</colorMode></LineStyle></Style>\n",sep=""))
cat(paste("<Style id=\"track-",ids$Collar[ints$LEADER[int_num]],"\"><IconStyle><scale>0.5</scale><Icon><href>baboon_blue.png</href></Icon></IconStyle><LineStyle><color>ffff0000</color><colorMode>normal</colorMode></LineStyle></Style>\n",sep=""))



hour <- as.numeric(substr(mydata$TIME,12,13))
hours <- unique(as.numeric(substr(mydata$TIME,12,13)))
# FOR EACH HOUR
for (h in hours) {
	cat("<Folder>\n")
	cat(paste("<name>",substr(mydata$TIME[1],1,10)," ",sprintf("%02d", h),":00:00</name>\n",sep=""))

	inds <- unique(mydata$ID[which(hour==h)])
	# FOR EACH INDIVIDUAL PRESENT
	for (i in inds) {
		cat("<Folder>\n")
		cat(paste("<name>",ids$Collar[i],"</name>\n",sep=""))
		cat("<Placemark>")
		cat(paste("<styleUrl>#track-",ids$Collar[i],"</styleUrl>\n",sep=""))
		cat("<visibility>0</visibility>\n")
		cat("<gx:Track>\n")
		cat("<gx:altitudeMode>relativeToGround</gx:altitudeMode>\n")
		
		times <- which(hour==h & mydata$ID == i)
		# FOR EACH TIME
		#for (t in times) {
			#cat(paste("<when>",mydata$TIME[t],"</when>\n",sep=""))
			cat(sprintf("<when>%s</when>\n",mydata$TIME[times]),sep="")
		#}
		# FOR EACH TIME
		#for (t in times) {
			#cat(paste("<gx:coord>",mydata$LONLAT[t]," 0</gx:coord>\n",sep=""))
			cat(sprintf("<gx:coord>%s 0</gx:coord>\n",mydata$LONLAT[times]),sep="")
		#}
		cat("</gx:Track>\n")
		cat("</Placemark>\n")
		cat("</Folder>\n")
	}
	cat("</Folder>\n")
}
cat("</Document>\n")
cat("</kml>\n")
sink()
}
}