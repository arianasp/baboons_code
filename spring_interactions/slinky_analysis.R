#Slinky analysis

library(psych)
library(fields)
library(heR.Misc)

#-----------------------------FUNCTIONS-----------------------------------------------

#EVENT.COUNT.MAT
#make a matrix M, where M(i,j)=number of time i lead j in events of type event.type
#if event.type=='pull.frac', then M(i,j) = # of times i pulled j / (# of times i pulled j + # of times i anchored j)
#Inputs:
#	ints: data frame containing all of the interaction data
#	day: which day to compute the event matrix for
#	event.type: 'pull','push','anchor', or 'repel'
#	min.strength: events with strengths below this threshold will not be counted
#	min.disp: events with disparities below this threshold will not be counted
#Outputs:
#	event.mat: matrix of results (# of events where i lead j during that day)
event.count.mat<-function(ints, day, event.type, min.strength,min.disp){
	
	N <- length(union(unique(ints$leader),unique(ints$follower)))
	
	#need to treat pull.frac differently than regular events - use recursion
	if(event.type == 'pull.frac'){
		pull.mat <- event.count.mat(ints,day,'pull',min.strength,min.disp)
		anchor.mat <- event.count.mat(ints,day,'anchor',min.strength,min.disp)
		event.mat = pull.mat / (pull.mat + t(anchor.mat))
		event.mat[which(pull.mat)]
	}
	
	else{
		#pull out events on the correct day and that exceed the threholds
		data<-ints[which(ints$day==day & ints$event.type==event.type & ints$strength >= min.strength & ints$disp >= min.disp),]
	
		#put them into a matrix where event.mat(i,j) = # of time where i lead j
		event.mat<-matrix(0,nrow=N,ncol=N)
	
		for(i in 1:dim(data)[1]){
			event.mat[data$leader[i],data$follower[i]] <- event.mat[data$leader[i],data$follower[i]] + 1
		}
	}
	
	event.mat[is.nan(event.mat)]<-NA
	
	invisible(event.mat)
	
}

#DIREC.MAT
#convert an event matrix into a directionality matrix
#direc.mat(i,j) = (event.mat[i,j]-event.mat[j,i])/(event.mat[i,j]+event.mat[j,i])
#Inputs:
#	event.mat: matrix giving number of events where i lead j
#	min.tot.events: minimum total number of events needed in order to compute the directionality (otherwise that cell will get a NA)
#Outputs:
#	direc.mat: directionality matrix
direc.mat<-function(event.mat,min.tot.events){
	
	N<-dim(event.mat)[1]
	direc.mat<-matrix(NA,nrow=N,ncol=N)
	for(i in seq(1,N-1,1)){
		for(j in seq(i+1,N,1)){
			if(!is.na(event.mat[i,j]) & !is.na(event.mat[j,i])){
				if(event.mat[i,j]+event.mat[j,i]>=min.tot.events){
					direc.mat[i,j]<-(event.mat[i,j] - event.mat[j,i]) / (event.mat[i,j] + event.mat[j,i])
					direc.mat[j,i]<--direc.mat[i,j]
				}
			}
		}
	}
	
	direc.mat[which(is.nan(direc.mat))]<-NA
	
	invisible(direc.mat)
	
}


#RANK.MAT
#rank the individuals in a matrix based on a function of the rows or columns (can be 'sum' or 'mean')
#Inputs:
#	mat: an NxN matrix
#	dimension: 'row' or 'col' (whether to rank based on row or column sums)
#	direction: 'descend' or 'ascend' (whether to rank them from smallest to largest (ascend) or largest to smallest (descend))
#Outputs:
#	ranks: rank of each individual (where 1 = most leading and N = most following)
rank.mat <- function(mat,func='mean',dimension='col',direction='descend'){
	
	#number of rows/cols
	N <- dim(mat)[1]
	
	#get values to rank based on (mean or sum of rows or columns)
	if(dimension == 'row'){
		if(func == 'mean'){
			vals <- rowMeans(mat,na.rm=T)
		}
		if(func == 'sum'){
			vals <- rowSums(mat,na.rm=T)
		}
	}
	if(dimension == 'col'){
		if(func == 'mean'){
			vals <- colMeans(mat,na.rm=T)
		}
		if(func == 'sum'){
			vals <- colSums(mat,na.rm=T)
		}
	}
	
	vals[which(is.nan(vals))]<-NA
	
	#put into a data frame with the individual ids
	dat <- data.frame(id=seq(1,N,1),vals=vals)
	
	
	#sort the data frame
	if(direction == 'ascend'){
		dat <- dat[order(-dat$vals),]
	}
	if(direction == 'descend'){
		dat <- dat[order(dat$vals),]
	}
	
	
	#get the ranks from the order
	dat$ranks <- 1:N
	
	#reorder based on the original ids
	dat<-dat[order(dat$id),]
	
	if(any(is.na(dat$vals))){
		dat[which(is.na(dat$vals)),]$ranks <- NA
	}
	
	ranks <- dat$ranks
	invisible(ranks)
	
}

#RE.RANK.SUBSET
#re-rank a group of individuals that is a subset of the originally ranked set
#(so your new rank will be out of only the individuals present)
#Inputs:
#	ranks: vector or list containing the original ranks of each individual
#	inds.to.include: vector or list of the individuals to include in the new ranking
re.rank.subset <- function(ranks, inds.to.include=1:length(ranks)){
	N <- length(ranks)
	dat <- data.frame(inds = 1:N, orig.ranks = ranks)
	dat <- dat[inds.to.include,]
	
	dat <- dat[order(dat$orig.ranks),]
	dat$new.ranks <- 1:length(inds.to.include)
	
	dat <- dat[order(dat$inds),]
	
	new.ranks <- dat$new.ranks
	
	invisible(new.ranks)
	
	
}

#RANKS.ACROSS.DAYS
#compute the intraclass correlation coefficient for the rankings across days for leadership hierarchies
ranks.across.days <- function(ints, inds.to.include=c(1:7,9:11,15,18,19,21,22,26), days=seq(1,14,1), event.type='pull', min.strength=0.1,min.disp=0.1,min.tot.events=0){
	
	# number of individuals
	N <- length(union(unique(ints$leader),unique(ints$follower)))
	
	#number of days
	n.days <- length(days)
	
	#array to hold all matrices for each day
	event.mats <- array(NA,dim=c(N,N,n.days))
	direc.mats <- array(NA,dim=c(N,N,n.days))
	direc.mats.subset <- array(NA,dim=c(length(inds.to.include),length(inds.to.include),n.days))
	rank.vecs <- array(NA,dim=c(N,n.days))
	rank.vecs.subset <- array(NA,dim=c(length(inds.to.include),n.days))
	
	#get data matrix and ranks over all days for all individuals and for a subset
	event.mat.tot <- event.count.mat(ints,day=days,event.type,min.strength,min.disp)
	direc.mat.tot <- direc.mat(event.mat.tot,min.tot.events)
	ranks.tot <- rank.mat(direc.mat.tot,func='mean',dimension='col',direction='descend')
	
	direcs.tot.subset <- direc.mat.tot[inds.to.include,inds.to.include]
	ranks.tot.subset <- rank.mat(direcs.tot.subset,func='mean',dimension='col',direction='descend')
	
	#get data matrices and ranks for each day
	for(d in 1:length(days)){
		day <- days[d]
		events <- event.count.mat(ints,day,event.type,min.strength,min.disp)
		direcs <- direc.mat(events,min.tot.events)
		ranks <- rank.mat(direcs,func='mean',dimension='col',direction='descend')

		#include only a subset of individuals in the rankings and icc computation
		direcs.subset <- direcs[inds.to.include,inds.to.include]
		ranks.subset <- rank.mat(direcs.subset,func='mean',dimension='col',direction='descend')
		
		#store data
		event.mats[,,d] <- events
		direc.mats[,,d] <- direcs
		rank.vecs[,d] <- ranks
		rank.vecs.subset[,d]<-ranks.subset
		direc.mats.subset[,,d] <- direcs.subset
	}
	
	#compute icc on ranks
	if(any(is.na(rank.vecs.subset))){
		icc <- NA
		icc.lower.bound <- NA
		icc.upper.bound <- NA
	}
	else{
		icc.data <- ICC(rank.vecs.subset)
		icc<- icc.data$results$ICC[6]
		icc.lower.bound <- icc.data$results$lower[6]
		icc.upper.bound <- icc.data$results$upper[6]
	}
	
	#compute icc on directionality coefficients
	direcs.by.day <- matrix(NA,nrow=(length(inds.to.include)^2-length(inds.to.include))/2,ncol=n.days)
	for(d in 1:n.days){
		curr.direcs <- direc.mats.subset[,,d]
		direcs.by.day[,d]<-curr.direcs[upper.tri(curr.direcs,diag=FALSE)]
	}
	
	if(any(is.na(direcs.by.day))){
		direcs.icc <- NA
		direcs.icc.lower.bound <- NA
		direcs.icc.upper.bound <- NA
	}
	else{
		icc.data <- ICC(direcs.by.day)
		direcs.icc <- icc.data$results$ICC[6]
		direcs.icc.lower.bound <- icc.data$results$lower[6]
		direcs.icc.upper.bound <- icc.data$results$upper[6]
	}
	
	out<-list()
	out$event.mat.tot <- event.mat.tot
	out$direc.mat.tot <- direc.mat.tot
	out$ranks.tot <- ranks.tot
	out$ranks.tot.subset <- ranks.tot.subset
	out$event.mats <- event.mats
	out$direc.mats <- direc.mats
	out$rank.vecs <- rank.vecs
	out$inds.included.in.subset <- inds.to.include
	out$rank.vecs.subset <- rank.vecs.subset
	out$icc <- icc
	out$icc.lower.bound <- icc.lower.bound
	out$icc.upper.bound <- icc.upper.bound
	out$direcs.icc <- direcs.icc
	out$direcs.icc.lower.bound <- direcs.icc.lower.bound
	out$direcs.icc.upper.bound <- direcs.icc.upper.bound
	
	invisible(out)
	
}

#RANK.CONSISTENCY.PLOT
rank.consistency.plot<-function(ranks.tot,rank.vecs,ids){
	quartz()
	n.days <- dim(rank.vecs)[2]
	n.inds <- dim(rank.vecs)[1]
	colors <- ids$color
	plot(NA,xlim=c(-1,n.days),ylim=c(n.inds,1),ylab='rank',xlab='day',xaxt='n',yaxt='n')
	par(las=1)
	axis(1,at=seq(-1,n.days,1),labels=c('','overall',seq(1,n.days,1)),cex.axis=0.75)
	axis(2,at=seq(1,n.inds,1),labels=seq(1,n.inds,1),cex.axis=0.75)
	for(i in 1:n.inds){
		lines(0:n.days,c(ranks.tot[i],rank.vecs[i,]),col=colors[i])
		points(0,ranks.tot[i],pch=19,col=colors[i])
		text(-1,ranks.tot[i],col=colors[i],labels=ids$Collar[i])
	}
	
}

#AGE.SEX.COLORS
age.sex.colors <- function(ids){
	N <- nrow(ids)
	ids$color <- rgb(1,0,0)
	ids[which(ids$Sex=='M' & ids$Age == 'A'),]$color <- rgb(0,0,1)
	ids[which(ids$Sex=='M' & ids$Age == 'SA'),]$color <- rgb(102/255,204/255,1)
	ids[which(ids$Sex=='F' & ids$Age == 'A'),]$color <- rgb(1,0,0)
	ids[which(ids$Sex=='F' & ids$Age == 'SA'),]$color <- rgb(1,204/255,102/255)
	ids[which(ids$Age == 'J'),]$color <- rgb(102/255,102/255,102/255)
	
	invisible(ids)
	
}

#CONSISTENCY.ACROSS.DAYS.OVER.PARAM.RANGE
#consistency across days for various combinations of parameters (strength and disparities)
#computes the ICC across the first 14 days (for the individuals that are present over all those days) for various combos of parameters
consistency.across.days.over.param.range <- function(ints,event.type='pull',min.strengths=seq(0,0.9,.1),min.disps=seq(0,0.9,.1)){
	n.strengths <- length(min.strengths)
	n.disps <- length(min.disps)
	rank.iccs <- matrix(NA,nrow=n.strengths,ncol=n.disps)
	direc.iccs <- matrix(NA,nrow=n.strengths,ncol=n.disps)
	rank.icc.lower.bounds <- matrix(NA,nrow=n.strengths,ncol=n.disps)
	rank.icc.upper.bounds <- matrix(NA,nrow=n.strengths,ncol=n.disps)
	direc.icc.lower.bounds <- matrix(NA,nrow=n.strengths,ncol=n.disps)
	direc.icc.upper.bounds <- matrix(NA,nrow=n.strengths,ncol=n.disps)
	for(i in 1:n.strengths){
		for(j in 1:n.disps){
			print(paste('min strength = ',min.strengths[i], ', min disp = ',min.disps[j],sep=''))
			dat <- ranks.across.days(ints,event.type=event.type,min.strength=min.strengths[i],min.disp=min.disps[j])
			rank.iccs[i,j]<-dat$icc
			rank.icc.lower.bounds[i,j] <-dat$icc.lower.bound
			rank.icc.upper.bounds[i,j] <-dat$icc.upper.bound
			direc.iccs[i,j]<-dat$direcs.icc
			direc.icc.lower.bounds[i,j]<-dat$direcs.icc.lower.bound
			direc.icc.upper.bounds[i,j]<-dat$direcs.icc.upper.bound
			
		}
	}
	
	out<-list()
	out$min.strengths <- min.strengths
	out$min.disps <- min.disps
	out$rank.iccs <- rank.iccs
	out$rank.icc.lower.bounds <- rank.icc.lower.bounds
	out$rank.icc.upper.bounds <- rank.icc.upper.bounds
	out$direc.iccs <- direc.iccs
	out$direc.icc.lower.bounds <- direc.icc.lower.bounds
	out$direc.icc.upper.bounds <- direc.icc.upper.bounds
	
	
	invisible(out)
	
}

#ICC.RAND.RANKS
icc.rand.ranks <- function(n.inds,n.days,n.reps){
	
	iccs <- rep(NA,n.reps)
	for(r in 1:n.reps){
		ranks <- matrix(NA,n.inds,n.days)
		for(i in 1:n.days){
			ranks[,i] <- sample(seq(1,n.inds,1))
		}
	
		icc.data <- ICC(ranks)
		icc<- icc.data$results$ICC[6]
		iccs[r] <- icc
	}
	invisible(iccs)
	
}

#CONSISTENCY.ACROSS.PARAMS
#computes the correlation the overall ranks defined by every combination of strengths and disparities and a given base strength/disparity combo
consistency.across.params <- function(ints,event.type='pull',inds.to.include=c(1:7,9:11,15,18,19,21,22,26),day.range=c(1,14),strengths=seq(0,.9,.1),disps=seq(0,0.9,.1),strength0=0,disp0=0,min.tot.events = 0){
	
	#get number of individuals and number of parameter combos
	n.inds <- length(inds.to.include)
	n.strengths <- length(strengths)
	n.disps <- length(disps)
	
	#get ranks to compare everything to (based on base strength/disparity combo)
	events0 <- event.count.mat(ints,day.range,event.type,strength0,disp0)
	direcs0 <- direc.mat(events0,min.tot.events)
	direcs0 <- direcs0[inds.to.include,inds.to.include]
	ranks0 <- rank.mat(direcs0)
	n.events0 <- sum(events0,na.rm=T)
	
	#array to store ranks for all parameter combos
	ranks.across.params <- array(NA,dim=c(n.strengths,n.disps,n.inds))
	
	#array to store correlation values between ranks for base strength/disparity combo and ranks for each other parameter combo
	corrs.across.params <- array(NA,dim=c(n.strengths,n.disps))
	
	n.events.across.params <- array(NA,dim=c(n.strengths,n.disps))
	
	for(i in 1:n.strengths){
		for(j in 1:n.disps){
			print(paste('min strength = ',strengths[i], ', min disp = ',disps[j],sep=''))
			strength <- strengths[i]
			disp <- disps[j]
			events <- event.count.mat(ints,day.range,event.type,strength,disp)
			direcs <- direc.mat(events,min.tot.events)
			direcs <- direcs[inds.to.include,inds.to.include]
			ranks <- rank.mat(direcs)
			ranks.across.params[i,j,] <- ranks
			corrs.across.params[i,j] <- cor(ranks,ranks0,method='spearman')
			n.events.across.params[i,j] <- sum(events,na.rm=T)
		}
	}
	
	out<-list()
	out$strength0 <- strength0
	out$disp0 <- disp0
	out$strengths <- strengths
	out$disps <- disps
	out$inds.included <- inds.to.include
	out$ranks.across.params <- ranks.across.params
	out$corrs.across.params <- corrs.across.params
	out$n.events.across.params <- n.events.across.params
	
	invisible(out)
	
	
}


#ANGLE.DIFFS
#for pull events: computes and plots the angle between the displacement vector of the leader during t1 --> t2
#and the displacement vector of the follower during t2 --> t3
#for anchor events: computes and plots the angle between the displacement vector of the follower (anchored) during t1 --> t2
#and the displacement vector of the follower (anchored) during t2 --> t3
angle.diffs <- function(ints,event.type,leader,follower,day){
	data <- ints[which(ints$event.type == event.type & ints$leader == leader & ints$follower == follower & ints$day==day ),]
	
	if(event.type == 'pull'){
		data$l1 <- sqrt(data$lead.dx.12^2 + data$lead.dy.12^2)
		data$l2 <- sqrt(data$foll.dx.23^2 + data$foll.dy.23^2)
		data$norm.x1 <- data$lead.dx.12 / data$l1
		data$norm.x2 <- data$foll.dx.23 / data$l2
		data$norm.y1 <- data$lead.dy.12 / data$l1
		data$norm.y2 <- data$lead.dy.23 / data$l2
	}
	
	if(event.type == 'anchor'){
		data$l1 <- sqrt(data$foll.dx.12^2 + data$foll.dy.12^2)
		data$l2 <- sqrt(data$foll.dx.23^2 + data$foll.dy.23^2)
		data$norm.x1 <- data$foll.dx.12 / data$l1
		data$norm.x2 <- data$foll.dx.23 / data$l2
		data$norm.y1 <- data$foll.dy.12 / data$l1
		data$norm.y2 <- data$foll.dy.23 / data$l2
		
	}
	
	data$dx <- data$norm.x2 - data$norm.x1
	data$dy <- data$norm.y2 - data$norm.y1
	
	data$dtheta <- atan2(data$dy,data$dx)
	
	rose.diag(data$dtheta,bins=20)
	
	
	
	
}



























#--------------------------------------MAIN-------------------------------------------

#parameters:
datadir		 	<- '/Users/arianasp/Desktop/Baboons/output/push_pull/interactions_data'
idsdir		 	<- '/Users/arianasp/Desktop/Baboons/data/csv_raw'
noise.thresh 	<- 5
min.strength 	<- 0.1
min.disp 	 	<- 0.1
day.range 	 	<- 1:14
inds.to.include <- c(1:7,9:11,15,18,19,21,22,26)

#load data
filename <- paste(datadir,'/dyad_interactions_noise_thresh_',noise.thresh,'.csv',sep='')
ints <- read.table(file=filename,header=TRUE,sep=',')
ids <- read.table(file=paste(idsdir,'/IDs.csv',sep=''),header=TRUE,sep=',')
ids <- age.sex.colors(ids)

#subset data
ids.subset <- ids[inds.to.include,]

#plotting an image plot of repeatability across parameters:
#image.plot(x=pull.iccs$min.strengths,y=pull.iccs$min.disps,z=pull.iccs$iccs,zlim=c(0.9,1),col=colorpanel(256,'black','red'),xlab='strength threshod',ylab='disparity threshold',asp=1,main='repeatability of pull hierarchy ranks across days (ICC)')


#load interactions