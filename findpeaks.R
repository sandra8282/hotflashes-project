check.if.not.peak<-function(notpeak, peakpts) {
  if (length(notpeak)>0) {peakpts<-peakpts[-notpeak,]}
  return(peakpts)}

order.pts<-function(p) {
  if (class(p)=="matrix" && nrow(p)>1) {p <- p[order(p[,2]),]
  } else {p <- NULL} 
  return(p)}

# FUNCTION:find.best.in.group
# given one of the peak.groups we define the peak as the point in the group 
# with max.umho and return its umho and time
# we say that peak starts at the first point in the group 
# and return its begin umho and time
# we say that peak ends at the last point in the group 
# and return its end umho and time
# PARAMETERS: one element of the list peak.groups
# RETURNS: a vector with the information for the best peak (begin, max, end info)

          find.best.in.group<-function(temp) {
            
            best.peak<-rep(NA, (ncol(temp)))
            
            highest.umho<-which(temp$peak.umho==max(temp$peak.umho))
            end.time<-which(temp$peak.end.time==max(temp$peak.end.time))
            
            best.peak[1]<-temp[highest.umho, 1]             #peak.index
            best.peak[2]<-temp$begin.index[1]                #begin.index
            best.peak[3]<-temp$end.index[highest.umho]      #end.index
            best.peak[4]<- max(temp$peak.umho)    #peak.umho
            best.peak[5]<- temp[highest.umho, 5]  #peak.max.time
            best.peak[6]<- temp$peak.begin.umho[1]  #peak.begin.umho
            best.peak[7]<- temp$peak.begin.time[1] #peak.begin.time
            best.peak[8]<- temp$peak.end.umho[end.time] #peak.end.umho
            best.peak[9]<- temp$peak.end.time[end.time] #peak.end.time
            best.peak[10] <-best.peak[4]-best.peak[6]
            best.peak[11] <- best.peak[4]-best.peak[8]
            best.peak<-as.data.frame(t(best.peak))
            colnames(best.peak)<-colnames(temp)
            return(best.peak[1:11])
          }

#FUNCTION find.peaks: Find peaks in kmooth curve

#Parameters: 
# filterwoman.list.element is a dataframe with columns 
# "timeID","observed", "fitted.normal","MA.median", "Diff.Zero.Median"
# observed - the measure vector we gave to the function 
# fitted.normal is the fitted values of the ksmooth (gaussian filter)
# MA.median is the median from a moving average smooth version of the data
# Diff.Zero.Median - median for the points where the first difference is zero

#Returns: either a data frame (if more than 1 peak) or a vector (if only 1 peak) 
# columns in data set or vector are: "peak max umho", "peak index", 
# "begin index", "end index", "peak begin time", "peak max time", "peak end time"

find.peaks<- function(filterwoman) {
  #filterwoman<-BM.filterwoman.list[[1]]
  
  library(pracma)
  largest.peak.dist<-max(filterwoman$fitted.normal)-filterwoman$Diff.Zero.Median[1]
  p<-findpeaks(filterwoman$fitted.normal, minpeakheight=largest.peak.dist/4, 
               sortstr=FALSE)
  p<-order.pts(p)
  if (class(p)=="matrix") {
    peakpts<-as.data.frame(cbind(p[,2:4], p[,1]))
    colnames(peakpts)<-c("peak.index", "begin.index",
                         "end.index", "peak.umho")
    peakpts$peak.max.time<-filterwoman[peakpts[,1], 1] 
    peakpts$peak.begin.umho<-filterwoman[peakpts[,2], 3] 
    peakpts$peak.begin.time<-filterwoman[peakpts[,2], 1] 
    peakpts$peak.end.umho<-filterwoman[peakpts[,3], 3] 
    peakpts$peak.end.time<-filterwoman[peakpts[,3], 1] 
    peakpts$begin2peak.umho<-peakpts$peak.umho-peakpts$peak.begin.umho
    peakpts$peak2end.umho<-peakpts$peak.umho-peakpts$peak.end.umho
    
    #Check for very small differences in 
    # peak max umho - peak begin umho
    notpeak<- which(peakpts$begin2peak.umho<0.1)  #identify very small differences
    peakpts<-check.if.not.peak(notpeak, peakpts) #remove pts w/ very small diff
    
    #plot to see what we have
      # plot(filterwoman[[1]]/3600, filterwoman[[3]], type = "l")
      # points(peakpts$peak.max.time/3600, peakpts$peak.umho, pch=20, col="red")
      # points(peakpts$peak.begin.time/3600, peakpts$peak.begin.umho, pch=20, col="green")
    # points(peakpts$peak.end.time/3600, peakpts$peak.end.umho, pch=20, col="blue")
    # 
    #Check difference between the begin2peak distances 
    # peakpts$end.diff<-c(NA,diff(peakpts$peak2end.umho))
    # peakpts$begins.diff<-c(NA, diff(peakpts$begin2peak.umho))
    
    if (nrow(peakpts)==0) { #no peaks found return NA
      semifinal.best.peaks<-NA
    } else {
      if (nrow(peakpts)==1) { #only one peak exists so return that one
        if (peakpts$peak.umho>=1) {semifinal.best.peaks<-peakpts} else {semifinal.best.peaks<-NA}
      } else {
        
          # First find peaks that are isolated and save those
          leads.from.begin<-which(peakpts$begin2peak.umho>1)
          leads.from.ends<-which(peakpts$peak2end.umho>1)
          if (length(leads.from.begin)>=length(leads.from.ends)) {
            same.lead <- sapply(leads.from.ends, function(x) which(leads.from.begin==x))
            same.lead <- unlist(same.lead)
            lone.peaks.index<-leads.from.begin[same.lead]
          } else {
            same.lead <- sapply(leads.from.begin, function(x) which(leads.from.ends==x))
            same.lead <- unlist(same.lead)
            lone.peaks.index<-leads.from.ends[same.lead]
          }
          lone.peaks<-peakpts[lone.peaks.index,]
          
          # Second from non-isolated peaks find a grouping
          if (length(lone.peaks.index)>0){ #Take out lone peaks if they exist
            newpeakpts<-peakpts[-lone.peaks.index,]
          } else {newpeakpts<-peakpts} 
        
          #If there are any peaks look at:
          # 1) which diff of begin2peak levels are >1
          # 2) which consecutive peak points 
          if (nrow(newpeakpts)>1){ 
            newpeakpts$beg.diff<-c(NA, diff(newpeakpts$begin2peak.umho))
            newpeakpts$end.diff<-c(NA, diff(newpeakpts$peak2end.umho))
            group.leads.frombegin <-c(1, which(abs(newpeakpts$beg.diff)>=0.5))
            group.leads.fromend <- c(1, which(abs(newpeakpts$end.diff)>=0.5))
            group.leads <- unique(c(group.leads.frombegin, group.leads.fromend))
                if (length(group.leads)>=2) {
                  peak.groups<-list()
                  if (length(group.leads)==2) {
                    peak.groups[[1]]<-newpeakpts[1:(group.leads[2]-1),]
                    peak.groups[[2]]<-newpeakpts[group.leads[2]:nrow(newpeakpts),]
                  } else {
                    for (i in 1:(length(group.leads)-1)) {
                      peak.groups[[i]]<-newpeakpts[group.leads[i]:(group.leads[i+1]-1),]
                    }
                  }
                  best.peak.list<-lapply(peak.groups, function(x) find.best.in.group(x))
                  best.peakpts<-rep(NA, nrow(best.peak.list[[1]]))
                  for (d in 1:length(best.peak.list)) {
                    best.peakpts<-rbind(best.peakpts, best.peak.list[[d]])
                  }
                  best.peakpts<-as.data.frame(best.peakpts[-1,])
                  colnames(best.peakpts)<-colnames(peakpts)[1:ncol(best.peakpts)]
                } else {
                  best.peakpts<-find.best.in.group(newpeakpts)}
                  semifinal.best.peaks<-rbind(best.peakpts,lone.peaks)
          } else {semifinal.best.peaks<-lone.peaks}
          semifinal.best.peaks<-semifinal.best.peaks[order(semifinal.best.peaks$peak.max.time),]
          semifinal.best.peaks$time.diff <- c(NA, diff(semifinal.best.peaks$peak.max.time)/60)
          # peaks.too.close <- which(na.omit(semifinal.best.peaks$time.diff)<=20)
          # peaks.not.close <-semifinal.best.peaks[-c(peaks.too.close, peaks.too.close-1),]
          # if (length(peaks.too.close)>1) {
          #     best.close.peaks<-list()
          #     c.count=2
          #     while (c.count<=length(peaks.too.close)) {
          #       close.peaks.group<-rbind(semifinal.best.peaks[(peaks.too.close[c.count])-1,],
          #                                           semifinal.best.peaks[peaks.too.close[c.count],])
          #       best.close.peaks[[c.count]]<-find.best.in.group(close.peaks.group)
          #       c.count<-c.count+1
          #     }
          #     for (d2 in 1:length(best.close.peaks)) {
          #       final.best.peaks<-rbind(peaks.not.close[,1:11], best.close.peaks[[d2]])
          #     }
          #  } else {final.best.peaks<-semifinal.best.peaks[,1:11]}
        }
      }
         
    return(semifinal.best.peaks)} else {return(p)}
}



    #check what we have
    # plot(filterwoman[[1]]/3600, filterwoman[[3]], type = "l")
    # points(semifinal.best.peaks$peak.max.time/3600, semifinal.best.peaks$peak.umho, pch=20, col="red")
    # points(final.best.peaks$peak.begin.time/3600, final.best.peaks$peak.begin.umho, pch=20, col="green")
    # points(final.best.peaks$peak.end.time/3600, final.best.peaks$peak.end.umho, pch=20, col="blue")
    
    