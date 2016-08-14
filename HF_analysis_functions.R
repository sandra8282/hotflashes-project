
#FUNCTION getmode: Find the mode 
    #Parameters: v - a vector of data
    #Returns: the mode of the given vector
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }

#FUNCTION getsmooth: Smooth each woman's data Using K-smooth Filter

  #Parameters: 
    # time - a vector that contains time in seconds for each observation
    # measure - a vector that contains the skin conductance observations
  #Returns: 
    # a data frame with columns "timeID","observed", "fitted.normal". 
    # timeID and observed - the time and measure vector we gave to the function  
    # fitted.normal is the fitted values of the ksmooth (gaussian filter)
    
    get.smooth<-function(time, measure) {
      filterwoman <- na.omit(data.frame(time,measure))
      filterwoman[,3] <- ksmooth(filterwoman[,1], filterwoman[,2], 
                                 kernel = "normal", bandwidth = 120)[[2]]
      filterwoman[,4] <- rep(median(filterwoman[,3]), nrow(filterwoman))
      colnames(filterwoman)<-c("timeID","observed", "fitted.normal", "median")# "ma.flat"
      # Uncomment following if we want to use MA filters
      # filterwoman[,5] <- filter(filterwoman[,2], sides=2, filter=rep(1/10800,10800), 
      #                           method = "convolution")
      
      
      filterwoman<-na.omit(filterwoman)
      return(filterwoman)
    }

#FUNCTION reduce: Reduces one column of data with 1 obs per second to 1 obs per 10 secs 
                # (the average of 10 seconds of raw data)

    #Parameters: 
      # row.from.data - a vector indicating one row of raw data (one woman)
      # index - a vector indicating sequences of every ten until end of data (1, 11, 21,...)
    
    #Returns: temp - a vector of data with averages from 10 seconds of raw data

      reduction<-function(row.from.data, index) {
        temp<-rep(NA, length(index))
        for (i in 1:(length(temp)-1)) {
          temp[i]<-mean(row.from.data[index[i]:index[i+1]])
        }
        return(temp)
      }

#FUNCTION find.peaks: Find peaks in kmooth curve

    #Parameters: 
      # filterwoman.list.element is a dataframe with columns 
      # "timeID","observed","fitted.ma", "ma.flat", "fitted.normal"
      # observed - the measure vector we gave to the function 
      # fitted.normal - the fitted values of the gaussian filter (ksmooth)
    
    #Returns: either a data frame (if more than 1 peak) or a vector (if only 1 peak) 
    # columns in data set or vector are: "peak max umho", "peak index", 
    # "begin index", "end index", "peak begin time", "peak max time", "peak end time"
    find.peaks<- function(filterwoman.list.element) {   
        library(pracma)
        filterwoman<-filterwoman.list.element
        p<-findpeaks(filterwoman$fitted.normal, minpeakheight=1.1, nups = 6, ndowns = 6,
                     minpeakdistance = 90, sortstr=FALSE)
        
        if (class(p)=="matrix" && nrow(p)>1) 
        {
          p <- p[order(p[,2]),]
        } else {p <- NULL}
        
        if (class(p)=="matrix") {
          peakpts<-as.data.frame(p)
          peakpts[,5]<-filterwoman[peakpts[,3], 1] #peak begin time
          peakpts[,6]<-filterwoman[peakpts[,2], 1] #peak max time
          peakpts[,7]<-filterwoman[peakpts[,4], 1] #peak end time
          peakpts[,8]<-filterwoman[peakpts[,3], 3] #peak.begin.umho
          peakpts[,9]<-filterwoman[peakpts[,4], 3] #peak.end.umho
          colnames(peakpts)<-c("peak.umho", "peak.index", "begin.index",
                               "end.index", "peak.begin.time", "peak.max.time",
                               "peak.end.time", "peak.begin.umho", "peak.end.umho")
          peakpts[,10]<-filterwoman[peakpts[,3], 4] #median of fitted umho 
          peakpts[,11]<-peakpts$peak.umho - peakpts[,10]
          notpeak<- which(peakpts$V11<0.5)
          if (length(notpeak)>0) {peakpts<-peakpts[-notpeak,]}
          result<-cbind(peakpts[,1], peakpts[,5:9])
          colnames(result)<-c(colnames(peakpts)[1], colnames(peakpts)[5:9])
          return(result)} else {return(p)}
      }  