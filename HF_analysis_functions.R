
#FUNCTION getmode: Find the mode 
    #Parameters: v - a vector of data
    #Returns: the mode of the given vector, the # of times the mode repeated
    getmode <- function(v) {
      uniqv <- unique(v)
      mode<-round(uniqv[which.max(tabulate(match(v, uniqv)))], digit=1)
      repeats<-as.integer(max(tabulate(match(v, uniqv))))
      return(c(mode, repeats))
    }

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#FUNCTION getsmooth: Smooth each woman's data Using K-smooth Filter
    #and find a moving average median

  #Parameters: 
    # time - a vector that contains time in seconds for each observation
    # measure - a vector that contains the skin conductance observations
  #Returns: 
    # a data frame with columns "timeID","observed", "fitted.normal","MA.median"
    # and "Diff.Zero.Median". 
    # timeID and observed - the time and measure vector we gave to the function  
    # fitted.normal is the fitted values of the ksmooth (gaussian filter)
    # MA.median is the median from a moving average smooth version of the data
    # Diff.Zero.Median - median for the points where the first difference is zero
    
    get.smooth<-function(time, measure) {
      filterwoman <- na.omit(data.frame(time,measure))
      filterwoman[,3] <- ksmooth(filterwoman[,1], filterwoman[,2], 
                                 kernel = "normal", bandwidth = 120)[[2]]
      filterwoman[,4] <- rep(median(na.omit(filterwoman[,3])), nrow(filterwoman))
      d<-diff(filterwoman[,3])
      flat.umho<-filterwoman[which(d==0),3]
      filterwoman[,5]<-rep(median(flat.umho), nrow(filterwoman))
      colnames(filterwoman)<-c("timeID","observed", "fitted.normal", "MA.median", "Diff.Zero.Median")
      return(filterwoman)
    }

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
    
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


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#FUNCTION:Plot data, ksmooth, and peak points given 1) list of data&ksmooth
      #2) list of peak points for each group (each element in both lists is a woman)
      #3) the reduced data frame for all women that was used to get data$kmooth list
#PARAMETERS: filterwoman.list, peak.from.ksmooth, reduced.data
#RETURNS: plots for each woman
      
plot.smooth.and.peaks<-function(filterwoman.list, peak.from.ksmooth, reduced.data) {
  # filterwoman.list<-BM.filterwoman.list
  # peak.from.ksmooth<-BM.peak.from.ksmooth
  # reduced.data<-reduced.BM
  plots<-list()
  library(gridExtra)
  library(Rmisc)
  library(ggplot2)
  for (i in 1:length(filterwoman.list)) {
    if (class(peak.from.ksmooth[[i]])=="data.frame"&&nrow(peak.from.ksmooth[[i]]>1)) {
      plots[[i]]<-
        ggplot(data=filterwoman.list[[i]], aes(x=timeID/3600, y=observed)) +
        theme_bw() + geom_line(colour="grey") +
        geom_line(data=filterwoman.list[[i]], aes(x=timeID/3600,y=fitted.normal),
                  colour="black") +
        # geom_line(data=filterwoman.list[[i]], aes(x=timeID/3600,y=hourly.median),
        #          colour="blue") +      
        # geom_line(data=filterwoman.list[[i]], aes(x=timeID/3600,y=ma.flat),
        #     colour="green") +
        geom_point(data=peak.from.ksmooth[[i]], aes(x = peak.max.time/3600,
                                                    y = peak.umho, colour = "red")) +
        geom_point(data=peak.from.ksmooth[[i]], aes(x = peak.begin.time/3600,
                                             y = peak.begin.umho, colour = "green")) +
        geom_point(data=peak.from.ksmooth[[i]], aes(x = peak.end.time/3600,
                                 y = peak.end.umho, colour = "blue")) +
        labs(x = "", y = "") +
        ggtitle(colnames(reduced.data)[i+1]) +
        theme(legend.position="none")
    } else {
      plots[[i]]<-
        ggplot(data=filterwoman.list[[i]], aes(x=timeID/3600, y=observed)) +
        theme_bw() + geom_line(colour="grey") +
        geom_line(data=filterwoman.list[[i]], aes(x=timeID/3600,y=fitted.normal),
                  colour="black") +
        # geom_line(data=filterwoman.list[[i]], aes(x=timeID/3600,y=hourly.median),
        #           colour="blue") +        
        labs(x = "", y = "") +
        ggtitle(colnames(reduced.data)[i+1]) +
        theme(legend.position="none")
      
    }
  }
  return(plots)
}


#--------------
#BINNED MODES FUNCTION
get.binned.modes<-function(filterwoman){
  max.fw<-max(filterwoman$fitted.normal)
  rounded.down.min<-floor(min(filterwoman$fitted.normal))
  split.umho<-seq(rounded.down.min, max.fw+0.4, by=0.2) 
  split.list<-list()
  for (n in 1:(length(split.umho)-1)) {
    split.list[[n]]<-filterwoman[which(filterwoman$fitted.normal<split.umho[n+1]),]
    split.list[[n]]<-split.list[[n]][which(split.list[[n]]$fitted.normal>split.umho[n]),]
  }
  nrow.vect<-unlist(lapply(split.list, function(x) nrow(x)))
  max.nrow<-max(nrow.vect)
  max.obs<-which(nrow.vect==max.nrow)
  if (length(max.obs)>1) {max.obs<-min(max.obs)}
  best.mode.data=split.list[[max.obs]]
  best.mode<-getmode(best.mode.data$fitted.normal)
  # plot(filterwoman$timeID, filterwoman$fitted.normal)
  # lines(filterwoman$timeID, rep(best.mode[1], length(filterwoman$timeID)), col="blue")
  return(best.mode[1])
}

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

plot.STATS<-function(filterwoman.list, summary.stats, reduced.data) {
  plots.of.stats<-list()
  library(gridExtra)
  library(Rmisc)
  library(ggplot2)
  names<-colnames(reduced.data)
  for (i in 1:length(filterwoman.list)) {
  plots.of.stats[[i]]<-
     plot(filterwoman.list[[i]]$timeID/3600, 
          filterwoman.list[[i]]$fitted.normal, 
          ylab="Umho",xlab="Time (Hrs)", type = "l", main = names[i+1])
     lines(filterwoman.list[[i]]$timeID/3600, 
           rep(summary.stats[i,1], length(filterwoman.list[[i]]$timeID)), col="blue")
     lines(filterwoman.list[[i]]$timeID/3600, 
           rep(summary.stats[i,2], length(filterwoman.list[[i]]$timeID)), col="red")
     lines(filterwoman.list[[i]]$timeID/3600, 
           rep(summary.stats[i,3], length(filterwoman.list[[i]]$timeID)), col="green")
     lines(filterwoman.list[[i]]$timeID/3600, 
           rep(summary.stats[i,4], length(filterwoman.list[[i]]$timeID)), col="purple")
  }
  return(plots.of.stats)
}

