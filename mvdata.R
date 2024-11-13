mvdata <- function(data,studlab,treat,n,event,reference.group,sm,...){
  
  # Read the data from a data frame
  
nulldata <- is.null(data)

sfsp <- sys.frame(sys.parent())

mc <- match.call()

if (nulldata) {
  
  data <- sfsp
  
}

studlab <- catch("studlab", mc, data, sfsp)

treat <- catch("treat", mc, data, sfsp)

event <- catch("event", mc, data, sfsp)

#mean <- catch("mean", mc, data, sfsp)

#sd <- catch("sd", mc, data, sfsp)

n <- catch("n", mc, data, sfsp)

args <- list(...) 
# 
nam.args <- names(args)

dat1 <- cbind.data.frame(studlab,treat,n,event[[1]])

dat2 <- cbind.data.frame(studlab,treat,n,event[[2]])

names(dat1)[ncol(dat1)] <- "event1"

names(dat2)[ncol(dat2)] <- "event2"

dat1$event2 <- dat2$event2

data_new <- dat1

dat <- structure(data = data_new, sm = sm)

jags_data <- jags_data(dat)

return(jags_data)
# 
}
