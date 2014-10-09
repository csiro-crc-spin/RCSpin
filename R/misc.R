organize_results<-function(object){
    if(class(object)=="CrcSpinModel"){
        number.M<-sum(unlist (lapply(object$study_group,f<-function(x){x$sex}))=="M")
        number.F<-sum(unlist(lapply(object$study_group,f<-function(x){x$sex}))=="F")
        results.F<-data.frame(object$study_results[(1:99)*2-1,] )
        results.M<-data.frame(object$study_results[(1:99)*2,]   )
        names(results.F)<-c(
            "adenoma.state.adenoma",
            "adenoma.state.large.adenoma",
            "adenoma.state.pre.symptomatic.CRC",
            "person.colon_clinical.CRC",
            "colon.state.clear",
            "colon.state.adenoma",
            "colon.state.large.adenoma",
            "colon.state.pre.symptomatic.CRC",
            "colon.state.CRC",
            "colon.cancer_site.cecum",
            "colon.cancer_site.ascending",
            "colon.cancer_site.transverse",
            "colon.cancer_site.descending",
            "colon.cancer_site.sigmoid",
            "colon.cancer_site.rectum",
            "person.state.deceased",
            "person.colon_clinical.CRC",
            "person.initiateCRCTreatment.called"
            )
        names(results.M)<-names(results.F)
        descriptions<-c(
            "A count of adenomas in all patients of the study group that have state 'adenoma'",
            "A count of adenomas in all patients of the study group that have state 'large adenoma'",
            "A count of adenomas in all patients of the study group that have state 'pre symptomatic CRC'",
            "A count of patients that have colon_clinical state 'CRC'",
            "A count of how many peoples colons were in state 'clear'",
            "A count of how many peoples colons were in state 'adenoma'",
            "A count of how many peoples colons were in state 'large adenoma'",
            "A count of how many peoples colons were in state 'pre symptomatic CRC'",
            "A count of how many peoples colons were in state 'CRC'",
            "A count of how many people's colons overall that have the majority of their adenomas in the 'cecum' location",
            "A count of how many people's colons overall that have the majority of their adenomas in the 'ascending' location",
            "A count of how many people's colons overall that have the majority of their adenomas in the 'transverse' location",
            "A count of how many people's colons overall that have the majority of their adenomas in the 'descending' location",
            "A count of how many people's colons overall that have the majority of their adenomas in the 'sigmoid' location",
            "A count of how many people's colons overall that have the majority of their adenomas in the 'rectum' location",
            "A count of people in the study that have died (i.e. their state is 'deceased')",
            "A count of patients that have colon_clinical state 'CRC'",
            "A count of people put into a treatment program THIS ITERATIN"
            )
        descriptions<-data.frame(names(results.F),descriptions)
        res<-list(number.M,number.F,results.M,results.F,descriptions)
        names(res)<-c("number.M","number.F","results.M","results.F","descriptions")
    }
    if(class(object)=="DukesCrcSpinModel"){
        number.M<-sum(unlist (lapply(object$study_group,f<-function(x){x$sex}))=="M")
        number.F<-sum(unlist(lapply(object$study_group,f<-function(x){x$sex}))=="F")
        results.F<-data.frame(object$study_results[(1:99)*2-1,] )
        results.M<-data.frame(object$study_results[(1:99)*2,]   )
        names(results.F)<-c("adenoma.state.adenoma",
                            "adenoma.state.large.adenoma",
                            "adenoma.state.CRC",
                            "adenoma.state.deceased",
                            "person.colon_clinical.symptomatic.CRC",
                            "colon.state.clear",
                            "colon.state.adenoma",
                            "colon.state.large.adenoma",
                            "colon.state.CRC",
                            "colon.state.symptomatic.CRC",
                            "colon.state.deceased",
                            "colon.state.CRC.adenoma.stage.A",
                            "colon.state.CRC.adenoma.stage.B",
                            "colon.state.CRC.adenoma.stage.C",
                            "colon.state.CRC.adenoma.stage.D",
                            "colon.state.CRC.adenoma.stage.deceased",
                            "colon.state.CRC.colon.cancer_site.cecum",
                            "colon.state.CRC.colon.cancer_site.ascending",
                            "colon.state.CRC.colon.cancer_site.transverse",
                            "colon.state.CRC.colon.cancer_site.descending",
                            "colon.state.CRC.colon.cancer_site.sigmoid",
                            "colon.state.CRC.colon.cancer_site.rectum",
                            "colon.state.symptomatic.CRC.adenoma.stage.A",
                            "colon.state.symptomatic.CRC.adenoma.stage.B",
                            "colon.state.symptomatic.CRC.adenoma.stage.C",
                            "colon.state.symptomatic.CRC.adenoma.stage.D",
                            "colon.state.symptomatic.CRC.adenoma.stage.deceased",
                            "colon.state.symptomatic.CRC.colon.cancer_site.cecum",
                            "colon.state.symptomatic.CRC.colon.cancer_site.ascending",
                            "colon.state.symptomatic.CRC.colon.cancer_site.transverse",
                            "colon.state.symptomatic.CRC.colon.cancer_site.descending",
                            "colon.state.symptomatic.CRC.colon.cancer_site.sigmoid",
                            "colon.state.symptomatic.CRC.colon.cancer_site.rectum",
                            "colon.cancer_site.cecum",
                            "colon.cancer_site.ascending",
                            "colon.cancer_site.transverse",
                            "colon.cancer_site.descending",
                            "colon.cancer_site.sigmoid",
                            "colon.cancer_site.rectum",
                            "colon.stage.A",
                            "colon.stage.B",
                            "colon.stage.C",
                            "colon.stage.D",
                            "colon.stage.deceased",
                            "person.state.deceased.person.state.deceased.CRC",
                            "person.colon_clinical.symptomatic.CRC",
                            "not.implemented",
                            "not.implemented",
                            "not.implemented",
                            "not.implemented",
                            "not.implemented",
                            "person.initiateCRCTreatment.colonoscopy_performed",
                            "person.initiateCRCTreatment.colon.state.large.adenoma",
                            "person.initiateCRCTreatment.colon.stage.A",
                            "person.initiateCRCTreatment.colon.stage.B",
                            "person.initiateCRCTreatment.colon.stage.C",
                            "person.initiateCRCTreatment.colon.stage.D",
                            "person.initiateCRCTreatment.person.in_treatment_program",
                            "person.initiateCRCTreatment.colonoscopy_caused_bleeding",
                            "person.initiateCRCTreatment.colonoscopy_caued_perforation"
                            )
        names(results.M)<-names(results.F)
        descriptions<-c("A count of adenomas in all patients of the study group that have state \"adenoma\"",
                        "A count of adenomas in all patients of the study group that have state \"large adenoma\"",
                        "A count of adenomas in all patients of the study group that have state \"CRC\"",
                        "A count of adenomas in all patients of the study group that have state \"deceased\"",
                        "A count of patients that have colon_clinical state \"symptomatic CRC\"",
                        "A count of how many peoples colons were in state \"clear\"",
                        "A count of how many peoples colons were in state \"adenoma\"",
                        "A count of how many peoples colons were in state \"large adenoma\"",
                        "A count of how many peoples colons were in state \"CRC\"",
                        "A count of how many peoples colons were in state \"symptomatic CRC\"",
                        "A count of how many peoples colons were in state \"deceased\"",
                        "A count of the adenomas that are at stage \"A\" in colons with state \"CRC\" across all people in the study group",
                        "A count of the adenomas that are at stage \"B\" in colons with state \"CRC\" across all people in the study group",
                        "A count of the adenomas that are at stage \"C\" in colons with state \"CRC\" across all people in the study group",
                        "A count of the adenomas that are at stage \"D\" in colons with state \"CRC\" across all people in the study group",
                        "A count of the adenomas that are at stage \"deceased\" in colons with state \"CRC\" in all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"cecum\" location across all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"ascending\" location across all people in the study grou",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"transverse\" location across all people in the study gro",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"descending\" location across all people in the study gro",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"sigmoid\" location across all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"rectum\" location across all people in the study group",
                        "A count of the adenomas that are at stage \"A\" in people of the study group with colon_clinical set to \"symptomatic CRC\"",
                        "A count of the adenomas that are at stage \"B\" in people of the study group with colon_clinical set to \"symptomatic CRC\"",
                        "A count of the adenomas that are at stage \"C\" in people of the study group with colon_clinical set to \"symptomatic CRC\"",
                        "A count of the adenomas that are at stage \"D\" in people of the study group with colon_clinical set to \"symptomatic CRC\"",
                        "A count of the adenomas that are at stage \"deceased\" in people of the study group with colon_clinical set to \"symptomatic CRC\"",
#                        "A count of the adenomas that are at stage \"deceased\" in people of the study group with colon_clinical set to \"symptomatic CRC\"",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"cecum\" location across all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"ascending\" location across all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"transverse\" location across all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"descending\" location across all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"sigmoid\" location across all people in the study group",
                        "A count of colons in \"state\" that have the majority of their adenomas in the \"rectum\" location across all people in the study group",
                        "A count of how many people's colons overall that have the majority of their adenomas in the \"cecum\" location",
                        "A count of how many people's colons overall that have the majority of their adenomas in the \"ascending\" location",
                        "A count of how many people's colons overall that have the majority of their adenomas in the \"transverse\" location",
                        "A count of how many people's colons overall that have the majority of their adenomas in the \"descending\" location",
                        "A count of how many people's colons overall that have the majority of their adenomas in the \"sigmoid\" location",
                        "A count of how many people's colons overall that have the majority of their adenomas in the \"rectum\" location",
                        "A count of how many people's colons overall were in stage \"A\"",
                        "A count of how many people's colons overall were in stage \"B\"",
                        "A count of how many people's colons overall were in stage \"C\"",
                        "A count of how many people's colons overall were in stage \"D\"",
                        "A count of how many people's colons overall were in stage \"deceased\"",
 #                       "A count of people in the study that have died (i.e. their state is \"deceased\")",
                        "A count of people in the study that have died (i.e. their state is \"deceased\" or \"deceased CRC\")",
                        "A count of patients that have colon_clinical state \"symptomatic CRC\"",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "A count of people that where treated for CRC IN THIS ITERATION, in which a colonoscopy was performed",
                        "A count of people that where treated for CRC IN THIS ITERATION, who had colons in state \"adenoma\" or \"large adenoma\"",
                        "A count of people that where treated for CRC IN THIS ITERATION, who had colons in state \"symptomatic CRC\" and in stage \"A\"",
                        "A count of people that where treated for CRC IN THIS ITERATION, who had colons in state \"symptomatic CRC\" and in stage \"B\"",
                        "A count of people that where treated for CRC IN THIS ITERATION, who had colons in state \"symptomatic CRC\" and in stage \"C\"",
                        "A count of people that where treated for CRC IN THIS ITERATION, who had colons in state \"symptomatic CRC\" and in stage \"D\"",
                        "A count of people put into a treatment program THIS ITERATION!!!",
                        "A count of people that where treated for CRC IN THIS ITERATION, whose colonoscopy caused bleeding",
                        "A count of people that where treated for CRC IN THIS ITERATION, whose colonoscopy caused perforation")
        
        descriptions<-data.frame(names(results.F),descriptions)
        names(descriptions)<-c("variable.name", "description")
        res<-list(number.M,number.F,results.M,results.F,descriptions)
        names(res)<-c("number.M","number.F","results.M","results.F","descriptions")
    }
    res
}

###############################################################################
## http://www.r-bloggers.com/data-mountains-and-streams-stacked-area-plots-in-r/
##http://menugget.blogspot.com.au/2013/12/data-mountains-and-streams-stacked-area.html
##plot.stacked makes a stacked plot where each y series is plotted on top
##of the each other using filled polygons
##
##Arguments include:
## 'x' - a vector of values
## 'y' - a matrix of data series (columns) corresponding to x
## 'order.method' = c("as.is", "max", "first")
## "as.is" - plot in order of y column
## "max" - plot in order of when each y series reaches maximum value
## "first" - plot in order of when each y series first value > 0
## 'col' - fill colors for polygons corresponding to y columns (will recycle)
## 'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
## 'lwd' - border line width for polygons corresponding to y columns (will recycle)
## '...' - other plot arguments

plot.stacked <- function(x, y,order.method = "as.is",ylab="", xlab="",
                         border = NULL, lwd=1,
                         col=rainbow(length(y[1,])),
                         ylim=NULL,
                         ...
                         ){
    
    if(sum(y < 0) > 0) error("y cannot contain negative numbers")
    
    if(is.null(border)) border <- par("fg")
    border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
    col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
    lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))
    
    if(order.method == "max") {
        ord <- order(apply(y, 2, which.max))
        y <- y[, ord]
        col <- col[ord]
        border <- border[ord]
    }
    
    if(order.method == "first") {
        ord <- order(apply(y, 2, function(x) min(which(r>0))))
        
        
        y <- y[, ord]
        col <- col[ord]
        border <- border[ord]
    }
    
    top.old <- x*0
    polys <- vector(mode="list", ncol(y))
    for(i in seq(polys)){
        top.new <- top.old + y[,i]
        polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
        top.old <- top.new
    }
    
    if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
    plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
    for(i in seq(polys)){
        polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
    }
    
}
#plot.stream makes a "stream plot" where each y series is plotted
#as stacked filled polygons on alternating sides of a baseline.
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first")
# "as.is" - plot in order of y column
# "max" - plot in order of when each y series reaches maximum value
# "first" - plot in order of when each y series first value > 0
#'center' - if TRUE, the stacked polygons will be centered so that the middle,
#i.e. baseline ("g0"), of the stream is approximately equal to zero.
#Centering is done before the addition of random wiggle to the baseline.
#'frac.rand' - fraction of the overall data "stream" range used to define the range of
#random wiggle (uniform distrubution) to be added to the baseline 'g0'
#'spar' - setting for smooth.spline function to make a smoothed version of baseline "g0"
#'col' - fill colors for polygons corresponding to y columns (will recycle)
#'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#'lwd' - border line width for polygons corresponding to y columns (will recycle)
#'...' - other plot arguments
 
plot.stream <- function( x, y, order.method = "as.is", frac.rand=0.1,
                        spar=0.2, center=TRUE, ylab="", xlab="", border = NULL, lwd=1,
                        col=rainbow(length(y[1,])), ylim=NULL, ...  ){
    
    if(sum(y < 0) > 0) error("y cannot contain negative numbers")
    
    if(is.null(border)) border <- par("fg")
    border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
    col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
    lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))
    
    if(order.method == "max") {
        ord <- order(apply(y, 2, which.max))
        y <- y[, ord]
        col <- col[ord]
        border <- border[ord]
    }
    
    if(order.method == "first") {
	ord <- order(apply(y, 2, function(x) min(which(r>0))))
        y <- y[, ord]
        col <- col[ord]
        border <- border[ord]
    }
    
    bottom.old <- x*0
    top.old <- x*0
    polys <- vector(mode="list", ncol(y))
    for(i in seq(polys)){
        if(i %% 2 == 1){ #if odd
            top.new <- top.old + y[,i]
            polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
            top.old <- top.new
        }
        if(i %% 2 == 0){ #if even
            bottom.new <- bottom.old - y[,i]
            polys[[i]] <- list(x=c(x, rev(x)), y=c(bottom.old, rev(bottom.new)))
            bottom.old <- bottom.new
        }
}
 
    ylim.tmp <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
    outer.lims <- sapply(polys, function(r) rev(r$y[(length(r$y)/2+1):length(r$y)]))
    mid <- apply(outer.lims, 1, function(r) mean(c(max(r, na.rm=TRUE), min(r, na.rm=TRUE)), na.rm=TRUE))
                                        #center and wiggle
    if(center) {
        g0 <- -mid + runif(length(x), min=frac.rand*ylim.tmp[1], max=frac.rand*ylim.tmp[2])
    } else {
        g0 <- runif(length(x), min=frac.rand*ylim.tmp[1], max=frac.rand*ylim.tmp[2])
    }
    fit <- smooth.spline(g0 ~ x, spar=spar)
    
    for(i in seq(polys)){
        polys[[i]]$y <- polys[[i]]$y + c(fit$y, rev(fit$y))
    }
    
    if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
    plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
    for(i in seq(polys)){
        polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
    }
    
}


#this function converts a vector of values("z") to a vector of color
#levels. One must define the number of colors. The limits of the color
#scale("zlim") or the break points for the color changes("breaks") can 
#also be defined. when breaks and zlim are defined, breaks overrides zlim.
val2col<-function(z, zlim, col = heat.colors(12), breaks){
    if(!missing(breaks)){
        if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
    }
    if(missing(breaks) & !missing(zlim)){
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
    }
    if(missing(breaks) & missing(zlim)){
        zlim <- range(z, na.rm=TRUE)
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    CUT <- cut(z, breaks=breaks)
    colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
    return(colorlevels)
}


## set.seed(1)
## m <- 500
## n <- 30
## x <- seq(m)
## y <- matrix(0, nrow=m, ncol=n)
## colnames(y) <- seq(n)
## for(i in seq(ncol(y))){
## mu <- runif(1, min=0.25*m, max=0.75*m)
## SD <- runif(1, min=5, max=20)
## TMP <- rnorm(1000, mean=mu, sd=SD)
## HIST <- hist(TMP, breaks=c(0,x), plot=FALSE)
## fit <- smooth.spline(HIST$counts ~ HIST$mids)
## y[,i] <- fit$y
## }
## y <- replace(y, y<0.01, 0)
 
 
## #Plot Ex. 1 - Color by max value
## pal <- colorRampPalette(c(rgb(0.85,0.85,1), rgb(0.2,0.2,0.7)))
## BREAKS <- pretty(apply(y,2,max),8)
## LEVS <- levels(cut(1, breaks=BREAKS))
## COLS <- pal(length(BREAKS )-1)
## z <- val2col(apply(y,2,max), col=COLS)
 
## #plot.stacked(x,y, xlim=c(100, 400), ylim=c(0, 1.2*max(apply(y,1,sum), na.rm=TRUE)),
## #             yaxs="i", col=z, border="white", lwd=0.5)
