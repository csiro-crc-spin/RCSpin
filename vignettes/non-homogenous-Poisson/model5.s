#Tuesday, October 22, 2013  model4.s -- modified to set input parameters to the values in
#column 2 (posterior means) fo table 2 for Rutter et al 2009

#Tuesday, September 17, 2013
# just realized htat I have not used Carolyn's sojourn time function. Forgotten
# this entirely. Should document more.

# here we will try to do a strict CRC-Spin model

#sojourn time is set as 
#sojourn.time<- sojourn.time.A+sojourn.time.B+ sojourn.time.C+ sojourn.time.D
# This is done in transition.adenomas.dukes.stages

## the standard sequence is 
## initiate.adenomas
## update.adenoma.growth.curve
## update.transitions calls  transition.adenomas.dukes.stages on the list of sites
## now calls transition.adenomas on the list of sites
## get.patient.state

## update.transitions.pre.clinical.to.clinical
## my model
## update.transitions.pre.clinical.to.clinical moves to clinical symptome with a probability based
## on the Dukes stage

## CRC-spin does
## year of transition.to.stage.A + sojourn.time where sojourn.time has a distribution described in
## Rutter wt al 2012.



## death.rate.other.causes
## result.person.4

## sojourn.time.A etc is use move between stages in transition.adenomas.dukes.stages
## update.transitions.pre.clinical.to.clinical is just done on probability of transition given stage.

## Here is what we will change.
## rename transition.adenomas.dukes.stages to transition.adenomas
## We will set sojourn.time in transition.adenomas (will transition adenomas to "pre symptomatic CRC")
## update.transitions.pre.clinical.to.clinical will be modified to reflect the CRC-Spin model
## ie 


## added 
## mu.colon and tau.colon
## to class adenoma
## and initialized them in "initiate.adenomas"


setClass("clinical.history",representation(status="character",
                                           events="list"))

setClass("adenoma",representation(initiated.in.year="numeric",size="numeric",location="character",
                                  state="character",
                                  stage="character",
                                  beta1="numeric",beta2="numeric",
                                  gamma1="numeric",
                                  gamma2="numeric",
                                  gamma3="numeric",
                                  d10="numeric",
                                  lambda="numeric",
                                  sojourn.time="numeric",
                                  transition.to.preclinical.crc.year="numeric",
                                  transition.to.crc.year="numeric",
                                  crc.survival.time="numeric",
                                  sojourn.time.A="numeric",
                                  sojourn.time.B="numeric",
                                  sojourn.time.C="numeric",
                                  sojourn.time.D="numeric",
                                  transition.to.stage.A="numeric",
                                  transition.to.stage.B="numeric",
                                  transition.to.stage.C="numeric",
                                  transition.to.stage.D="numeric",
                                  nu="numeric",
                                  xi="numeric"
                                  ))

setClass("risk.parameters",representation(baseline.risk="numeric",
                                          sex.linked.risk="numeric",
                                          age.risks="numeric"))


setClass("colon",representation(state="character",stage="character",cancer.site="character",
                                risk.parameters="risk.parameters",sites="list"))
# state takes the values "clear", "adenoma","large  adenoma",
# "non symptomatic CRC","CRC","deceased"
# "deceased" means "deceased CRC"
#  stage takes the values A", "B", "B", "D", "deceased"
# cancer.site takes the values "cecum","ascending","transverse","descending","sigmoid","rectum"

setClass("person",representation(age="numeric",sex="character",
                                 state="character", #"deceased","living",  "deceased CRC"
                                 colon.clinical="character", #"clear", "CRC"
                                 in.treatment.program="character",
                                 clinical.history="clinical.history",colon="colon",
                                 system.data="list"
                                 ))
# if colon.clinical = "CRC" then we have a diagnosed CRC. This can be the result of a
# test or is can be the backgroud level of diagnosis.

setClass("test",representation(age="numeric",test.type="character",compliance="character",
                               test.result="character",test.state="character"))

setClass("symptomatic_presentation",representation(age="numeric",cancer.stage="character"))



initialize.person<-function(sex=sample(c("F","M"),1,prob=c(0.5,0.5)),
                            risk="standard"){
  sex.linked.risk= -0.24 
  if(sex=="M"){
    sex.linked.risk<- -sex.linked.risk
  }
  if (risk=="high"){ 
    baseline.risk <- rnorm(1, mean = -4, sd = (0.27))
  }else {
    mu <-  -6.6 
    sigma <-  1.1
    baseline.risk <- rnorm(1, mean = mu, sd = sigma)
  }
  age.risks<-rep(0,4)
  age.risks[1] <- 0.037
  age.risks[2] <- 0.031
  age.risks[3] <- 0.029
  age.risks[4] <- 0.030
  
  person1<-new("person",age=1 ,sex=sex,state="living",colon.clinical="clear",
               in.treatment.program="no",
               clinical.history=new("clinical.history",events=vector("list", length=0)),
               colon=new("colon",
                 state="clear",
                 risk.parameters=new("risk.parameters",
                   baseline.risk= baseline.risk,
                   sex.linked.risk=sex.linked.risk,
                   age.risks=age.risks
                   )
                 )
               )
  person1
}

risk.of.an.adenoma<-function(person){
  age<-person@age
  sex<-person@sex
  #  if (sex=="M"){
  #    sex.factor<- -1
  #  }else{
  #    sex.factor<- 1
  #  }
  r1<-person@colon@risk.parameters@baseline.risk
  aa<-  person@colon@risk.parameters@age.risks
#  aa<-ifelse(aa>0,aa,0)
  #  r1<-r1+ sex.factor*person@colon@risk.parameters@sex.linked.risk
  r1<-r1+ person@colon@risk.parameters@sex.linked.risk
  r2<-0
  r3<-0
  if (age >= 20){
    # k<-as.numeric(cut(age,breaks=c(20,50,60,70,200),right=FALSE))
    # Replace inefficient cut code above with following if statement
    if (age>=70){
      k<-4
    } else if (age>=60){
      k<-3
    } else if (age>=50){
      k<-2
    } else {
      k<-1
    }
    A<-c(20,50,60,70,120)
    r2<-age*aa[k]
    
    r3<-0
    if ( k >= 2 ){
      for (j in 2:k){
        r3<-r3+A[j]*(aa[j-1]-aa[j])
      }
    }
  }
  aa<-exp(r1+r2+r3)
  ifelse((aa>1),1,aa)
}

initiate.adenomas<-function(person,study.year){
  sites<-person@colon@sites
  age<-person@age
  sex<-person@sex
  #  temp<-rpois(1, risk.of.an.adenoma(person))
  temp<-risk.of.an.adenoma(person)
  temp<-sample(c(0,1),1,prob=c(1-temp,temp))
  
  if (temp>0){
    for ( i in 1:temp){
      a1<-sample(c("cecum","ascending","transverse","descending","sigmoid","rectum"),
                 1, prob =c(0.08,0.23,0.24,0.12,0.24,0.09))
      if (a1=="rectum"){
        beta1<-10.3
        beta2<-2.7
        mu<-2.7
        tau<-0.84
        if(sex=="M"){
          gamma1<-0.035
          gamma2<-0.010
        }else{ #F
          gamma1<-0.043
          gamma2<-0.015
        }
      } else{ #colon
        beta1<-28.6
        beta2<-2.7
        mu<-1.9
        tau<-0.8
        if(sex=="M"){
          gamma1<-0.045
          gamma2<-0.008
        }else{ #F
          gamma1<-0.048
          gamma2<-0.005
        }
      }
      gamma3<-0.5
      d0<-1
      d10<-beta1*((-log(runif(1,0,1)))^(-1/beta2))
      dinfinity<-50
      nu<-sqrt(log(tau^2+1))
      xi<-log(mu)-0.5*nu^2
      lambda<- (-1/d10)*(log( (dinfinity-10)/(dinfinity-d0)))
      sites<-lappend(sites,new("adenoma",initiated.in.year=study.year,
                               size=1,state="adenoma",
                               location=a1,
                               beta1=beta1,
                               beta2=beta2,
                               d10=d10,
                               lambda=lambda,
                               gamma1=gamma1,
                               gamma2=gamma2,
                               gamma3=gamma3,
                               nu=nu,
                               xi=xi
                               )
                     )
    }
  }
  person@colon@sites<-sites
  person
}


update.adenoma.growth.curve<-function(object){
  object@colon@sites<-lapply(object@colon@sites,adenoma.growth.curve,object@age)
  object
}


adenoma.growth.curve<-function(object,age){
  #time is since the initiation of the adenoma i.e
  # person@age-object@initiated.in.year or
  if( class(object)!="adenoma"){cat("not an adenoma","\n")}
  age.of.adenoma <- age-object@initiated.in.year
  dinfinity<-50
  d0<-1
  size<-  dinfinity-(dinfinity-d0)*exp(-object@lambda*age.of.adenoma)
  object@size<-size
  object
}



update.transitions<-function(object){
  object@colon@sites<-lapply(object@colon@sites,transition.adenomas,object@age)
  object
}


get.patient.state<-function(object){
  #gets the list of colon sites
  #object@colon@state  c("adenoma","large adenoma","pre symptomatic CRC","CRC","deceased")
  # object@colon@stage<-stage
  
  ss<-object@colon@sites
  if (length(ss)>0){
    tt<-lapply(ss,f<-function(x){x@state})
    tt<-unlist(tt) # c("adenoma","large adenoma","pre clinical CRC","CRC","deceased")
    if (  sum(tt=="deceased")>0 ) {
      state<-"deceased"
    } else if ( sum(tt=="CRC")>0 ) {  
      state<-"CRC"
    } else if ( sum(tt=="pre symptomatic CRC")>0 ) {  
      state<-"pre symptomatic CRC"
    } else if ( sum(tt=="large adenoma")>0 ) {
      state<-"large adenoma"
    } else
    state<-"adenoma"
  }
  else{
    state<-"clear"
  }
  
  if (length(ss)>0){
    tt<-lapply(ss,f<-function(x){x@stage})
    tt<-unlist(tt) # c("A","B","C","D","deceased")
    if (  sum(tt=="deceased")>0 ) {
      stage<-"deceased"
    }else if (  sum(tt=="D")>0 ) {
      stage<-"D"
    } else if ( sum(tt=="C")>0 ) {
      stage<-"C"
    } else if ( sum(tt=="B")>0 ) {
      stage<-"B"
    } else if ( sum(tt=="A")>0 ) {
      stage<-"A"
    } else  stage<-"clear"
  }   else{
    stage<-"clear"
  }
  
  
  
  object@colon@state<-state
  object@colon@stage<-stage
  
  
  if (state=="CRC" | state=="pre symptomatic CRC"){
    tt<-lapply(ss,f<-function(x){x@location})
    tt<-unlist(tt) 
    object@colon@cancer.site<-tt[1]
  }
  #what if they have cancer in more than one site?
  
  if (object@colon@state=="deceased"){
    object@state<-"deceased CRC"
  }
  
  object
}



#http://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/3302.0.55.0012008-2010?OpenDocument
#
#lx the number of persons surviving to exact age x;
#qx the proportion of persons dying between exact age x and exact age x+1. 
#It is the mortality rate, from which other functions of the life table are derived;
#Lx the number of person years lived within the age interval x to x+1; and
#ex life expectancy at exact age x.

#aa<-read.table(file="./script_CR_figs/3302055001DO001_20082010.csv",sep=",",skip=5,header=TRUE,as.is=TRUE)
#aa<-aa[-c(1,103,104,105),]
#death.rate.male<-as.numeric(aa$qx)
#death.rate.female<-as.numeric(aa$qx.1)

load("../data/death_rates.RData")
death.rate.male<-death_rate_male
death.rate.female<-death_rate_female

death.rate.other.causes <- function(object){
  #check that the person is not dead already
  if(!is.element(object@state,c("deceased","deceased CRC"))){
    
    age<-object@age
    sex<-object@sex
    if(sex=="F"){
      state<-sample(c("deceased","living"),1, prob =c( death.rate.female[age], 1-death.rate.female[age]))
    }
    if(sex=="M"){
      state<-sample(c("deceased","living"),1, prob =c( death.rate.male[age], 1-death.rate.male[age]))
    }
    object@state<-state
  }
  object
}


lappend <- function(my.list, obj) {
  my.list[[length(my.list)+1]] <- obj
  return(my.list)
}


summary.person<-function(object,all.adenomas=FALSE){
  cat("age is",object@age,"\n")
  cat("number of adenomas", length(object@colon@sites),"\n")
  if (all.adenomas==TRUE){
    ss<-object@colon@sites
    if (length(ss)>0){
      #    for (i in 1:length(ss)){
      #      temp<-ss[[i]]
      #      cat(temp@location, " ", temp@stage," ", "years ",study.year-temp@initiated.in.year," ","\n")
      #    }
      #  }
      temp<-lapply(ss,f<-function(x){x@size})
      cat("size of largest adenoma",  max(unlist(temp)),"\n")
      aa<-lapply(ss,f<-function(x){x@location})
      cat("locations of adenomas:",unlist(aa),"\n")
      aa<-lapply(ss,f<-function(x){x@state})
      cat("state of adenomas:",unlist(aa),"\n")
    }
  }
  cat("colon state:",object@colon@state,"\n")
  cat("colon stage:",object@colon@stage,"\n")
  cat("clinical:",object@colon.clinical,"\n")
  cat("person state:",object@state,"\n")
  
}

setMethod("summary","person",summary.person)



transition.adenomas<-function(object,age){
  ## object is of class colon@site.
  #The function looks after 
  # adenoma to large adenoma
  # large adenoma to  "pre symptomatic CRC"
  # the transtition to "CRC" is done in the function
  #  update.transitions.pre.clinical.to.clinical
  
  size<- object@size
  initiated.at.age <- object@initiated.in.year
  
  if (size >=10 & object@state=="adenoma"){
    object@state<-"large adenoma"
  }
  
  #if it is already a "pre symptomatic CRC" then there is nothing to do
  #if it is an adenoma (or large adenoma)  we check to see if it transitions to a "pre clinical CRC"
  if ( is.element(object@state , c("adenoma", "large adenoma"))){
    p1<-pnorm( (log(object@gamma1*size)+object@gamma2*(initiated.at.age-50))/object@gamma3)
    dice<-sample(c("transition","no transition"),1,prob=c(p1,1-p1))
    if (dice =="transition"){
      object@state<-"pre symptomatic CRC"
      object@transition.to.preclinical.crc.year<-study.year
      object@transition.to.stage.A<-study.year
      object@location
      object@sojourn.time <-exp(rnorm(1, mean = object@xi, sd =  object@nu))
    }
  } #end of transition to 'pre symptomatic CRC"
  object
}





## update.transitions.pre.clinical.to.clinical<-function(object){
  
##   # if the patient is showing symptoms then 
##   # changes object@colon.clinical to "CRC"
##   # object@colon@state to "CRC"
##   # add a new event to object@clinical.history@events  
  
##   #check to  see that the colon.clinical is not already "CRC". If it is there is nothing to do
##   if (object@colon.clinical=="CRC"){
##     return
##   }
  
  
##   if (object@colon@state =="pre symptomatic CRC"){
    
##     ss<-object@colon@sites
##     if (length(ss)>0){
##       for ( i in 1:length(ss)){
##         if (object@transition.to.stage.A+object@sojourn.time  <= study.year){
##           object@state<-"CRC"
##         }
##       }
##     }
    
    
##   }
##   object
## }


update.transitions.pre.clinical.to.clinical<-function(object){
  
  # if the patient is showing symptoms then 
  # changes object@colon.clinical to "CRC"
  # object@colon@state to "CRC"
  # add a new event to object@clinical.history@events  
  
  #check to  see that the colon.clinical is not already "CRC". If it is there is nothing to do
  if (object@colon.clinical=="CRC"){
    return
  }
  
  
  if (object@colon@state =="pre symptomatic CRC"){
    
    ss<-object@colon@sites
    if (length(ss)>0){
      for ( i in 1:length(ss)){
        temp<-ss[[i]]
        if(temp@state=="pre symptomatic CRC"){
          if (temp@transition.to.stage.A+temp@sojourn.time  <= study.year){
            object@state<-"CRC"
            object@colon.clinical <- "CRC"
            object@colon@state <-   "CRC"
          }
        }
      }
    }
    
    
  }
  object
}




## ################################################################################






## mean.rectum<-2.5
## tau.rectum<-1.25
## aa<-rep(0,10000)
## for ( i in 1:10000){

## aa[i]<-rlnorm(1, meanlog = mean.rectum, sdlog =  tau.rectum*mean.rectum)
## }
## quantile(aa,probs = c(0.05,0.5,0.95)) #7.745306e-04 1.530951e+02 2.905131e+07 
## #median should be in (0.3,5) for 
## #  mean.rectum<-runif(1,0.5,5)
## #  tau.rectum<-runif(1,0.1,1.5)
## #something is wrong here
## # so we make
## mean.rectum<-runif(1,1,6)
## tau.rectum<-runif(1,1,1.5)
## aa<-rep(0,10000)
## for ( i in 1:10000){
## mean.rectum<-5
## tau.rectum<-1.5
##   aa[i]<-rlnorm(1, meanlog = log(mean.rectum), sdlog =  log(tau.rectum*mean.rectum))
## }
## quantile(aa,probs = c(0.05,0.5,0.95))  #0.1723114   5.0042261 134.3565766 
## better

## for ( i in 1:10000){
## #  mean.rectum<-runif(1,0.5,5)
## #  tau.rectum<-runif(1,0.1,1.5)
##   mean.rectum<-5
##   tau.rectum<-1.5
##   nu<-sqrt(log(tau.rectum^2+1))
##   xi<-log(mean.rectum)-0.5*nu^2
##   aa[i]<-sojourn.time<-exp(rnorm(1, mean = xi, sd =  nu))
## }


## for ( i in 1:50000){
##   mean.rectum<-5
##   tau.rectum<-1.5
##   nu<-sqrt(log(tau.rectum^2+1))
##   xi<-log(mean.rectum)-0.5*nu^2
##   aa[i]<-sojourn.time<-exp(rnorm(1, mean = xi, sd =  nu))
## }
## quantile(aa,probs = c(0.05,0.5,0.95)) 
## #0.4732065  2.9050271 17.1558212 

## for ( i in 1:50000){
##   mean.rectum<-5
##   tau.rectum<-0.1
##   nu<-sqrt(log(tau.rectum^2+1))
##   xi<-log(mean.rectum)-0.5*nu^2
##   aa[i]<-sojourn.time<-exp(rnorm(1, mean = xi, sd =  nu))
## }
## quantile(aa,probs = c(0.05,0.5,0.95)) 
## #4.227341 4.975202 5.865645 

## for ( i in 1:10000){
##   mean.rectum<-0.5
##   tau.rectum<-1.5
##   nu<-sqrt(log(tau.rectum^2+1))
##   xi<-log(mean.rectum)-0.5*nu^2
##   aa[i]<-sojourn.time<-exp(rnorm(1, mean = xi, sd =  nu))
## }
## quantile(aa,probs = c(0.05,0.5,0.95)) 
## #0.04787084 0.27418899 1.62336030 

## for ( i in 1:10000){
##   mean.rectum<-0.5
##   tau.rectum<-0.1
##   nu<-sqrt(log(tau.rectum^2+1))
##   xi<-log(mean.rectum)-0.5*nu^2
##   aa[i]<-sojourn.time<-exp(rnorm(1, mean = xi, sd =  nu))
## }
## quantile(aa,probs = c(0.05,0.5,0.95)) 
## #0.4218187 0.4977533 0.5859797
                                  
## sojourn.time<-exp(rnorm(1, mean = xi, sd =  nu))







