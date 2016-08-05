###############################
#
# Duke's CRC spin model classes
#
###############################

# forward definitions
DukesAdenomaParams <- setRefClass( "DukesAdenomaParams")
DukesPersonWithColon <- setRefClass( "DukesPersonWithColon")

# Duke's CRC model class
########################
#' Class orchestrating the Duke's CRC  Spin model
#'
#' This class acts as the entry point and driver for an implementation of the
#' Dukes CRC Spin model, described by X et al., which is implemented using a
#' collection of classes.
#'
#' It orchestrates the modeling on a population represented as a
#' \code{study_group} of \code{\link{PersonWithColon}} class objects. And models
#' the development of \code{\link{Adenoma}}s, and ultimately CRC, within the
#' subjects' \code{\link{Colon}}.
#'
#' Users of the model should only call this class's \code{new}
#' (\code{\link{initialize}}), \code{\link{set_adenoma_modeling_parameters}},
#' \code{\link{set_crcrisk_modeling_parameters}}, and \code{\link{run}}
#' methods. And directly access it's \code{\link{study_results field}}.
#'
#' @field study_results A two dimensional matrix containing an 18 column summary of model
#'        state at each year of the study. Inherited from \code{\link{GenericModel}}
#' 
#' @export DukesCrcSpinModel
#' @exportClass DukesCrcSpinModel
#' 
#' @family DukesCrcSpinModel_classes


DukesCrcSpinModel <- setRefClass( "DukesCrcSpinModel",
            contains="CrcSpinModel",

    methods = list(

        initialize = function (...) {

            callSuper(...)

        },

        # Overwrite the functions that tell the superclass model
        # which kind of AdenomaParams and PersonWithColon
        # derived/subclasses to use
        adenomaParamsType = function () {
            return(getRefClass("DukesAdenomaParams"))
        },

        personWithColonType = function () {
            return(getRefClass("DukesPersonWithColon"))
        },


        screening.colonoscopy = function (person) {
            temp1<-rep(FALSE,person$NBCSPRecordSize())
            up.to.date<-FALSE
            do.test <- "decline"
            test.result <- "negative"
            test.state <- "decline"
            
#            if (person$age==61){browser()}
            
            ##has the peson had a colonoscopy on the past 10 years
            if (length(person$clinical_history$events) >0) {
                cc<-rev(unlist(lapply(person$clinical_history$events,f<-function(x){x$compliance})))
                                        #                cc<-match("accept",rev(unlist(lapply(person$clinical_history$events,f<-function(x){x$compliance}))))
                                        #                if(!is.na(cc)){
                aa<-rev(unlist(lapply(person$clinical_history$events,f<-function(x){x$type})))
                which((cc=="accept")&(aa=="colonoscopy"))[1]
                bb<-unlist(rev(lapply(person$clinical_history$events,f<-function(x){x$age}))[which((cc=="accept")&(aa=="colonoscopy"))[1]])
                up.to.date <- (person$age - bb < 9)
            }

#            print(paste(up.to.date," ",person$age,sep=""))
            if (!up.to.date){
                uu<-person$BSA.propensity
                ww<-age.specific.compliance.rates.for.BSA(person)*400
                mm<-min(1,max(0,qlnorm(uu,mean=log(ww),sd=1.1)))
                aa1<-sample(c(1,0),1,prob=c(mm,1-mm )) 
                do.test<-sample(c("accept","decline"),1, prob =c(aa1,1-aa1))
                temp1[12]<-1  #test is offered
                
            }
            
 #           print(do.test)
            
            
            if ( (do.test=="accept") & ( person$colon_clinical=="clear") #
                &(person$in_treatment_program=="no")){
                person$updateState()  #object<-get.patient.state(object)
                state<-person$colon$state    #object@colon@state
                if( state=="symptomatic CRC" ){
                    sensitivity<-1
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                } else if ( state== "CRC" ){
                    sensitivity<-1
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                } else if ( state=="large adenoma" ){
                    sensitivity<-0.8
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                } else if ( state=="adenoma" ){
                    specificity<-0.6
                    test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
                    if(test.result=="positive"){
                        test.state<-"FP"
                    }
                    else{
                        test.state<-"TN"
                    }
                } else if ( state=="clear" ){
                    specificity<-0.99
                    test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
                    if(test.result=="positive"){
                        test.state<-"FP"
                    }
                    else{
                        test.state<-"TN"
                    }
                }#end state =clear
            


                temp1[6]<-1 #person has a colonoscopy woth probability 1
                temp1[13]<-sample(c(0,1),1,prob=c(0.9997,0.0003)) #probability of bleeding
                temp1[14]<-sample(c(0,1),1,prob=c(0.9999,0.0001)) #probability of perforation

            
                person$clinical_history$events<-lappend(person$clinical_history$events,
                                                        Test$new(
                                                                 age=person$age,
                                                                 type="colonoscopy",
                                                                 compliance=do.test,
                                                                 result=test.result,
                                                                 state= test.state)
                                                        )
                
                
              }               
            if(test.result=="positive"){  #if it is
              if ( (person$colon$state=="adenoma") | (person$colon$state=="large adenoma")){   #this may be wrong. 
                temp1[7]<-1
                temp1[12]<-1
                        person$in_treatment_program<-"yes"
                    } else if (person$colon$state=="CRC"){   #has this been changed to just "CRC" ??  yes **Changed
                        if (person$colon$stage=="A"){
                            temp1[8]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="B"){
                            temp1[9]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="C"){
                            temp1[10]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if(person$colon$stage=="D"){
###            person$colon.clinical<-"CRC"
                            temp1[11]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                                        #we do nothing.
                        }
                    }
            } #end test.result=="positive"
            temp1
        },
        

        BSA = function (person) {
            temp1<-rep(FALSE,person$NBCSPRecordSize())
            
            
            not.up.to.date<-TRUE
            if (length(person$clinical_history$events) >0){
            aa<-rev(lapply(person$clinical_history$events,f<-function(x){x$type}))
            bb<-rev(lapply(person$clinical_history$events,f<-function(x){x$age}))
            not.up.to.date <- (person$age - unlist(bb[match("iFOBT",aa)]) > 1)
            }
            
            if (  ( person$colon_clinical=="clear")            &(person$in_treatment_program=="no") & (not.up.to.date)) {
                uu<-person$BSA.propensity
                ww<-age.specific.compliance.rates.for.BSA(person)*10.0   
                mm<-min(1,max(0,qlnorm(uu,mean=log(ww),sd=1.1)))
                aa1<-sample(c(1,0),1,prob=c(mm,1-mm )) 
#              do.test<-sample(c("accept","decline"),1, prob =c(aa1,1-aa1))
#                temp1[1]<-1 # a test was offered
#            }
            
#                if (do.test=="accept"){
                    
                    iFOBTscreening(person,aa1) #same as .self$iFOBT.screening(person)
                    
                                        #offer iFOBT. Relevant parameters are the compliance rate and the
                                        #sensitivity and specificity, depending on the person1@colon@state and stage
                                        #The test results are retained in an object of class "test", appended to the list
                                        #person1@clinical.history@events
                    
                    test.outcome<-tail(person$clinical_history$events,1)[[1]] #returns a list -- the first item of which is a test
                    temp1[1]<-ifelse(is.element(test.outcome$type,c("iFOBT")),1,0)
                    temp1[2]<-ifelse(is.element(test.outcome$compliance,c("accept")),1,0)
                    temp1[3]<-ifelse(is.element(test.outcome$type,c("blood")),1,0)
                    temp1[4]<-ifelse(is.element(test.outcome$compliance,c("screen")),1,0)
###assumes that they only have one test. Needs to be changed
###we are also assuming that if the test is positive than the person has a colonoscopy. This
###is not the case --  0.938 go on to a colonoscopy (Cronin et al 2010)
                    if(test.outcome$result=="positive"){  #if it is a false positive, then we will skip over all the
                                        #conditions on person1@colon@state and move on
                        temp1[5]<-1 #test.result=="positive"
                        temp1[6]<-1 #person has a colonoscopy woth probability 1
                        temp1[13]<-sample(c(0,1),1,prob=c(0.9997,0.0003)) #probability of bleeding
                        temp1[14]<-sample(c(0,1),1,prob=c(0.9999,0.0001)) #probability of perforation
                        if ( (person$colon$state=="adenoma") | (person$colon$state=="large adenoma")){   #this may be wrong. 
                            temp1[7]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        } else if (person$colon$state=="CRC"){   #has this been changed to just "CRC" ??  yes **Changed
                            if (person$colon$stage=="A"){
                                temp1[8]<-1
                                temp1[12]<-1
                                person$in_treatment_program<-"yes"
                            }
                            if (person$colon$stage=="B"){
                                temp1[9]<-1
                                temp1[12]<-1
                                person$in_treatment_program<-"yes"
                            }
                            if (person$colon$stage=="C"){
                                temp1[10]<-1
                                temp1[12]<-1
                                person$in_treatment_program<-"yes"
                            }
                            if(person$colon$stage=="D"){
###            person$colon.clinical<-"CRC"
                                temp1[11]<-1
                                temp1[12]<-1
                                person$in_treatment_program<-"yes"
                                        #we do nothing.
                            }
                        }
                    } #end test.result=="positive"
            } #end person$colon==clear
            temp1
        },

        NBCSP = function (person) {
            temp1<-rep(FALSE,person$NBCSPRecordSize())
            if ( (person$age %in% c(55,60,65,70,72)) & ( person$colon_clinical=="clear") #
                &(person$in_treatment_program=="no")){
                                        #the current screening scheme offers iFOBT to people at the ages 50,55,60,65,70.
                                        #We do not offer it if the person already has a diagnosis of "CRC"
                uu<-person$NBCSP.propensity
                if (person$age ==55){
                    ww<- .359
                }else if (person$age == 60){
                    ww<- .427
                }  else {
                    ww<- 0.4
                }
                mm<-min(1,max(0,qlnorm(uu,mean=log(ww),sd=0.7)))
                aa1<-sample(c(1,0),1,prob=c(mm,1-mm )) 
                
                iFOBTscreening(person,aa1) #same as .self$iFOBT.screening(person)

                                        #offer iFOBT. Relevant parameters are the compliance rate and the
                                        #sensitivity and specificity, depending on the person1@colon@state and stage
                                        #The test results are retained in an object of class "test", appended to the list
                                        #person1@clinical.history@events

                test.outcome<-tail(person$clinical_history$events,1)[[1]] #returns a list -- the first item of which is a Test
                temp1[1]<-ifelse(is.element(test.outcome$type,c("iFOBT")),1,0)
                temp1[2]<-ifelse(is.element(test.outcome$compliance,c("accept")),1,0)
                temp1[3]<-ifelse(is.element(test.outcome$type,c("blood")),1,0)
                temp1[4]<-ifelse(is.element(test.outcome$compliance,c("screen")),1,0)
###assumes that they only have one test. Needs to be changed
###we are also assuming that if the test is positive than the person has a colonoscopy. This
###is not the case --  0.938 go on to a colonoscopy (Cronin et al 2010)
                if(test.outcome$result=="positive"){  #if it is a false positive, then we will skip over all the
                                        #conditions on person1@colon@state and move on
                    temp1[5]<-1 #test.result=="positive"
                    temp1[6]<-1 #person has a colonoscopy woth probability 1
                    temp1[13]<-sample(c(0,1),1,prob=c(0.9997,0.0003)) #probability of bleeding
                    temp1[14]<-sample(c(0,1),1,prob=c(0.9999,0.0001)) #probability of perforation
                    if ( (person$colon$state=="adenoma") | (person$colon$state=="large adenoma")){   #this may be wrong. 
                        temp1[7]<-1
                        temp1[12]<-1
                        person$in_treatment_program<-"yes"
                    } else if (person$colon$state=="CRC"){   #has this been changed to just "CRC" ??  yes **Changed
                        if (person$colon$stage=="A"){
                            temp1[8]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="B"){
                            temp1[9]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="C"){
                            temp1[10]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if(person$colon$stage=="D"){
###            person$colon.clinical<-"CRC"
                            temp1[11]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                                        #we do nothing.
                        }
                    }
                } #end test.result=="positive"
                
            }
            temp1
        },
        
        iFOBTscreening = function(person,aa1) {#person is offered iFOBT
            age<-person$age
            test.result<-"none"
            test.state<-"none"
#            uu<-person$NBCSP.propensity
#            if (age ==55){
#                ww<- .359
#            }else if (age == 60){
#                ww<- .427
#            }  else {
#                ww<- 0.4
#            }
#            mm<-min(1,max(0,qlnorm(uu,mean=log(ww),sd=0.7)))
#            aa1<-sample(c(1,0),1,prob=c(mm,1-mm )) 
            compliance<-sample(c("accept","decline"),1, prob =c(aa1,1-aa1))
            if (compliance=="accept"){
                person$updateState()  #object<-get.patient.state(object)
                state<-person$colon$state    #object@colon@state
                if( state=="symptomatic CRC" ){
                    sensitivity<-0.86
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                } else if ( state== "CRC" ){
                    sensitivity<-0.86
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                } else if ( state=="large adenoma" ){
                    sensitivity<-0.474
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                }else if ( state=="adenoma" ){
                    specificity<-0.9585
                    test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
                    if(test.result=="positive"){
                        test.state<-"FP"
                    }
                    else{
                        test.state<-"TN"
                    }
                } else if ( state=="clear" ){
                    specificity<-0.968
                    test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
                    if(test.result=="positive"){
                        test.state<-"FP"
                    }
                    else{
                        test.state<-"TN"
                    }
                }#end state =clear
            }
            person$clinical_history$events<-lappend(person$clinical_history$events,
                                                    Test$new(
                                                        age=age,
                                                        type="iFOBT",
                                                        compliance=compliance,
                                                        result=test.result,
                                                        state= test.state)
                                                    )
### The figures for sensitivity and specificity are from the Goodall report.
### Table from the report is not clear
### specificity 0.9585 P(test.result="negative"| state="adenoma")
### specificity 0.9680 P(test.result="negative"| state="clear")
### sensitivity 0.86   P(test.result="positive"| state="CRC")
### sensitivity 0.47  P(test.result="positive"| state="large adenoma")
            
        },

        modelSubjectTesting = function (person) {
            # test for and treat CRC
            #   Call function containing what was originally:
            #       if ((person1@colon.clinical=="CRC")
            #           & (person1@in.treatment.program=="no")){
            #               ...some code...
            #       }
            treatment_record.1 <- testForAndTreatCRC(person)


            if (screening_flag=="none"){
                treatment_record.2<-rep(0,14)
            } else if (screening_flag=="NBCSP"){
                       treatment_record.2<-NBCSP(person)
            } else if (screening_flag=="screening.colonoscopy"){
                treatment_record.2<-screening.colonoscopy(person)
            } else if (screening_flag=="gemini.screening"){
                treatment_record.2<-gemini.screening(person)
            } else if (screening_flag=="BSA"){
                treatment_record.2<-BSA(person)
            } else if (screening_flag=="GP.screening"){
               treatment_record.2<-GP.screening(person)
            } 
            
            return(c(treatment_record.1,treatment_record.2))
        },

        gemini.screening = function(person){
            temp1<-rep(0,14)
            do.test<-sample(c("accept","decline"),1, prob =c(0.03,1-0.03))
            if ( (person$age %in% c(50:80)) & ( person$colon_clinical=="clear") &(person$in_treatment_program=="no")   &(do.test=="accept")){
                
                                        #the current screening scheme offers iFOBT to people at the ages 50,55,60,65,70.
                                        #We do not offer it if the person already has a diagnosis of "CRC"
                gemini.blood.screening(person) #  blood.test.screening(person) #
                
                                        #offer Gemini test. Relevant parameters are the compliance rate and the
                                        #sensitivity and specificity, depending on the person1@colon@state and stage
                                        #The test results are retained in an object of class "test", appended to the list
                                        #person1@clinical.history@events
                
                test.outcome<-tail(person$clinical_history$events,1)[[1]] #returns a list -- the first item of which is a Test
                temp1[1]<-ifelse(is.element(test.outcome$type,c("iFOBT")),1,0)
                temp1[2]<-ifelse(is.element(test.outcome$compliance,c("accept")),1,0)
                temp1[3]<-ifelse(is.element(test.outcome$type,c("blood")),1,0)
                temp1[4]<-ifelse(is.element(test.outcome$compliance,c("screen")),1,0)
###assumes that they only have one test. Needs to be changed
###we are also assuming that if the test is positive than the person has a colonoscopy. This
###is not the case --  0.938 go on to a colonoscopy (Cronin et al 2010)
                if(test.outcome$result=="positive"){  #if it is a false positive, then we will skip over all the
                                        #conditions on person1@colon@state and move on
                    temp1[5]<-1 #test.result=="positive"
                    temp1[6]<-1 #person has a colonoscopy woth probability 1
                    temp1[13]<-sample(c(0,1),1,prob=c(0.9997,0.0003)) #probability of bleeding
                    temp1[14]<-sample(c(0,1),1,prob=c(0.9999,0.0001)) #probability of perforation
                    #so these are the results of colonoscopy. 
                    if ( (person$colon$state=="adenoma") | (person$colon$state=="large adenoma")){   #this may be wrong. 
                        temp1[7]<-1
                        temp1[12]<-1
                        person$in_treatment_program<-"yes"
                    } else if (person$colon$state=="CRC"){   #has this been changed to just "CRC" ??  yes **Changed
                        if (person$colon$stage=="A"){
                            temp1[8]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="B"){
                            temp1[9]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="C"){
                            temp1[10]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if(person$colon$stage=="D"){
###            person$colon.clinical<-"CRC"
                            temp1[11]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                                        #we do nothing.
                        }
                    }
                } #end test.result=="positive"
                
            }
            temp1   
        }, #end gemini.testing
        
        gemini.blood.screening = function(person){
            #person is offered blood screening. Not everyone is offered the blood test every year so
            # this probability includes visiting the doctor, being offered the test and accepting the test
            #
            #where od I want to do -- offer and compliance? 
            test.result<-"none"
            test.state<-"none"
            person$updateState() 
            state<-person$colon$state
            stage<-person$colon$stage
            if (state=="CRC" | state=="pre symptomatic CRC"){
                if (stage =="A"){
                    sensitivity<-0.50
                    specificity<-0.0
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                }#end stage A
                if (stage =="B"){
                    sensitivity<-0.68
                    specificity<-0.0
                    sensitivity<-0.68
                    specificity<-0.0
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                }#end stage B
                if (stage =="C"){
                    sensitivity<-0.80
                    specificity<-0.0
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                }#end stage C
                if (stage =="D"){
                    sensitivity<-1.0
                    specificity<-0.0
                    test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                    if(test.result=="positive"){
                        test.state<-"TP"
                    }
                    else{
                        test.state<-"FN"
                    }
                }#end stage D
            } else if ( state=="large adenoma" ){
                sensitivity<-0.33
                specificity<-0.0
                test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
                if(test.result=="positive"){
                    test.state<-"TP"
                }
                else{
                    test.state<-"FN"
                }
            } else if ( state=="adenoma" ){
                sensitivity<-0.0
                specificity<-0.93
                test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
                if(test.result=="positive"){
                    test.state<-"FP"
                }
                else{
                    test.state<-"TN"
                }
            }  else if ( state=="clear" ){
                sensitivity<-0.0
                specificity<-0.93
                test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
                if(test.result=="positive"){
                    test.state<-"FP"
                }
                else{
                    test.state<-"TN"
                }
            }#end state =clear
            person$clinical_history$events<-lappend(person$clinical_history$events,
                                                    Test$new(
                                                        age=person$age,
                                                        type="blood",
                                                        compliance="accept",
                                                        result=test.result,
                                                        state= test.state)
                                                    )
        }, # blood.test.screening 


        age.specific.compliance.rates.for.BSA = function(person){
            compliance<- NA
            age<-person$age
            compliance.rates<-structure(c(25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 
                                          38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 
                                          54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 
                                          70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 
                                          86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 0.000112114391175735, 
                                          0.000144059006569091, 0.000140745953553835, 0.000178600284059499, 
                                          0.000210958002480866, 0.000206700493815974, 0.000256992454023768, 
                                          0.000276503597457331, 0.000341535126145256, 0.000449870277132416, 
                                          0.000416345996139337, 0.000560955558779438, 0.000581870844068312, 
                                          0.000723750319064036, 0.0010183234577173, 0.00139337414745452, 
                                          0.0013756079152236, 0.001501001641173, 0.00149055311606172, 0.0018664612442872, 
                                          0.002099961903346, 0.00238890000839343, 0.00284066511200809, 
                                          0.00352980948971881, 0.00435570972446687, 0.00331614443793592, 
                                          0.0051611514117174, 0.00652146344556464, 0.00769247906608148, 
                                          0.0100330702644304, 0.00575132800343403, 0.00795829796884718, 
                                          0.00923958919527293, 0.0120125261258177, 0.0145806052937455, 
                                          0.0165975418926515, 0.0195853065994219, 0.0199481276967024, 0.0206706123967076, 
                                          0.0181319021570854, 0.0113700676161581, 0.0134075236037529, 0.0136003269309358, 
                                          0.0164902793122832, 0.016571732361206, 0.0172908205144605, 0.0193048819532556, 
                                          0.0184263552717827, 0.0175446654112892, 0.0164079488006146, 0.015579392613731, 
                                          0.015072929791066, 0.0140010770059235, 0.0138120591896003, 0.0129765106350968, 
                                          0.0127438624792055, 0.0120982342541907, 0.0104820138172, 0.00967267930902822, 
                                          0.00856947975444906, 0.00762470568909818, 0.00664491784055204, 
                                          0.00605023217977833, 0.00478679610163326, 0.00413796188601593, 
                                          0.00321031545992412, 0.00246481157164169, 0.00166608547405944, 
                                          0.00113300492610837, 0.000945562609752803, 0.000977623289159244, 
                                          0.000863309352517986, 0.00118764845605701, 0.00107411385606874, 
                                          0.00107411385606874), .Dim = c(75L, 2L), .Dimnames = list(NULL, 
                                                                                                    c("age", "compliance")))
            
            temp<-match(age,compliance.rates[,1])
            if (is.na(temp)){
                compliance<-0
            }else {
                compliance<- compliance.rates[temp,2]
            }
            compliance
        },

       GP.screening = function(person){
             if (person$age < 45){
               return(NULL)
                 }
            temp1<-rep(FALSE,person$NBCSPRecordSize())
            
            
            not.up.to.date<-TRUE
             if (length(person$clinical_history$events) >0){
                 aa<-rev(lapply(person$clinical_history$events,f<-function(x){x$type}))
                 bb<-rev(lapply(person$clinical_history$events,f<-function(x){x$age}))
                 not.up.to.date <- (person$age - unlist(bb[match("iFOBT",aa)]) > 1)
             }
            
            if (  ( person$colon_clinical=="clear")            &(person$in_treatment_program=="no") & (not.up.to.date)) {
                
                if (person$age >= 85){
                    aa1 <- 0.01674272
                } else if (person$age >= 75){
                    aa1 <- 0.03297858
                } else if (person$age >= 65){
                    aa1 <- 0.03094770
                } else if (person$age >= 55){
                    aa1 <- 0.02790423 
                } else if (person$age >= 45){
                    aa1 <-  0.01605888
                }
                
                
                iFOBTscreening(person,aa1) #same as .self$iFOBT.screening(person)
                
                                        #offer iFOBT. Relevant parameters are the compliance rate and the
                                        #sensitivity and specificity, depending on the person1@colon@state and stage
                                        #The test results are retained in an object of class "test", appended to the list
                                        #person1@clinical.history@events
                
                test.outcome<-tail(person$clinical_history$events,1)[[1]] #returns a list -- the first item of which is a test
                temp1[1]<-ifelse(is.element(test.outcome$type,c("iFOBT")),1,0)
                temp1[2]<-ifelse(is.element(test.outcome$compliance,c("accept")),1,0)
                temp1[3]<-ifelse(is.element(test.outcome$type,c("blood")),1,0)
                temp1[4]<-ifelse(is.element(test.outcome$compliance,c("screen")),1,0)
                
                if(test.outcome$result=="positive"){  #if it is a false positive, then we will skip over all the
                    temp1[5]<-1 #test.result=="positive"
                    temp1[6]<-1 #person has a colonoscopy woth probability 1
                    temp1[13]<-sample(c(0,1),1,prob=c(0.9997,0.0003)) #probability of bleeding
                    temp1[14]<-sample(c(0,1),1,prob=c(0.9999,0.0001)) #probability of perforation
                    if ( (person$colon$state=="adenoma") | (person$colon$state=="large adenoma")){   #this may be wrong. 
                        temp1[7]<-1
                        temp1[12]<-1
                        person$in_treatment_program<-"yes"
                    } else if (person$colon$state=="CRC"){   #has this been changed to just "CRC" ??  yes **Changed
                        if (person$colon$stage=="A"){
                            temp1[8]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="B"){
                            temp1[9]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if (person$colon$stage=="C"){
                            temp1[10]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                        }
                        if(person$colon$stage=="D"){
                            temp1[11]<-1
                            temp1[12]<-1
                            person$in_treatment_program<-"yes"
                                        #we do nothing.
                        }
                    }
                } #end test.result=="positive"
                } #end person$colon==clear
                    temp1
               }

        ) #end method list
    )



# Duke's Adenoma Parameters class
#################################
#' @export DukesAdenomaParams
#' @exportClass DukesCrcSpinModel
#' 
#' @family DukesCrcSpinModel_classes
DukesAdenomaParams <- setRefClass( "DukesAdenomaParams",
                                  contains = "AdenomaParams",
    fields = list(
        sojournA_min = "numeric",
        sojournA_max = "numeric",
        sojournB_val = "numeric",
        sojournC_min = "numeric",
        sojournC_max = "numeric",
        sojournD_val = "numeric"
    ),

    methods = list (

        setParams = function(...){
            # Set all the parameters passed in
            # NOTE: by default this handles inherited
            # variables from AdenomaParams also
            initFields(...)

            # Call AdenomaParams$setParams() with no arguments
            # to ensure inherited variables are capped within
            # valid ranges
            callSuper()

            # Cap Duke's specific variables to valid ranges
            sojournA_max <<- max(sojournA_min, sojournA_max)
            sojournC_max <<- max(sojournC_min, sojournC_max)
        },

        initialize = function(){

            # Set Duke's specific defaults here
            .self$sojournA_min <<- 1
            .self$sojournA_max <<- 3
            .self$sojournB_val <<- 1
            .self$sojournC_min <<- 1
            .self$sojournC_max <<- 2
            .self$sojournD_val <<- 1

            # Call AdenomaParams superclass initializer to:
            #   - set defaults of inherited members
            #   - call setParams with no arguments to trigger
            #     valid range checks
            # NOTE: order is important here!
            #   callSuper MUST be called AFTER duke's variables
            #   are set, otherwise setParams validity checks will
            #   triggered prior to those variables receiving values.
            callSuper()

        }
    )
)


# Duke's Adenoma class
######################
#' @export DukesAdenoma
#' @exportClass DukesCrcSpinModel
#' 
#' @family DukesCrcSpinModel_classes
# Reference Class
DukesAdenoma <- setRefClass( "DukesAdenoma",
                contains="Adenoma",

    fields = list(
        adenoma_model_params               = "DukesAdenomaParams",
        stage                              = "character",
        transition_to_crc_year             = "numeric",
        crc_survival_time                  = "numeric",
        sojourn_time_A                     = "numeric",
        sojourn_time_B                     = "numeric",
        sojourn_time_C                     = "numeric",
        sojourn_time_D                     = "numeric",
        transition_to_stage_A              = "numeric",
        transition_to_stage_B              = "numeric",
        transition_to_stage_C              = "numeric",
        transition_to_stage_D              = "numeric"
        ),

    methods = list(

        initialize = function (adenoma_params=NULL,
                        stage="", ...){

            if(is.null(adenoma_params)){
                # no adenoma_model_parameters were passed in
                # create one with default values
                adenoma_params<-DukesAdenomaParams$new()
            }

            # Call the Adenoma superclass initialiser to do
            # the bulk of the work
            callSuper(adenoma_params=adenoma_params, ...)

            # Set necessary Duke's specific variables
            initFields(stage=stage, adenoma_model_params=adenoma_params)
        },

        transition = function(subject_age=0){
            # The function handles ADENOMA transitions
            #
            # it calls its super class implementation to handle
            # adenoma to large adenoma
            # large adenoma to CRC (no stage set)
            #
            # Then handles transition of CRC stages
            # A to B
            # B to C
            # C  to D
            # D to deceased
            # A,B,C,D, are all "CRC"

            # The transition of COLON to "symptomatic CRC" is done
            # in the function DukesColon::modelTransitionToClinicalCRC()

            # Call Adenoma superclass implementation to
            # handle Classic CRC transitions
            callSuper(subject_age=subject_age)

            if (state =="CRC"){
                if(transition_to_preclinical_crc_year==subject_age){

                    stage<<-"A"
                    transition_to_stage_A<<-subject_age

                    sojourn_time_A<<-runif(1,adenoma_model_params$sojournA_min,adenoma_model_params$sojournA_max)
                    sojourn_time_B<<-adenoma_model_params$sojournB_val
                    sojourn_time_C<<-runif(1,adenoma_model_params$sojournC_min,adenoma_model_params$sojournC_max)
                    sojourn_time_D<<-adenoma_model_params$sojournD_val

                    # next line overwrites the sojourn time saved by super/parent class
                    sojourn_time<<-sojourn_time_A+sojourn_time_B+sojourn_time_C+sojourn_time_D
                }
              
                if(stage =="A"){
                  if (transition_to_stage_A+sojourn_time_A  <= subject_age){
                    transition_to_stage_B<<-subject_age
                    stage<<-"B"
                  }
                }
                if(stage =="B"){
                    if (transition_to_stage_B+sojourn_time_B  <= subject_age){
                      transition_to_stage_C<<-subject_age
                      stage<<-"C"
                    }
                }
                if(stage =="C"){
                      if (transition_to_stage_C+sojourn_time_C  <= subject_age){
                        transition_to_stage_D<<-subject_age
                        stage<<-"D"
                      }
                }
                if(stage =="D"){
                    if (transition_to_stage_D+sojourn_time_D  <= subject_age){
                        state<<-"deceased"
                        stage<<-"deceased"
                    }
                }
            }
        },

        checkCRCSymptomPresentation = function(subject_age=0){
            # In a Duke's CRC Spin model, presentation of symptoms
            # is modeled as a function of the colon state/stage
            # and does not involve a direct check of adenomas
            #
            # To avoid accidental errors/bugs, this dummy function
            # explicitly overwrites that inherited from the Classic
            # CRC spin models Adenoma Class, which DOES use a direct
            # check on adenomas
            return(FALSE)
        }

    )
)


# Duke's Colon Class
####################
#' @export DukesColon
#' @exportClass DukesColon
#' 
#' @family DukesCrcSpinModel_classes
# Reference Class
DukesColon <- setRefClass( "DukesColon",
                contains="Colon",

    fields = list(
        stage="character"
        ),

    methods = list(

        initialize = function (...) {

            # Call the Colon superclass initializer to
            # set colon generic variables
            callSuper(...)

            # Set DukesColon specific variables
            initFields(stage="clear")

        },

        getAdenomaType = function () {
            return(getRefClass("DukesAdenoma"))
        },

        treatCRC = function () {
            #This is a very simple treatment model. If Dukes A or B then all adenomas and tumors are removed.
            #If Dukes C or D then treatment is ineffective and the cancer in incurable.

            if ( (state=="adenoma") | (state=="large adenoma")){
              #I dont think this ever happens --
              #if state=="adenoma" "large adenoma then colon_clinical will not be CRC
              sites <<- vector("list",length=0)
              state<<-"clear"
              stage<<-"clear"
            } else if (state=="CRC"){
              if (stage=="A"){
                sites <<- vector("list",length=0)
                state<<-"clear"
                stage<<-"clear"
              }
              if (stage=="B"){
                sites <<- vector("list",length=0)
                state<<-"clear"
                stage<<-"clear"
              }
              if (stage=="C"){
                # we do nothing
              }
              if(stage=="D"){
                #we do nothing.
              }
            }
        },

        modelTransitionToClinicalCRC = function (subject_age=0) {
            if (state=="CRC"){
                p<-c(NA,NA)
                if (stage=="A"){
                    p=c(0.04,1-0.04)
                } else if (stage=="B"){
                    p=c(0.18,1-0.18)
                } else if (stage=="C"){
                    p=c(0.37,1-0.37)
                } else if (stage=="D"){
                    p=c(0.74,1-0.74)
                }

                if(!is.na(p[1])){
                    if (sample(c("symptoms","no symptoms"),1,prob=p)=="symptoms"){
                        state<<-"symptomatic CRC"
                        return(TRUE)
                    }
                }
            }
            return(FALSE)
        },

        getMostAdvancedAdenomaState = function () {

            # Call Colon superclass implementation to do
            # Classic CRC spin checks
            adv_state<-callSuper()

            # Now do additional Duke's updates
            if (length(sites)>0){

                tt<-lapply(sites,f<-function(x){x$state})
                tt<-unlist(tt)
                if (  sum(tt=="deceased")>0 ) {
                    adv_state<-"deceased"
                }
            }

            return(adv_state)

        },

        getMostAdvancedAdenomaDukesStage = function() {

            if (length(sites)>0){
                tt<-lapply(sites,f<-function(x){x$stage})
                tt<-unlist(tt) # c("A","B","C","D","deceased")
                if (  sum(tt=="deceased")>0 ) {
                adv_stage<-"deceased"
                }else if (  sum(tt=="D")>0 ) {
                adv_stage<-"D"
                } else if ( sum(tt=="C")>0 ) {
                adv_stage<-"C"
                } else if ( sum(tt=="B")>0 ) {
                adv_stage<-"B"
                } else if ( sum(tt=="A")>0 ) {
                adv_stage<-"A"
                } else {
                adv_stage<-"clear"
                }

            } else {
              adv_stage<-"clear"
            }

            return(adv_stage)
        },

        updateState = function () {

            # Call the Colon superclass implementation
            # to set state and cancer_site
            callSuper()

            # Set the Duke's stage
            stage<<-getMostAdvancedAdenomaDukesStage()
        }
    )
)

# Duke's Person with colon class
################################
#' @export DukesPersonWithColon
#' @exportClass DukesPersonWithColon
#' 
#' @family DukesCrcSpinModel_classes
# Reference Class
DukesPersonWithColon <- setRefClass( "DukesPersonWithColon",
            contains="PersonWithColon",

    methods = list(

        initialize = function(...){
            # call PersonWithColon superclass initializer
            # to initialize generic PersonWithColon values
            # BUT set colon argument to a new DukesColon
            # to override/replace the standard Colon object
            # otherwise created in the superclass initializer
            callSuper(colon=DukesColon$new(), ...)
        },

        modelDeathFromCRC = function () {
            if (colon$state=="deceased"){
                state<<-"deceased CRC"
            }
        },

        isDead = function () {
            return(state=="deceased" | state=="deceased CRC")
        },

        updateState = function () {

            # call the PersonWithColon superclass implementation
            # to do Classic CRC updates
            callSuper()

            # Do updates specific to Duke's Model
            modelDeathFromCRC()
        },

        modelPreClinicalToClincalTransition = function () {

            if(colon_clinical!="symptomatic CRC"){
                # Call our PersonWithColon superclass
                # implementation to do the transition
                callSuper()
     
                # If a transition occurred, add a clinical event
                if(colon_clinical=="symptomatic CRC"){
                    clinical_history$events<<-c(clinical_history$events,
                                SymptomaticPresentation$new(age=age,
                                                            cancer_stage=colon$stage))
                }
            }
        },

        treatmentRecordSize = function () {
            return(28)
        },
        testForAndTreatSize = function () {
            return(14)
        },
       NBCSPRecordSize = function () {
            return(14)
        },

        initiateCRCTreatment = function (prob_colonoscopy_performed=1) {
            # model chance of having a colonoscopy
            #   Originally:
            #       temp1[6]<-1
            colonoscopy_performed<-
                    sample(
                        c(FALSE,TRUE),
                        1,
                        prob=c(1-prob_colonoscopy_performed,
                               prob_colonoscopy_performed))

            colonoscopy_caused_bleeding<-FALSE
            colonoscopy_caused_perforation<-FALSE
            if (colonoscopy_performed==TRUE) {
                # model chance of complications
                #   Originally:
                #       temp1[13]<-sample(c(0,1),1,prob=c(0.9997,0.0003))
                #       temp1[14]<-sample(c(0,1),1,prob=c(0.9999,0.0001))
                colonoscopy_caused_bleeding<-
                    sample(
                        c(FALSE,TRUE),
                        1,
                        prob=c(0.9997,0.0003))
                colonoscopy_caused_perforation<-
                    sample(
                        c(FALSE,TRUE),
                        1,
                        prob=c(0.9999,0.0001))
            }

            # Collect more information for the pre-treatment record/chart
            #   Originally:
            #       this code was spread around as temp[x]=1 throughout code
            colon_has_small_to_large_adenomas<-
                (colon$state=="adenoma" | colon$state=="large adenoma")
            colon_state_CRC<-(colon$state=="symptomatic CRC"||colon$state=="CRC")
            colon_in_pre_symtomatic_stage_A<-(colon_state_CRC&&colon$stage=="A")
            colon_in_pre_symtomatic_stage_B<-(colon_state_CRC&&colon$stage=="B")
            colon_in_pre_symtomatic_stage_C<-(colon_state_CRC&&colon$stage=="C")
            colon_in_pre_symtomatic_stage_D<-(colon_state_CRC&&colon$stage=="D")


            # model treatment of CRC
            #   Originally:
            #       #This is a very simple treatment model. If Dukes A or B then all adenomas and tumors are removed.
            #       #If Dukes C or D then treatment is ineffective and the cancer in incurable.
            #       if ( (person1@colon@state=="adenoma") | (person1@colon@state=="large adenoma")){
            #           ...etc...
#            colon$treatCRC()

            if(length(colon$sites)==0 &&
               colon$stage=="clear" &&
               colon$state=="clear"){
                colon_clinical<<-"clear"
            }

            in_treatment_program<<-"yes"

            return(c(FALSE, # offered iFOBT
                     FALSE, # accepted iFOBT
                     FALSE, # offered blood test
                     FALSE, # accepted blood test
                     FALSE, # test result positive
                     colonoscopy_performed,
                     colon_has_small_to_large_adenomas,
                     colon_in_pre_symtomatic_stage_A,
                     colon_in_pre_symtomatic_stage_B,
                     colon_in_pre_symtomatic_stage_C,
                     colon_in_pre_symtomatic_stage_D,
                     (in_treatment_program=="yes"),
                     colonoscopy_caused_bleeding,
                     colonoscopy_caused_perforation))

        },

        medicalSnapshotSize = function () {
            return(46)
        },

        getMedicalSnapshot = function () {
            # This function collects a tally of various study parameters
            # related to a DukesPersonWithColon
            #
            # It calls the PersonWithColon superclass's implementation of this
            # function to retrieve some parameters, then augments these results

            aa<-rep(FALSE,46)


            ##############################################################
            # Copy the stuff we want from our super class results into the
            # location we want them

            # Call the PersonWithColon superclass implementation
            # to do some of the work
            res<-callSuper()

            # Copy over COUNTS of:
            #   adenomas, large adenomas, & CRC adenomas
            aa[1:3]<-res[1:3]

            # Copy over FLAGS for:
            #   colon_clinical=="symptomatic CRC"
            #   colon$state=="clear", =="adenoma", =="large adenoma",
            #              =="CRC", =="symptomatic CRC"
            aa[5:10]<-res[4:9]

            # Copy FLAGS indicating location of colon$cancer_site
            aa[34:39]<-res[10:15]

            # Copy FLAG for colon_clinical=="symptomatic CRC"
            aa[46]<-res[17]



            ##############################################################
            # Add new information important to Duke's CRCSpin Model

            # Add a count of deceased adenomas
            if (length(colon$sites)>0){
                tt<-lapply(colon$sites,f<-function(x){(x$state=="deceased")})
                tt<-unlist(tt)
                aa[4]<-sum(tt)
            }

            # Add flag for colon$state=="deceased"
            aa[11]<-(colon$state=="deceased")

            # Add in some additional info for when colon$state is 
            # symptomatic CRC or CRC
            if (length(colon$sites)>0
                & (colon$state=="CRC" | colon$state=="symptomatic CRC")
                ){

                # rather than duplicating value counting code
                # within the following two if statements (which
                # never execute together) to operate on separate
                # segments of the output vector, we simply use an
                # offset to the start of the effected range to achieve
                # the same results.
                start_range<-12
                if(colon$state=="CRC"){
                  start_range<-12
                } else if (colon$state=="symptomatic CRC"){
                  start_range<-23
                }

                # Add COUNTS of adenomas in stages A,B,C,D,or deceased
                tt<-lapply(colon$sites,f<-function(x){x$stage})
                tt<-unlist(tt)
                tt<-factor(tt, levels=c("A","B","C", "D","deceased"))
                aa[(0:4)+start_range]<-table(tt)

                # Add FLAGS indication location of colon$cancer_site
                tt<-factor(colon$cancer_site,
                        levels=c("cecum","ascending","transverse","descending","sigmoid","rectum"))
                aa[(5:10)+start_range]<-table(tt)

            }

            # Add FLAGS indicating colon$stage
            tt<-factor(colon$stage, levels=c("A","B","C", "D","deceased"))
            aa[40:44]<-table(tt)

            # Add FLAG indicating if person is dead
            aa[45]<-ifelse(is.element(state,c("deceased","deceased CRC")),1,0)

            return(aa)

        }
    )
)

# Misc
######

# Reference Class
SymptomaticPresentation <- setRefClass( "SymptomaticPresentation",
            fields = list(
                age="numeric",
                cancer_stage="character"
                )
)
