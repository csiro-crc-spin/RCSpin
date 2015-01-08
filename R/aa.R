# vim: set filetype=r fileformat=unix : #

############################
#
# Generic spin model classes
#
############################


# the following variables should not be global...but i got lazy

#http://www.abs.gov.au/AUSSTATS/abs@.nsf/DetailsPage/3302.0.55.0012008-2010?OpenDocument
#
#lx the number of persons surviving to exact age x;
#qx the proportion of persons dying between exact age x and exact age x+1
#It is the mortality rate, from which other functions of the life table are derived;
#Lx the number of person years lived within the age interval x to x+1; and
#ex life expectancy at exact age x.

#dr_data<-read.table(file="../data/sources/3302055001DO001_20082010.csv",sep=",",skip=5,header=TRUE,as.is=TRUE)
#dr_data<-dr_data[-c(1,103,104,105),]
#death_rate_male<-as.numeric(dr_data$qx)
#death_rate_female<-as.numeric(dr_data$qx.1)

# Model class
#############
#' Class implementing a simple natural aging model on a group of people
#'
#' This class implements a natural aging model on a group of people
#' represented as a \code{study_group} of \code{\link{Person}} class
#' objects.
#'
#' The main aim of this Class is to act as a framework for building more 
#' complex models through extension of this class. The simple natural
#' aging model implemented in \code{\link{updateSubject}} acts as an 
#' example and place holder for subclasses to override with a more
#' complex implementation of the function
#'
#' @field study_results A two dimensional matrix containing summary of model
#'        state at each year of the study. Each column describes a particular
#'        property or tally of the study population's state. Results for
#'        females and males of the population are recorded separately in odd
#'        and even numbered row respectively. So the number of rows is 2x the
#'        number of iterations (years).
#' 
#' @export GenericModel
#' @exportClass GenericModel
#' 
#' @family GenericModel_classes
GenericModel <- setRefClass( "GenericModel",

    fields = list(
        iterations="numeric",
        iteration_resolution="numeric",
        study_group="list",
        study_results="matrix",
        commencement_age="numeric"
    ),

    methods = list(

        initialize = function(study_group=list(),
                        iterations=0,
                        iteration_resolution=1,
                        num_subjects=0,
                        base_seed=NA,
                        commencement_age=1,
                        study_results=NA,
                        set.sex=NA_character_,
                        ...){

        "Create and initialize a new instance of a GenericModel
        
         This function creates and initializes a new instance of the
         GenericModel class. Parameters relating to the length of the
         modelled study and the number of subjects in the study group
         are passed to this function. See \\code{\\link{run}} to run the
         created model.
        
         @section Model details:
        
         The GenericModel classes present a framework or set of base
         classes on which to build more complex models. As a place holder
         the GenericModel implements a simple natural aging model, and
         and records the number of non-deceased subjects remaining at each
         iteration.
        
         @param study_group list containing objects of a class derived from 
                \\code{\\link{Person}}. Used to specify an existing group you
                wish to run/continue the model on

         @param iterations the number of iterations the model should run to

         @param iteration_resolution how many iterations are modelled in one
                call to the subject model step function \\code{\\link{updateSubject}}

         @param num_subjects integer indicating number of subjects in the study
                group. If \\code{study_group} is not specified, a new random \\code{study_group}
                of size \\code{sum_subjects}, containing objects of class \\code{\\link{Person}},
                is created. Otherwise, it is taken as the size of the specified
                \\code{study_group} and may result in \\code{study_group} being
                truncated

         @param base_seed integer RNG seed used for model

         @param commencement_age integer indicating age at which to start
                running the model of subjects of the \\code{study_group}

         @param study_results an integer matrix with \\code{iterations*2} rows
                and enough columns to fit the model's iteration results.
                Used to specify existing results when an existing \\code{study_group}
                is specified.

         @param ... additional model values/parameters

         @return a new object of type GenericModel with specified configuration

         @seealso \\code{\\link{run}} to run the model

         @family model_initializers

         @callGraph

         @author Luke Domanski \\email{luke.domanski@@csiro.au}

         @aliases GenericModel$new
         "
            if (class(study_results)!="matrix" &&
                length(study_results)==1 &&
                is.na(study_results)) {
                study_res<-matrix(0,iterations*2,.self$getModelResultSize())
            } else if (dim(study_results)[1] < iterations*2){
                # Handle the situation where provided
                # study_results matrix does not have enough rows
                # to hold an aggregated population result for each
                # iteration/step of the study

                # If you do not want this granularity of study result
                # (e.g. you want a single aggregated across all
                # subjects AND iterations/steps), this code and
                # run()/doIteration() will need to be overwritten in
                # a derived class
            } else {
                study_res<-study_results
            }

            initFields(iterations=iterations,
                iteration_resolution=iteration_resolution,
                study_group=study_group,
                study_results=study_res,
                commencement_age=commencement_age)

            if (length(study_group)==0) {
                if (num_subjects>0) {
                    for (i in 1:num_subjects){
                        .self$study_group <<-
                            c(.self$study_group,
                              Person$new(study_id=i,
                                         base_seed=base_seed,
                                         age=commencement_age,
                                         sex=set.sex)
                              )
                    }
                }
            } else {
                # checks & processing should be done here to:
                #   1) Check all elements of study_group are objects
                #      deriving from Person
                #   2) truncate study_group length to num_subjects
                .self$study_group <<- study_group

            }

        },

        run = function(){

        "Runs the model
        
         Runs the model for the number of iterations specified at model
         initialization, defined as \\code{num_iterations-commencement_age},
         and sets the male and female result rows for each iteration in
         the model's study_results matrix.
        
         The start iteration is \\code{commencement_age}, rows of the \\code{study_result}
         matrix prior to \\code{study_result[commencement_age]} are unchanged.
        
         @callGraph
         @author Luke Domanski \\email{luke.domanski@@csiro.au}
        "

            for (i in commencement_age:iterations){
                res<-doIteration()
                study_results[i*2-1,]<<-res[1,]
                study_results[i*2,]<<-res[2,]
            }

        },

        getModelResultSize = function () {

        "Returns the number of columns in the model's study result matrix
        
         @author Luke Domanski \\email{luke.domanski@@csiro.au}
        "
            return(1)
        },

        updateSubject = function (subject) {

        "Applies a model step to an individual subject
        
         This function applies a model step to an individual subject
         updating the subject's state. It implements the main, overall,
         structure and flow of the model.
        
         @param subject an object of class Person on which to apply
                 a model step

         @return an integer vector of length equal to the width of the
                 model's study_result matrix, containing the iteration
                 results for the subject

         @callGraph

         @author Luke Domanski \\email{luke.domanski@@csiro.au}
        "

            if (subject$state=="living") {
                subject$restoreRNGState()
                subject$age = subject$age + iteration_resolution
                subject$modelDeathFromOtherCauses()
                subject$saveRNGState()
            }
            return((subject$state=="living"))
        },

        doIteration = function(){

        "Performs one model iteration over all subjects in the study_group
        
         Performs one model iteration over all subjects in the model's
         study_group, sums the iteration results of all subjects into
         separate aggregated results for females and males, returns the
         aggregated results.
        
         @return A two row matrix the same width as the model's
                 study_result matrix. The first row contains aggregated
                 iteration results for females, the second, for males.

         @callGraph

         @author Luke Domanski \\email{luke.domanski@@csiro.au}
        "

            study_record<-matrix(0,2,dim(study_results)[2])
            for (i in study_group){
                sex_row<-1
                if(i$sex=='F'){
                    sex_row<-1
                } else {
                    sex_row<-2
                }
                study_record[sex_row,]<-study_record[sex_row,] + updateSubject(i)
            }
            return(study_record)
        },

        show = function(){

        "Display information about the object.
        "

            cat("Spin Model object of class",
                classLabel(class(.self)), "\n")
            cat("Iterations:")
            methods::show(iterations)
            cat("Iteration resolution:")
            methods::show(iteration_resolution)
            cat("Commencement age:")
            methods::show(commencement_age)
            cat("Study group:\n")
            #methods::show(study_group)
            cat("Study group size:")
            cat(length(study_group), "\n")
        }

    )
)


# Clinical Test classes
#######################

#' Class implementing records of clinical tests
#'
#' @field age        Age at which the test was performed.
#' @field type       String identifier of the test.
#' @field compliance String describing subjects acceptance of the test.
#' @field result     String describing the result of test.
#' @field state      String describing the state of health the person is in at
#'                   time of test.
#'  @export Test
#'  @exportClass Test
#' @family GenericModel_classes
Test <- setRefClass( "Test",
    fields = list(
        age="numeric",
        type="character",
        compliance="character",
        result="character",
        state="character"
    ),

    methods = list(

        summary = function(){},

        show = function(){

        "Display information about the object.
        "

            cat("Spin Test object of class",
                classLabel(class(.self)), "\n")
            cat("Age:")
            methods::show(age)
            cat("Type of test:")
            methods::show(type)
            cat("Compliance:")
            methods::show(compliance)
            cat("Result:")
            methods::show(result)
            cat("State:")
            methods::show(state)
        }

    )
)


# Clinical History class
########################
#' Class implementing records of clinical history
#'
#' This class is used to keep a list of test and events which occur through a
#' subject's clinical history. The records are stored in a list field
#' \code{events} which can contain an object of virtually any kind or class. It
#' is the responsibility of the class or function using this clinical history
#' class to filter out \code{event} objects that they do not understand how to
#' interpret.
#'
#' It is used by an associated \code{\link{Person}} class representing the
#' subjects studied in the \code{\link{GenericModel}} clinical model. However,
#' it aims is to act as generic event container that remains relevant across a
#' range of more complex clinical model derived from from those base classes.
#'
#' @field status .
#' @field events List of records on a subject's medical tests or events.
#'
#' @export ClinicalHistory
#' @exportClass ClinicalHistory
#' 
#' @family GenericModel_classes
ClinicalHistory <- setRefClass( "ClinicalHistory",

    fields = list(
        status="character",
        events="list"
        ),

    methods = list(

        summary = function(){},

        show = function(){

        "Display information about the object.
        "

            cat("Spin ClinicalHistory object of class",
                classLabel(class(.self)), "\n")
            cat("Status:")
            methods::show(status)
            cat("Events:")
            methods::show(events)
        }

        )
)


# Person class
##############
#' Class representing a person involved in a clinical study
#'
#' This class represents a person involved in a clinical study and
#' provides a mechanism to record the person's clinical history, their
#' involvement in a treatment program, and their life state.
#'
#' The class implements a function \code{\link{modelDeathFromOtherCauses}}
#' modeling general death of the person.
#'
#' @field age The current age of the person being modeled.
#' 
#' @field sex The gender of the person.
#' 
#' @field state String describing the state of the person e.g. "deceased",
#'        "living" etc. as required by dervived classes and their drivers 
#' 
#' @field in_treatment_program String "yes" or "no" etc. indicating whether the
#'        person is undergoing medical treatment that may alter the way they
#'        clinical progression is modeled.
#' 
#' @field clinical_history A \code{\link{ClinicalHistory}} object holding the
#'        persons medical record.
#'
#' @field study_id Study ID number assigned to the person
#'
#' @section Note:
#'
#' The main aim of this class is to act as a framework for building person
#' representations operated on by more complex extensions of the
#' \code{\link{GenericModel}} model driver class.
#'
#' @export Person
#' @exportClass Person
#' @family GenericModel_classes
Person <- setRefClass( "Person",

    fields = list(
        age="numeric",
        sex="character",
        state="character", #"deceased","living", etc (as needed by subclasses)
        in_treatment_program="character",
        clinical_history="ClinicalHistory",
        study_id="numeric",
        random_seed_state="integer"
    ),

    methods = list(

        initialize = function(study_id=0,
                            base_seed=NA,
                            age=1,
                            sex=NA_character_,
                            state="living",
                            in_treatment_program="no",
                            clinical_history=ClinicalHistory$new(),
                            random_seed_state=NA,
                            ...){

        "Create and initialize a new instance of a Person
        
         This function creates and initializes a new instance of the
         Person class. 
        
         @param age initial age of the person.

         @param sex the person's sex, if not specific, will be randomly chosen.

         @param in_treatment_program whether the person is currently in a treatment
                program.

         @param clinical_history the person's existing clinical history/record,
                if not specified a new empty/blank \\code{\\link{ClinicalHistory}}
                is created.

         @param study_id an integer assigned as their identifier within the study.

         @param base_seed integer RNG seed used for modeling

         @param random_seed_state the pseudo random number stream state to start
                if continuing from a previously saved study state 

         @param ... additional model values/parameters

         @return a new object of type Person with specified configuration

         @seealso \\code{\\link{GenericModel}} for usage

         @family model_initializers

         @callGraph

         @author Luke Domanski \\email{luke.domanski@@csiro.au}

         @aliases GenericModel$new
         "

            if(is.na(base_seed)){
                base_seed=125
            }

            if (is.na(random_seed_state)) {
                set.seed(base_seed+study_id)
                random_seed<-.Random.seed
            } else {
                random_seed<-random_seed_state
                .Random.seed<-random_seed
            }

            if (is.na(sex) |
                !(sex=="F" | sex=="M")) {
                s <- sample(c("F","M"),1,prob=c(0.5,0.5))
            } else {
                s<-sex
            }

            random_seed<-.Random.seed

            initFields(age=age,
                sex=s,
                state=state,
                in_treatment_program=in_treatment_program,
                clinical_history=clinical_history,
                study_id=study_id,
                random_seed_state=random_seed)

        },

        saveRNGState = function () {

        "Saves the pseudo random number stream state for the person"

            random_seed_state<<-.Random.seed
        },

        restoreRNGState = function () {

        "Saves the pseudo random number stream state for the person"

            .Random.seed<-random_seed_state
        },

        modelDeathFromOtherCauses = function () {

        "Models whether the person dies from natural/common causes
        
         Models whether the person dies from natural/common causes
         based on probabilities recorded in the \\code{death_rate_female}
         and \\code{death_rate_male} arrays upon package loading.
        
         @author Rob Dunne \\email{rob.dunne@@csiro.au}
                 Luke Domanski \\email{luke.domanski@@csiro.au}
        "
            if(sex=="F"){
                state<<-sample(
                            c("deceased","living"),1,
                            prob =c( death_rate_female[age], 1-death_rate_female[age]))
            }
            if(sex=="M"){
                state<<-sample(
                            c("deceased","living"),1,
                            prob =c( death_rate_male[age], 1-death_rate_male[age]))
            }
        },

        summary = function(){
            cat("age is",age,"\n")
            # cat("number of adenomas", length(colon$sites),"\n")
            # if (all.adenomas==TRUE){
                # ss<-colon$sites
                # if (length(ss)>0){
            # #    for (i in 1:length(ss)){
            # #      temp<-ss[[i]]
            # #      cat(temp@location, " ", temp@stage," ", "years ",study.year-temp@initiated.in.year," ","\n")
            # #    }
            # #  }
                    # temp<-lapply(ss,f<-function(x){x$size})
                    # cat("size of largest adenoma",  max(unlist(temp)),"\n")
                    # aa<-lapply(ss,f<-function(x){x$location})
                    # cat("locations of adenomas:",unlist(aa),"\n")
                    # aa<-lapply(ss,f<-function(x){x$state})
                    # cat("state of adenomas:",unlist(aa),"\n")
                # }
            # }
            # cat("colon state:",colon$state,"\n")
            # cat("colon stage:",colon$stage,"\n")
            # cat("clinical:",colon.clinical,"\n")
            # cat("person state:",state,"\n")

        },

        show = function(){

        "Display information about the object.
        "

            cat("Spin Person object of class",
                classLabel(class(.self)), "\n")
            cat("Age:")
            methods::show(age)
            cat("Sex:")
            methods::show(sex)
            cat("State:")
            methods::show(state)
            cat("In treatment program:")
            methods::show(in_treatment_program)
            cat("Study id:")
            methods::show(study_id)
#            cat("System data:")
 #           methods::show(system_data)
            cat("Clinical history:\n")
            clinical_history$show()
        }

        )
)



# ancillory function
#' append to a list
#'
#' This function appends on item to a list
#'
#' @param my.list a list
#' 
#' @param obj an item to append to the list
#' 
#' @export
#' 
#' @examples
#' my.list<-list("a","b","c")
#' lappend(my.list,"d")
lappend <- function(my.list, obj) {
    my.list[[length(my.list)+1]] <- obj
    return(my.list)
}
# vim: set filetype=r fileformat=unix : #
#source("genericSpin.R")


########################
#
# CRC spin model classes
#
########################

# forward definitions
CrcRiskParams <- setRefClass( "CrcRiskParams")
AdenomaParams <- setRefClass( "AdenomaParams")

# CRC model class
#################
#' Class orchestrating the CRC Spin model
#'
#' This class acts as the entry point and driver for an implementation of the
#' CRC Spin model, described by X et al., which is implemented using a
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
#' @export CrcSpinModel
#' @exportClass CrcSpinModel
#' 
#' @family CrcSpinModel_classes
CrcSpinModel <- setRefClass( "CrcSpinModel",
            contains="GenericModel",

    fields = list(
        crcrisk_model_params = "CrcRiskParams",
        adenoma_model_params = "AdenomaParams"
    ),

    methods = list(


#################################################################
# This block of functions is not part of the model.
# Although the initialize function does provide an interface to
# create a model instance and specify important parameters.
# They exist mainly to provide a generic and extensible Spin style
# model framework to build off, without the need to replicate this
# functionality in derived classes.
        initialize = function (study_group=list(),
                        iterations=0,
                        num_subjects=1,
                        base_seed=NA,
                        commencement_age=20,
                        set.sex=NA_character_,
                        ...) {

            crcrisk_model_params<<-crcRiskParamsType()$new()
            adenoma_model_params<<-adenomaParamsType()$new()

            temp_study_group<-list()

            if (num_subjects<1) {
                cat("WARNING: CrcSpinModel$new(): num_subjects<1, Number of subjects requested must be greater than 1: Forcing a study group of 1 person\n")
                num_subjects<-1
            }

            if (length(study_group)==0) {
                for (i in 1:num_subjects){
                    temp_study_group <-
                        c(temp_study_group,
                          personWithColonType()$new(study_id=i,
                                        base_seed=base_seed,
                                        crcrisk_params=crcrisk_model_params,
                                        age=commencement_age,
                                        sex=set.sex ))
                }
            } else {
                # checks & processing should be done here to:
                #   1) Check all elements of study_group are objects
                #      deriving from PersonWithColon
                #   2) truncate study_group length to num_subjects
                # TODO...
                temp_study_group<-study_group
            }

            callSuper(iterations=iterations,
                num_subjects=num_subjects,
                study_results=matrix(0,iterations*2,
                                    (temp_study_group[[1]]$treatmentRecordSize()+
                                     temp_study_group[[1]]$medicalSnapshotSize())
                              ),
                base_seed=base_seed,
                commencement_age=commencement_age,
                study_group=temp_study_group,
                ...)

        },

        # Overwriting these *Type() functions in subclasses
        # allows the CrcSpinModel$new()/initialize()
        # function to act like a templated function in
        # other languages. It will cause the initializer of
        # the specified classes to be called during CrcSpinModel
        # initialization, and they MUST be subclasses of
        # CrcRiskParams, AdenomaParams, and PersonWithColon
        # reference classes providing a superset of their
        # initialize functions interface
        crcRiskParamsType = function () {
            return(getRefClass("CrcRiskParams"))
        },

        adenomaParamsType = function () {
            return(getRefClass("AdenomaParams"))
        },

        personWithColonType = function () {
            return(getRefClass("PersonWithColon"))
        },
#################################################################



#################################################################
# Functions that implement the actual model start here

        subjectHasNotLeftStudy = function (person) {
            if(!person$isDead()
               && person$in_treatment_program=="no"){
                return(TRUE)
            } else {
                return(FALSE)
            }
        },

        modelSubjectDiseaseDevelopment = function (person){
                # model growth of adenomas
                #   Originally:
                #       person1<-initiate.adenomas(person1,study.year=study.year)
                #       person1<-update.adenoma.growth.curve(person1)
                #
                #   Alternate naming?:
                #       growAdenomas()
                person$colon$modelAdenomaGrowth(
                                        adenoma_model_params=adenoma_model_params,
                                        subject_age=person$age)

                # transition adenomas through different classification
                # states owing to their growth
                #   Originally:
                #       person1<-update.transitions(person1)
                #
                #   Alternate naming?:
                #       transitionAdenomas()
                person$colon$modelAdenomaTransitions(
                                            subject_age=person$age)

                # update the patient's and their colon's state based
                # on adenoma growth and classification
                #   Originally:
                #       person1<-get.patient.state(person1)
                person$updateState()

                # If the person has not already died
                #   Note:
                #       In the Classic CRC Spin model a
                #       person can never actually die while
                #       growing and transitioning adenomas, or
                #       updating a person's state. But
                #       including this check allows reuse in
                #       more complex models
                if(!person$isDead()){
                    # model transition of a person and their colon to 
                    # full CRC
                    #   Originally:
                    #       person1<-update.transitions.pre.clinical.to.clinical(person1)
                    #
                    #   Alternate naming?:
                    #       checkCRCSymptoms()
                    person$modelPreClinicalToClinicalTransition()
                }
        },

        testForAndTreatCRC = function (person) {

            treatment_record=rep(FALSE,person$testForAndTreatSize())

            # if person has symptoms of CRC but is not in treatment
            # program
            #   Originally:
            #       if ((person1@colon.clinical=="CRC")
            #           & (person1@in.treatment.program=="no"))
            if (person$showingCRCSymptoms()
                && person$in_treatment_program=="no"){

                # perform initial treatments if necessary
                #   Originally:
                #       temp1[12]<-1  #people entering in.treatment.program
                #       person1@in.treatment.program<-"yes"
                treatment_record<-person$initiateCRCTreatment()

            }

            return(treatment_record)

        },

        modelSubjectTesting = function (person) {

            # test for and treat CRC
            #   Call function containing what was originally:
            #       if ((person1@colon.clinical=="CRC")
            #           & (person1@in.treatment.program=="no")){
            #               ...some code...
            #       }
            treatment_record <- testForAndTreatCRC(person)

            return(treatment_record)

        },

        #' Returns the number of columns in the model's study result matrix
        #'
        #' @author Luke Domanski \email{luke.domanski@@csiro.au}
        getModelResultSize = function () {
            temp<-personWithColonType()$new()
            return(temp$treatmentRecordSize()+temp$medicalSnapshotSize())
        },

        updateSubject = function (person) {
            # This function defines the main/overall model sequence/progression

            treatment_record=rep(FALSE,person$treatmentRecordSize())
            medical_snapshot=rep(FALSE,person$medicalSnapshotSize())

            # Don't model persons progress if they have left
            # the study
            #   Originally:
            #       if(is.element(person1@state,c("deceased"))){
            #         break
            #       }
            #       if(person1@in.treatment.program=="yes"){
            #           break
            #       }
            if(subjectHasNotLeftStudy(person)){

                person$restoreRNGState()

                # age the person
                #   Originally:
                #       study.year<-i
                #       person1@age<-i
                person$age<-person$age + iteration_resolution

                # model disease in a person through various 
                # stages
                #   Call function containing what was originally:
                #       person1<-initiate.adenomas(person1,study.year=study.year)
                #       person1<-update.adenoma.growth.curve(person1)
                #       person1<-update.transitions(person1)
                #       person1<-get.patient.state(person1)
                #       person1<-update.transitions.pre.clinical.to.clinical(person1)
                modelSubjectDiseaseDevelopment(person)

                # If person has not died prior to CRC being
                # detected and treated
                #   Note:
                #       In the Classic CRC Spin model a
                #       person can never actually die during
                #       modelSubjectDiseaseDevelopment. But
                #       including this check allows reuse in
                #       more complex models
                if(!person$isDead()){

                    # if the person is not already in a treatment program
                    if (person$in_treatment_program=="no") {
                        treatment_record <- modelSubjectTesting(person)
                    }

                    # model natural/other death
                    #   Originally:
                    #       person1<-death.rate.other.causes(person1)
                    #
                    #   Alternate naming?
                    #       checkDeathFromOtherCauses()
                    person$modelDeathFromOtherCauses()

                }

                # get a snapshot of the persons health at
                # this point in time
                #   Originally:
                #       Obtained by call to result.person.4(person)
                medical_snapshot<-person$getMedicalSnapshot()

                person$saveRNGState()

            }

            return(c(medical_snapshot, treatment_record))
            # Originally:
            #   bb[i,1:17]<-result.person.4(person)
            #   bb[i,18:31]<-temp1
        },

        set_adenoma_modeling_parameters = function(...){
            adenoma_model_params$setParams(...)
        },

        propegate_model_parameters = function(){
            for (i in study_group){
                i$colon$risk$init(subject_sex=i$sex)
            }
        },

        set_crcrisk_modeling_parameters = function(...){
            crcrisk_model_params$setParams(...)
            propegate_model_parameters()
        }
    )
)

# CRC risk parameters class
###########################
#' Class that hold parameters used to model a persons risk factors
#'
#' This class contains a number of parameters utilised by \code{\link{CrCRisk}}
#' to model/generate a person's risk factors, and by other classes in the
#' \code{\link{CrcSpinModel}} family.
#'
#' Upon object initialisation via \code{$new()} the parameters take on hard
#' default values (see \code{\link{initialize}} which can later be changed via
#' the \code{\link{setParams}} method.
#'
#' Storing the risk factor modeling parameters independent to the risk factor
#' class allows driver classes to set population wide modeling parameters
#' without needing to interface the risk factor class of every individual
#' directly, and to pass the parameters to other classes if required.
#'
#' @field baseline_mean Mean used to draw base line risk factor.
#'
#' @field baseline_sd Standard deviation used to draw base line risk factor.
#'
#' @field sex_linked_mean Mean used to draw gender linked risk factor.
#'
#' @field sex_linked_sd Standard deviation used to draw gender linked risk factor.
#'
#' @field age_mean Mean used to draw age risk factor.
#'
#' @field age_sd Standard deviation used to draw age risk factor.
#'
#' @export CrcRiskParams
#' @exportClass CrcRiskParams
#' 
#' @family CrcSpinModel_classes
CrcRiskParams <- setRefClass( "CrcRiskParams",

    fields = list(
        sex_linked_mean = "numeric",
        sex_linked_sd   = "numeric",
        baseline_mean   = "numeric",
        baseline_sd     = "numeric",
        age_mean        = "numeric",
        age_sd          = "numeric"
    ),

    methods = list(

        setParams = function(...){
            initFields(...)
        },

        initialize = function(){

            .self$sex_linked_mean <<- -0.3
            .self$sex_linked_sd   <<- 0.04
            .self$baseline_mean   <<- -6.7
            .self$baseline_sd     <<- 0.27
            .self$age_mean        <<- 0.03
            .self$age_sd          <<- 0.003

        }
    )
)


# Risk class
###########
#' Class representing a person's risk from a number of factors
#'
#' This class represents a persons risk in some given context. It aims to be a
#' generic container of risk variable that can be reused in various models.
#'
#' @field baseline_risk
#'
#' @field sex_linked_risk Additional risk factor based on a persons gender
#'
#' @field age_risk Additional risk factor based on a persons current age
#'
#' @export Risk
#' @exportClass Risk
#' 
#' @family CrcSpinModel_classes
Risk <- setRefClass( "Risk",
            fields = list(
                baseline_risk="numeric",
                sex_linked_risk="numeric",
                age_risk="numeric"
                )
)


# CRC risk class
################
#' Class representing a person's risk of developing adenomas that lead to CRC
#'
#' This class stores a person's CRC risk factors for use in the
#' \code{\link{CrcSpinModel}}, as well as a reference to the initial risk
#' modeling parameters used to generate the person's risk factors
#'
#' In particular, it gives context to risk factors defined in \code{\link{Risk}}
#' (from which it derives) within the \code{\link{CrcSpinModel}} for a
#' \code{\link{PersonWithColon}} where it represents the risk of them
#' developing \code{\link{Adenoma}}s that lead to CRC developing.
#'
#' @field risk_level String "high"|"standard" that was specified for the
#'        associated person at initialisation
#'
#' @field crcrisk_model_params Reference to the CrcRiskParams object specified
#'        at initialisation and used to seed/generate the risk factor values
#'
#' @field baseline_risk Inherited from \code{\link{Risk}}
#'
#' @field sex_linked_risk Inherited from \code{\link{Risk}}
#'
#' @field age_risk A 4-vector representing risk at ages prior to 50, 50-60,
#'        60-70, 70+ respectively. Overwrites \code{\link{Risk}.age_risk}
#'
#' @export CrcRisk
#' @exportClass CrcRisk
#' 
#' @family CrcSpinModel_classes
CrcRisk <- setRefClass( "CrcRisk",
            contains="Risk",

    fields = list(
        risk_level           = "character",
        crcrisk_model_params = "CrcRiskParams"
    ),

    methods = list(

        initialize = function(subject_sex=sample(c("F","M"),1,prob=c(0.5,0.5)),
                        subject_risk_level="standard",
                        crcrisk_params=NULL,
                        ...){

            if(is.null(crcrisk_params)){
                # no crcrisk_model_parameters were passed in
                # create one with default values
                crcrisk_params<-CrcRiskParams$new()
            }

            .self$risk_level           <<- subject_risk_level
            .self$crcrisk_model_params <<- crcrisk_params
            .self$init(subject_sex)
        },

        init = function(subject_sex){

            # Set sex linked risk based on subject sex
            sex_linked_risk <<- rnorm(1, mean=crcrisk_model_params$sex_linked_mean, sd=(crcrisk_model_params$sex_linked_sd))
            if(subject_sex=="M"){
              sex_linked_risk <<- -sex_linked_risk
            }

            # Set the basline risk according to subject_risk_level
            if (risk_level=="high"){
              baseline_risk <<- rnorm(1, mean=-4, sd=(0.27))
            } else {
              mu <- rnorm(1, mean=crcrisk_model_params$baseline_mean, sd=(crcrisk_model_params$baseline_sd))
              sigma <- runif(1, 0.1, 2.5)
              baseline_risk <<- rnorm(1, mean=mu, sd=sigma)
            }

            # Set age link risk
            age_risk <<- rnorm(4, mean=crcrisk_model_params$age_mean,sd=(crcrisk_model_params$age_sd))
            age_risk[1] <<- max(age_risk[1], 0)
            age_risk[2] <<- max(age_risk[2], 0)
            age_risk[3] <<- max(age_risk[3], 0)
            age_risk[4] <<- max(age_risk[4], 0)

        },

        summary = function(){},

        show = function(){

        "Display information about the object.
        "

            callSuper()
        }

    )
)


# Adenoma parameters class
##########################
#' Class that hold parameters used to derive factors used in model adenoma transitions
#'
#' This class contains a number of parameters utilised by \code{\link{Adenoma}}
#' to derive value for factors it uses to model the transition of Adenomas from
#' one state to the next (see \code{\link{Adenoma::transition}})
#'
#' Upon object initialisation via \code{$new()} the parameters take on hard
#' default values (see \code{\link{initialize}} which can later be changed via
#' the \code{\link{setParams}} method.
#'
#' Storing the parameters used to derive modeling values independent to the
#' \code{\link{Adenomas}} allows driver classes to set population wide modeling
#' parameters without needing to interface every Adenoma track in the model,
#' and to pass the parameters to other classes if required.
#'
#' @field gamma1_min
#'
#' @field gamma1_max
#'
#' @field gamma2_min
#'
#' @field gamma2_max
#'
#' @field gamma3_val
#'
#' @export AdenomaParams
#' @exportClass AdenomaParams
#' 
#' @family CrcSpinModel_classes
AdenomaParams <- setRefClass( "AdenomaParams",

    fields = list(
        beta1_min  = "numeric",
        beta1_max  = "numeric",
        beta2_min  = "numeric",
        beta2_max  = "numeric",
        gamma1_min = "numeric",
        gamma1_max = "numeric",
        gamma2_min = "numeric",
        gamma2_max = "numeric",
        gamma3_val = "numeric"
    ),

    methods = list(

        setParams = function(...){
            # Set all the parameters passed in
            initFields(...)

            # Cap certain variables to valid ranges
            beta1_max  <<- max(beta1_min, beta1_max)
            beta2_max  <<- max(beta2_min, beta2_max)
            gamma1_max <<- max(gamma1_min, gamma1_max)
            gamma2_max <<- max(gamma2_min, gamma2_max)
        },

        initialize = function(){

            # Set default values here
            .self$beta1_min  <<- 1
            .self$beta1_max  <<- 100
            .self$beta2_min  <<- 1
            .self$beta2_max  <<- 4
            .self$gamma1_min <<- 0.02
            .self$gamma1_max <<- 0.05
            .self$gamma2_min <<- 0
            .self$gamma2_max <<- 0.02
            .self$gamma3_val <<- 0.5

            # call setParams() with no arguments to ensure
            # your defaults above meet valid ranges
            .self$setParams()

        }

    )
)

# Adenoma class
###############

risk_of_an_adenoma<-function(risk_params, subject_age){

    r1 <- risk_params$baseline_risk
    aa <- risk_params$age_risk
    r1 <- r1+risk_params$sex_linked_risk
    r2 <- 0
    r3 <- 0

    if (subject_age >= 20){
        if (subject_age>=70){
            k<-4
        } else if (subject_age>=60){
            k<-3
        } else if (subject_age>=50){
            k<-2
        } else {
            k<-1
        }
        A<-c(20,50,60,70,120)
        r2<-subject_age*aa[k]

        r3<-0
        if ( k >= 2 ){
            for (j in 2:k){
                r3<-r3+A[j]*(aa[j-1]-aa[j])
            }
        }
    }

    max(0, min(1, exp(r1+r2+r3)))
}


#' Class implementing Adenoma
#'
#' @export Adenoma
#' @exportClass Adenoma
#' @family CrcSpinModel_classes
Adenoma <- setRefClass( "Adenoma",

    fields = list(
        initiated_in_year                  = "numeric",
        size                               = "numeric",
        location                           = "character",
        state                              = "character",
        beta1                              = "numeric",
        beta2                              = "numeric",
        gamma1                             = "numeric",
        gamma2                             = "numeric",
        gamma3                             = "numeric",
        d10                                = "numeric",
        lambda                             = "numeric",
        sojourn_time                       = "numeric",
        transition_to_preclinical_crc_year = "numeric",
        nu_colon                           = "numeric",
        xi_colon                           = "numeric"
        ),

    methods = list(

        initialize = function (adenoma_params=NULL,
                        initiated_in_year=0,
                        size=1,
                        state="adenoma",
                        ...) {

            # variables prefixed with "l_" are function local versions of
            # class field variables, this avoids warnings being thrown by
            # the R source() function
            #
            # NOTE: this is probably not necessary now that the
            # "AdenomaParams" is being used

            if(is.null(adenoma_params)){
                # no adenoma_model_parameters were passed in
                # create one with default values
                adenoma_params<-AdenomaParams$new()
            }

            a1<-sample(c("cecum","ascending","transverse","descending","sigmoid","rectum"),
                       1, prob =c(0.08,0.23,0.24,0.12,0.24,0.09))
            l_beta1<-runif(1,adenoma_params$beta1_min,adenoma_params$beta1_max)
            l_beta2<-runif(1,adenoma_params$beta2_min,adenoma_params$beta2_max)
            dinfinity<-50
            d0<-1
            l_d10<-l_beta1*((-log(runif(1,0,1)))^(-1/l_beta2))

            if (a1=="rectum"){
                l_gamma1<-runif(1,adenoma_params$gamma1_min,adenoma_params$gamma1_max)
                l_gamma2<-runif(1,adenoma_params$gamma2_min,adenoma_params$gamma2_max)
                l_gamma3<-adenoma_params$gamma3_val
            } else {
                l_gamma1<-runif(1,adenoma_params$gamma1_min,adenoma_params$gamma1_max)
                l_gamma2<-runif(1,adenoma_params$gamma2_min,adenoma_params$gamma2_max)
                l_gamma3<-adenoma_params$gamma3_val
            }

            mean_colon <- runif(1,0.5,5)
            tau_colon  <- runif(1,0.1,1.5)
            l_nu_colon <- sqrt(log(tau_colon^2+1))
            l_xi_colon <- log(mean_colon)-0.5*l_nu_colon^2

            l_lambda<- (-1/l_d10)*(log( (dinfinity-10)/(dinfinity-d0)))

            initFields( initiated_in_year    = initiated_in_year,
                        size                 = size,
                        location             = a1,
                        state                = state,
                        beta1                = l_beta1,
                        beta2                = l_beta2,
                        gamma1               = l_gamma1,
                        gamma2               = l_gamma2,
                        gamma3               = l_gamma3,
                        d10                  = l_d10,
                        lambda               = l_lambda,
                        nu_colon             = l_nu_colon,
                        xi_colon             = l_xi_colon)

        },

        grow = function(subject_age=0){
            #time is since the initiation of the adenoma i.e
            # person@age-object@initiated.in.year originally
            adenoma_age <- subject_age-initiated_in_year
            dinfinity   <- 50
            d0          <- 1
            size       <<- dinfinity-(dinfinity-d0)*exp(-lambda*adenoma_age)
        },

        transition = function(subject_age=0){
            # The function looks after transition of ADENOMA:
            #   adenoma to large adenoma
            #   large adenoma to "CRC"
            #
            # The transition of COLON to "symptomatic CRC"
            # is done in the function Colon::modelTransitionToClinicalCrc()

              if (size >=10 & state=="adenoma"){
                state<<-"large adenoma"
              }

              #if it is already a "CRC" then there is nothing to do
              #if it is an adenoma (or large adenoma)  we check to see if it transitions to a "pre clinical CRC"
              if ( is.element(state , c("adenoma", "large adenoma"))){
                p1<-pnorm( (log(gamma1*size)+gamma2*(initiated_in_year-50))/gamma3)
                dice<-sample(c("transition","no transition"),1,prob=c(p1,1-p1))
                if (dice =="transition"){
                  state<<-"CRC"
                  transition_to_preclinical_crc_year<<-subject_age
                  sojourn_time<<-exp(rnorm(1, mean=xi_colon, sd=nu_colon))
                }
              } #end of transition to "CRC"

        },

        checkCRCSymptomPresentation = function(subject_age=0){
            # This function checks if an adenoma starts causing CRC symptoms to
            # present in a person
            if(state=="CRC"){
                if(transition_to_preclinical_crc_year+sojourn_time <= subject_age){
                        return(TRUE)
                }
            }
            return(FALSE)
        },

        summary = function(){},

        show = function(){

        "Display information about the object.
        "

            cat("Spin Adenoma object of class",
                classLabel(class(.self)), "\n")
        }

        )
)


# Colon Class
#############

# Reference Class
Colon <- setRefClass( "Colon",

        fields = list(
            state="character",
            cancer_site="character",
            risk="CrcRisk",
            sites="list"
            ),

    methods = list(

        initialize = function (host_risk_level="standard",
                            host_sex=sample(c("F","M"),1,prob=c(0.5,0.5)),
                            crcrisk_params=NULL,
                            ...) {

            initFields(
                state="clear",
                cancer_site="",
                risk=CrcRisk$new(subject_sex=host_sex,
                        subject_risk_level=host_risk_level,
                        crcrisk_params=crcrisk_params),
                sites=list())

        },

        # Overwriting this *Type() functions in subclasses
        # allows the Classic CRC Adenoma$initiateAdenomas()
        # function to act like a templated function in
        # other languages. It will cause the initializer of
        # the specified classes to be called when a new adenoma
        # is created. The specific class MUST be a subclass of
        # the Adenoma reference class providing a superset of
        # its initialize function's interface
        getAdenomaType = function () {
            return(getRefClass("Adenoma"))
        },

        initiateAdenomas = function(adenoma_params,
                                subject_age){

            ## temp<-rpois(1, risk_of_an_adenoma(
            ##                     subject_age=subject_age,
            ##                     subject_sex=subject_sex,
            ##                     risk_params=risk))
          temp<-risk_of_an_adenoma(
                                   risk_params=risk,
                                   subject_age=subject_age)

          temp<-sample(c(0,1),1,prob=c(1-temp,temp))

          if (temp>0){
                for ( i in 1:temp){
                    sites <<- c(sites,getAdenomaType()$new(
                                            adenoma_model_params=adenoma_params,
                                            initiated_in_year=subject_age,
                                            size=1,
                                            state="adenoma"
                                            )
                                )
                }
            }
        },

        growAdenomas = function (subject_age) {
            for (i in sites) {
                i$grow(subject_age=subject_age)
            }
        },

        modelAdenomaGrowth = function (adenoma_model_params,
                                subject_age) {

            initiateAdenomas(
                    adenoma_params=adenoma_model_params,
                    subject_age=subject_age)

            growAdenomas(
                    subject_age=subject_age)

        },

        modelAdenomaTransitions = function (subject_age) {
            for( i in sites ){
                i$transition(subject_age=subject_age)
            }
        },

        modelTransitionToClinicalCRC = function (subject_age) {
            if (state=="CRC"){
                if(length(sites)>0){
                    for ( i in sites){
                        if(i$checkCRCSymptomPresentation(subject_age)){
                            state<<-"symptomatic CRC"
                            return(TRUE)
                        }
                    }
                }
            }
            return(FALSE)
        },

        getMostAdvancedAdenomaState = function() {

            if (length(sites)>0){

                tt<-lapply(sites,f<-function(x){x$state})
                tt<-unlist(tt)
                if ( sum(tt=="CRC")>0 ) {
                  adv_state<-"CRC"
                } else if ( sum(tt=="large adenoma")>0 ) {
                  adv_state<-"large adenoma"
                } else {
                  adv_state<-"adenoma"
                }

            } else {
              adv_state<-"clear"
            }

            return(adv_state)

        },

        getCancerSite = function () {
            if (length(sites)>0){
              tt<-lapply(sites,f<-function(x){x$location})
              tt<-unlist(tt)
              cancer_location<-tt[1]
            } else {
                cancer_location<-""
            }

            return(cancer_location)
        },

        isCancerous = function () {
            return(state=="CRC" | state=="symptomatic CRC")
        },

        updateState = function () {

            state<<-getMostAdvancedAdenomaState()

            if (isCancerous()){
                cancer_site<<-getCancerSite()
            }

        },

        show = function(){

        "Display information about the object.
        "

            cat("Spin Colon object of class",
                classLabel(class(.self)), "\n")
        }
    )
)

# Person with colon class
#########################

# Reference Class
PersonWithColon <- setRefClass( "PersonWithColon",
            contains="Person",

    fields = list(
        colon_clinical="character", #"clear", "CRC"
        colon="Colon",
        risk_level="character"
        ),

    methods = list(

        initialize = function(colon=NULL,
                        crcrisk_params=NULL,
                        risk_level="standard",
                        ...){

            ##############################################################
            # call parent class Person$new() to initialise generic Person
            # values
            callSuper(...)

            ##############################################################
            # Do initialisation for PersonWithColon specific variables and
            # member objects
            .self$risk_level <<- risk_level

            # if a colon object was passed in, check if it that IS or was
            # derived from class Colon, then set .self$colon=colon
            if(!is.null(colon)
               && (class(colon)=="Colon"
                   || is.element("Colon", colon$getClass()@refSuperClasses))){
                .self$colon <- colon

            # Otherwise set .self$colon to new Colon object
            } else {
                .self$colon <- Colon$new(
                            host_risk_level=.self$risk_level,
                            host_sex=.self$sex,
                            crcrisk_params=crcrisk_params)
            }
            .self$colon_clinical <<- "clear"
        },

        isDead = function () {
            return(state=="deceased")
        },

        showingCRCSymptoms = function () {
            return(colon_clinical=="symptomatic CRC")
        },

        updateState = function () {
            # update our colon's state
            colon$updateState()
        },

        modelPreClinicalToClinicalTransition = function () {

            if(colon_clinical!="symptomatic CRC") {

                if (colon$modelTransitionToClinicalCRC(age)){
                    colon_clinical<<-"symptomatic CRC"
                }

            }
        },

        testForAndTreatSize = function () {
            return(1)
        },
        
        treatmentRecordSize = function () {
            return(1)
        },
        
        initiateCRCTreatment = function () {
            in_treatment_program<<-"yes"
            return(TRUE)
        },

        medicalSnapshotSize = function () {
            return(17)
        },

        getMedicalSnapshot = function () {

            aa<-rep(FALSE,medicalSnapshotSize())
            if (length(colon$sites)>0){
              tt<-lapply(colon$sites,f<-function(x){x$state})
              tt<-unlist(tt)
              tt<-factor(tt, levels=c("adenoma","large adenoma","CRC"))
              aa[1:3]<-table(tt)
            }

            aa[4]=(colon_clinical=="symptomatic CRC")

            tt<-factor(colon$state,
                  levels=c("clear"  , "adenoma"  , "large adenoma" , "CRC",
                      "symptomatic CRC"))
            aa[5:9]<-table(tt)

            tt<-factor(colon$cancer_site
                    , levels=c("cecum","ascending","transverse","descending","sigmoid","rectum"))
            aa[10:15]<-table(tt)

            aa[16]<-(state=="deceased")
            aa[17]<-(colon_clinical=="symptomatic CRC")
            return(aa)
        },

        show = function(){

        "Display information about the object.
        "

            callSuper()
            cat("Risk level:")
            methods::show(risk_level)
            cat("Colon clinical characteristic:")
            methods::show(colon_clinical)
            cat("Colon:\n")
            colon$show()
        }
    )
)

# vim: set filetype=r fileformat=unix : #
#source("crcSpin.R")

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

        NBCSP = function (person) {
            temp1<-rep(FALSE,person$NBCSPRecordSize())
            if ( (person$age %in% c(50,55,60,65,70)) & ( person$colon_clinical=="clear") #
                &(person$in_treatment_program=="no")){
                                        #the current screening scheme offers iFOBT to people at the ages 50,55,60,65,70.
                                        #We do not offer it if the person already has a diagnosis of "CRC"
                iFOBTscreening(person) #same as .self$iFOBT.screening(person)

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
        
        iFOBTscreening = function(person) {#person is offered iFOBT
            age<-person$age
            test.result<-"none"
            test.state<-"none"
             compliance<-sample(c("accept","decline"),1, prob =c(0.4,0.6))
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

#            treatment_record.2<-NBCSP(person)
            treatment_record.2<-gemini.screening(person)

            return(c(treatment_record.1,treatment_record.2))
        },

        gemini.screening = function(person){
            temp1<-rep(0,14)
            if ( (person$age %in% c(60)) & ( person$colon_clinical=="clear") &(person$in_treatment_program=="no")){
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
        } # blood.test.screening 

        
      ## blood.test.screening = function(person) {#person is offered blood test
      ##       age<-person$age
      ##       test.result<-"none"
      ##       test.state<-"none"
      ##       compliance<-sample(c("accept","decline"),1, prob =c(0.75,0.25))
      ##       if (compliance=="accept"){
      ##           person$updateState()  #object<-get.patient.state(object)
      ##           state<-person$colon$state    #object@colon@state
      ##           if( state=="symptomatic CRC" ){
      ##               sensitivity<-0.86
      ##               test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
      ##               if(test.result=="positive"){
      ##                   test.state<-"TP"
      ##               }
      ##               else{
      ##                   test.state<-"FN"
      ##               }
      ##           } else if ( state== "CRC" ){
      ##               sensitivity<-0.86
      ##               test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
      ##               if(test.result=="positive"){
      ##                   test.state<-"TP"
      ##               }
      ##               else{
      ##                   test.state<-"FN"
      ##               }
      ##           } else if ( state=="large adenoma" ){
      ##               sensitivity<-0.474
      ##               test.result<-sample(c("positive","negative"),1,prob=c(sensitivity,1-sensitivity))
      ##               if(test.result=="positive"){
      ##                   test.state<-"TP"
      ##               }
      ##               else{
      ##                   test.state<-"FN"
      ##               }
      ##           }else if ( state=="adenoma" ){
      ##               specificity<-0.9585
      ##               test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
      ##               if(test.result=="positive"){
      ##                   test.state<-"FP"
      ##               }
      ##               else{
      ##                   test.state<-"TN"
      ##               }
      ##           } else if ( state=="clear" ){
      ##               specificity<-0.968
      ##               test.result<-sample(c("positive","negative"),1,prob=c(1- specificity, specificity))
      ##               if(test.result=="positive"){
      ##                   test.state<-"FP"
      ##               }
      ##               else{
      ##                   test.state<-"TN"
      ##               }
      ##           }#end state =clear
      ##       }
      ##       person$clinical_history$events<-lappend(person$clinical_history$events,
      ##                                               Test$new(
      ##                                                   age=age,
      ##                                                   type="blood",
      ##                                                   compliance=compliance,
      ##                                                   result=test.result,
      ##                                                   state= test.state)
      ##                                               )
      ##   }  # blood.test.screening 
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
