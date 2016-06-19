# vim: set filetype=r fileformat=unix : #

###########################
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
Person <-  setRefClass( "Person",

    fields = list(
        age = "numeric",
        sex = "character",
        state = "character", #"deceased","living", etc (as needed by subclasses)
        in_treatment_program = "character",
        clinical_history = "ClinicalHistory",
        study_id = "numeric",
        random_seed_state = "integer",
        NBCSP.propensity = "numeric",
        BSA.propensity = "numeric"
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
                            NBCSP.propensity=NA,
                            BSA.propensity=NA,
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

            if (is.na(NBCSP.propensity)){
                NBCSP.propensity_l  <- rbeta(1, shape1=0.198,shape2=0.3, ncp = 0) 
            } else{
                NBCSP.propensity_l <- NBCSP.propensity
            }
            
            if (is.na(BSA.propensity)){
                BSA.propensity_l  <- rbeta(1, shape1=0.02, shape2=0.3, ncp = 0) 
            } else {
                BSA.propensity_l  <-  BSA.propensity
            }
                
            
            initFields(age = age,
                sex = s,
                state = state,
                in_treatment_program = in_treatment_program,
                clinical_history = clinical_history,
                study_id = study_id,
                random_seed_state = random_seed,
                NBCSP.propensity = NBCSP.propensity_l,
                BSA.propensity = BSA.propensity_l)
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
 #           methods::show(syste m_data)
            cat("Clinical history:\n")
            clinical_history$show()
            cat("BSA propensity to screen:")
            methods::show(BSA.propensity)
        }
        
        ## set.BSA.propensity =function(){
        ##     "set the propensity to screen with BSA
        ## "
        ##     u<- rbeta(1, shape1=0.02, shape2=0.3, ncp = 0) 
        ## }
        

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
