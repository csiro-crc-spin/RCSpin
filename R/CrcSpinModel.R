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
        adenoma_model_params = "AdenomaParams",
        screening_flag = "character"
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
                        screening_flag="none",
                        ...) {

            crcrisk_model_params<<-crcRiskParamsType()$new()
            adenoma_model_params<<-adenomaParamsType()$new()
            screening_flag<<- screening_flag

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

            .self$sex_linked_mean <<-   -0.3
            .self$sex_linked_sd   <<-   0.04
            .self$baseline_mean   <<-  -6.7
            .self$baseline_sd     <<-   0.27
            .self$age_mean        <<-   0.03
            .self$age_sd          <<-   0.003

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
#            age_risk <<- rep(0,4)
#            age_risk[1] <<- max(rnorm(1,0.037,0.003 ),0)
#            age_risk[2] <<- max(rnorm(1,0.031,0.003 ),0)
#            age_risk[3] <<- max(rnorm(1,0.029,0.003 ),0)
#            age_risk[4] <<- max(rnorm(1,0.03,0.003 ),0)
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
        beta1_min_colon   = "numeric",
        beta1_min_rectum   = "numeric",
        beta1_max_colon   = "numeric",
        beta1_max_rectum   = "numeric",
        beta2_min_colon   = "numeric",
        beta2_min_rectum   = "numeric",
        beta2_max_colon   = "numeric",
        beta2_max_rectum   = "numeric",
        gamma1_male_colon_min    = "numeric",
        gamma1_male_colon_max    = "numeric",
        gamma2_male_colon_min    = "numeric",
        gamma2_male_colon_max    = "numeric",
        gamma1_male_rectum_min   = "numeric",
        gamma1_male_rectum_max   = "numeric",
        gamma2_male_rectum_min   = "numeric",
        gamma2_male_rectum_max   = "numeric",
        gamma1_female_colon_min  = "numeric",
        gamma1_female_colon_max  = "numeric",
        gamma2_female_colon_min  = "numeric",
        gamma2_female_colon_max  = "numeric",
        gamma1_female_rectum_min = "numeric",
        gamma1_female_rectum_max = "numeric",
        gamma2_female_rectum_min = "numeric",
        gamma2_female_rectum_max = "numeric",
        gamma3_val = "numeric"
    ),

    methods = list(

        setParams = function(...){
                           # Set all the parameters passed in
            initFields(...)

            # Cap certain variables to valid ranges
            beta1_max_colon  <<- max(beta1_min_colon, beta1_max_colon)
            beta2_max_colon  <<- max(beta2_min_colon, beta2_max_colon)
            beta1_max_rectum  <<- max(beta1_min_rectum, beta1_max_rectum)
            beta2_max_rectum  <<- max(beta2_min_rectum, beta2_max_rectum)

            gamma1_male_colon_max     <<- max(gamma1_male_colon_min, gamma1_male_colon_max)
            gamma2_male_colon_max     <<- max(gamma2_male_colon_min, gamma2_male_colon_max)
            gamma1_male_rectum_max    <<- max(gamma1_male_rectum_min, gamma1_male_rectum_max)
            gamma2_male_rectum_max    <<- max(gamma2_male_rectum_min, gamma2_male_rectum_max)
            gamma1_female_colon_max   <<- max(gamma1_female_colon_min, gamma1_female_colon_max)  
            gamma2_female_colon_max   <<- max(gamma2_female_colon_min, gamma2_female_colon_max)  
            gamma1_female_rectum_max  <<- max(gamma1_female_rectum_min, gamma1_female_rectum_max)
            gamma2_female_rectum_max  <<- max(gamma2_female_rectum_min, gamma2_female_rectum_max)
        },

        initialize = function(){

            # Set default values here
        .self$beta1_min_colon   <<- 24.3
        .self$beta1_min_rectum  <<- 5.9
        .self$beta1_max_colon   <<- 34.2
        .self$beta1_max_rectum  <<- 13.8
        .self$beta2_min_colon   <<- 1.1
        .self$beta2_min_rectum  <<- 1.4
        .self$beta2_max_colon   <<- 3.9
        .self$beta2_max_rectum  <<- 3.9
        .self$gamma1_male_colon_min     <<- 0.04
        .self$gamma1_male_colon_max     <<-  0.049
        .self$gamma2_male_colon_min     <<- 0.002
        .self$gamma2_male_colon_max     <<- 0.016
        .self$gamma1_male_rectum_min    <<- 0.021
        .self$gamma1_male_rectum_max    <<- 0.049
        .self$gamma2_male_rectum_min    <<- 0.001
        .self$gamma2_male_rectum_max    <<- 0.019
        .self$gamma1_female_colon_min   <<- 0.044
        .self$gamma1_female_colon_max   <<- 0.05
        .self$gamma2_female_colon_min   <<- 0.000
        .self$gamma2_female_colon_max   <<- 0.013
        .self$gamma1_female_rectum_min   <<- 0.03
         .self$gamma1_female_rectum_max   <<- 0.05
        .self$gamma2_female_rectum_min   <<- 0.008
        .self$gamma2_female_rectum_max   <<- 0.019
        .self$gamma3_val                 <<- 0.5

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
    if (subject_age < 20){r1<-0}
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
        xi_colon                           = "numeric",
        div                                = "numeric",
        p1.i.minus.1                       = "numeric",
        temp                              = "character"
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
            dinfinity<-50
            d0<-1

            if (a1=="rectum"){
                mean_colon <- runif(1,1.1,4.7)
                tau_colon  <- runif(1,0.15,1.4)
                l_beta1<-runif(1,adenoma_params$beta1_min_rectum,adenoma_params$beta1_max_rectum)
                l_beta2<-runif(1,adenoma_params$beta2_min_rectum,adenoma_params$beta2_max_rectum)
                l_d10<-l_beta1*((-log(runif(1,0,1)))^(-1/l_beta2))
#                if (subject$sex=="M"){
                    l_gamma1<-runif(1,adenoma_params$gamma1_male_rectum_min,adenoma_params$gamma1_male_rectum_max)
                    l_gamma2<-runif(1,adenoma_params$gamma2_male_rectum_min,adenoma_params$gamma2_male_rectum_max)
                    l_gamma3<-adenoma_params$gamma3_val
#                } else {
#                    l_gamma1<-runif(1,adenoma_params$gamma1_female_rectum_min,adenoma_params$gamma1_female_rectum_max)
#                    l_gamma2<-runif(1,adenoma_params$gamma2_female_rectum_min,adenoma_params$gamma2_female_rectum_max)
#                    l_gamma3<-adenoma_params$gamma3_val
#                }
            } else {
                mean_colon <- runif(1,1.0,3.9)
                tau_colon  <- runif(1,0.15,1.4)
                l_beta1<-runif(1,adenoma_params$beta1_min_colon,adenoma_params$beta1_max_colon)
                l_beta2<-runif(1,adenoma_params$beta2_min_colon,adenoma_params$beta2_max_colon)
                l_d10<-l_beta1*((-log(runif(1,0,1)))^(-1/l_beta2))
 #               if (subject$sex=="M"){
                    l_gamma1<-runif(1,adenoma_params$gamma1_male_colon_min,adenoma_params$gamma1_male_colon_max)
                    l_gamma2<-runif(1,adenoma_params$gamma2_male_colon_min,adenoma_params$gamma2_male_colon_max)
                    l_gamma3<-adenoma_params$gamma3_val
 #               } else {
 #                   l_gamma1<-runif(1,adenoma_params$gamma1_female_colon_min,adenoma_params$gamma1_female_colon_max)
 #                   l_gamma2<-runif(1,adenoma_params$gamma2_female_colon_min,adenoma_params$gamma2_female_colon_max)
 #                   l_gamma3<-adenoma_params$gamma3_val
 #               }
            }

            
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
                        xi_colon             = l_xi_colon,
                        div                  = 1,
                        p1.i.minus.1         = 0,
                        temp                 = tempfile(fileext = ".error") 
                       )

        }, #end of initialize

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
                q1<-p1 - p1.i.minus.1
                div <<-div*(1-q1)
                q1 <- q1/div
                p1.i.minus.1 <<- p1
                dice<-sample(c("transition","no transition"),1,prob=c(q1,1-q1))
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
