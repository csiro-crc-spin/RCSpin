#' A collection of extensible classes implementing CRC-Spin and other Spin-alike models
#'
#' This package provides a set of R Reference Classes implementing a family of
#' CRC-Spin models. The classes are highly extensible, and can be used to build
#' other Spin-alike models.
#'
#' The package currently includes the following model classes alone with their
#' associated helper classes:
#' \describe{
#'   \item{GenericModel}{A simple natural aging model. Serves as a base class
#'                       for building more complex models that implement
#'                       stepwise clinical modeling of a study group of
#'                       subjects.}
#'   \item{CrcSpinModel}{An implementation of X et al. classic CRC-Spin model.
#'                       Derives from GenericModel.}
#'   \item{DukesCrcSpinModel}{An implementation of Y et al. extensions of the
#'                            classic CRC-Spin model to allow for transition of
#'                            adenomas through Duke's stages.}
#' }
#'
#' Each class utilises a set of helper classes that describe track and objects
#' such as a Person, an organ (e.g. Colon), a disease (e.g. Cancerous
#' Adenomas), risk parameters, and clinical events.
#'
#' Results of the model runs are available in a study_result matrix that
#' summarises the overall state of the model at each step.
#'
#' @section Usage:
#' Users of the models who do not wish to extend their functionality only need
#' to call a small number of functions to set up the parameters of a model,
#' then run it.
#'
#' To create a new model object, call the model's Reference Class constructor, passing the appropriate arguments. e.g.:
#' \code{m<-CrcSpinModel$new(iterations=99, num_subjects=1000)}
#'
#' If you wish to modify the modeling parameters away from the default, you can
#' call any available \code{set*Parameter} member functions of the new model
#' object. e.g.:
#' \code{...}
#'
#' To run the configured model, call it's \code{run} function:
#' \code{m$run()}
#'
#' Finally, you can extract results from the \code{\link{study_results}}
#' matrix:
#' \code{m$study_results}
#'
#' The meaning of the result matrices rows and columns for each model are
#' described in the respective model class documentation.
#'
#' For more information on the models implemented by this package, please see the following vignettes:
#' \itemize{
#'   \item \link[=]{Classic CRC Spin Model}
#'   \item \link[=]{Duke's CRC Spin Model}
#' }
#'
#' @section Extending the Model Classes:
#' Developers wishing to extend the models may wish to look at the detailed
#' documentation of the main model classes and all their helper class.
#' \emph{Please note that users of the models generally will not need to know
#' this level of detail}.
#'
#' For documentation of included model classes and their helper classes see:
#' \itemize{
#'    \item \code{\link{GenericModel}}
#'    \item \code{\link{CrcSpinModel}}
#'    \item \code{\link{DukesCrcSpinModel}}.
#' }
#'
#' @docType package
#' @name RCSpin
NULL

