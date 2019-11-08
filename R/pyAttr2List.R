#' Extract attribute information of Python object to R list
#'
#' Only the value like attributes will be extracted rather than private
#' functions, constants or methods (i.e. name starts with '_' will be skipped).
#' Attributes that can be converted to R objects by "reticulate" will just be
#' done like this. Those that are still a Python objects will be replaced by
#' their `str(a)` representation (whether the information makes sense depends on
#' how the class was defined by the original developer).
#'
#' @param py_object A Python object. Will be checked by whether
#' "python.builtin.object" is in `class(py_object)`. Otherwise, the object will
#' be returned untouched
#' @param silent Boolean, default FALSE. There can be some warning message from
#' Python side.
#' @return attrList list
#' @export
#' @examples
#' # In python 'testModule.py'
#' class testClass:
#'     def __init__(self, in1):
#'         self.attr1 = in1
#'         self.attr2 = [in1]
#'     def testFunc(self, in2):
#'         self.attr2.append(in2)
#' # In R
#' > library(reticulate)
#' > testModule <- import('testModule')
#' > obj <- testModule$testClass('hello')
#' > obj$testFunc('world')
#' > pyAttr2List(obj)
#' ## $attr1
#' ## [1] "hello"
#' ##
#' ## $attr2
#' ## [1] "hello" "world"
pyAttr2List <- function(py_object, silent = FALSE) {
    if (!"python.builtin.object" %in% class(py_object)) {
        if (!silent) {
            write('Input object is not identified as Python object. Skipped')
        }
        return(py_object)
    } else {
        pyBuiltins <- reticulate::import_builtins()
        allAttr <- pyBuiltins$dir(py_object)
        output <- list()
        for (i in 1:length(allAttr)) {
            if (!substring(allAttr[i], 1, 1) == "_") {
                # Condition for not an private method
                value <- py_object[[allAttr[i]]]
                if (!"python.builtin.object" %in% class(value)) {
                    output[[allAttr[i]]] <- value
                } else {
                    if (!"python.builtin.method" %in% class(value)) {
                        dataCoverted <- try(value$to_list(), silent = TRUE)
                        if (!class(dataCoverted) == 'try-error') {
                            output[[allAttr[i]]] <- dataCoverted
                        } else {
                            strRepr <- pyBuiltins$str(value)
                            if (!strRepr == "None") {
                            output[[allAttr[i]]] <- strRepr
                            }
                        }
                    }
                }
            }
        }
        return(output)
    }
}
