setClassUnion("numericOrNULL", c("numeric", "NULL"))

setClass("CellNumSol", slots=
             c(solution="numericOrNULL",
               evaluation_points = "numeric",
               AUC="numeric",
               regr_coefficients="numeric",
               r_squared="numeric",
               n_cells = "integer",
               last_iteration="integer",
               args="list")
)

