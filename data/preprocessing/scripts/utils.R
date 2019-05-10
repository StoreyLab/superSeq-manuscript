col_name <- function (x, default = stop("Please supply column name", call. = FALSE)) 
{
    if (is.character(x)) 
        return(x)
    if (identical(x, quote(expr = ))) 
        return(default)
    if (is.name(x)) 
        return(as.character(x))
    if (is.null(x)) 
        return(x)
    stop("Invalid column specification", call. = FALSE)
}

sample_unique_ <- function(.data, column, n = 1) {
    col <- .data[[column]]
    u <- sample(unique(col), n)
    .data[col %in% u, ]
}

sample_unique <- function(.data, column, ...) {
    sample_unique_(.data, col_name(substitute(column)))
}
