#' pdataSummary
#'
#'
#' @title pdataSummary
#' @param data pdata
#' @param group group
#' @param show_p show p
#' @importFrom gtsummary tbl_summary
#' @importFrom gtsummary add_n
#' @importFrom gtsummary modify_header
#' @importFrom gtsummary bold_labels
#' @importFrom gtsummary add_p
#' @importFrom magrittr %>%
#' @return a summary data frame
#' @export
#' @author Yuanlong Hu

pdataSummary <- function(data, group, show_p=T){

  data$group <- group
  res <- tbl_summary(
    data,
    by = "group",
    missing = "no"
  ) %>%
    add_n() %>%
    modify_header(label = "**Variable**") %>%
    bold_labels()

  if(show_p) res <- res %>% add_p()
  return(res)
}
