# Functions will live here

#' Title
#'
#' @param key
#'
#' @return
#' @export
#'
#' @examples
#' full_dataset(key = "SVD_REDCAP_API")
full_dataset <- function(key) {
  REDCapCAST::read_redcap_tables(
    uri = "https://redcap.au.dk/api/",
    token = keyring::key_get(key),
    fields = c(
      "record_id",
      "trial_name",
      "svd_perf",
      "svd_missing_reason"
    )
  )
}

#' Title
#'
#' @param data
#' @param export.path
#'
#' @return
#' @export
#'
#' @examples
#' data <- targets::tar_read("ls_data")
#' data |> create_flowchart()
create_flowchart <- function(data, export.path = NULL) {
  ds <- data[["svd_score"]] |>
    # Just filtering for repeat instance #1 would be much easier
    # This is now prepared for the future with other filtering
    (\(.x){
      split(.x, .x$record_id)
    })() |>
    purrr::map(function(.x) {
      .x[1, ]
    }) |>
    dplyr::bind_rows() |>
    dplyr::select(record_id, svd_perf, svd_missing_reason) |>
    dplyr::full_join(data[["basis"]]) |>
    dplyr::mutate(
      assessed = dplyr::if_else(!is.na(svd_perf), "Assessed", "Not yet assessed"),
      include = dplyr::if_else(is.na(svd_missing_reason), assessed, svd_missing_reason),
      include = dplyr::if_else(include == "Assessed", NA, include)
    )

  out <- ds |> consort::consort_plot(
    orders = c(
      record_id = "All subjects",
      include = "Excluded",
      record_id = "Final dataset"
    ),
    side_box = c("include"),
    labels = c(
      "1" = "Pooled subjects with AIS",
      "2" = "Subjects available for analysis"
    )
  )

  if (!is.null(export.path)) {
    out |> export_consort_dot(path = export.path)
  } else {
    plot(out)
  }
}


export_grViz_png <- function(data, path) {
  plot(data, grViz = TRUE) |>
    DiagrammeRsvg::export_svg() |>
    charToRaw() |>
    rsvg::rsvg_png(file = path)
}

#' Title
#'
#' @param data
#' @param path
#'
#' @return
#' @export
#'
#' @examples
export_consort_dot <- function(data, path) {
  ls <- data |>
    consort:::plot.consort(grViz = TRUE)
  writeLines(ls$x$diagram, con = here::here(path))
}
