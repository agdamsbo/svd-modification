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
      "trial_name"
    ), forms = c("svd_score", "consensus_score")
  )
}

#' Merge consensus assessments with raw scores and clean
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' data <- targets::tar_read("ls_data")
#' data |>
#'   final_score() |>
#'   (\(.x) summary(factor(.x$consensus)))()
final_score <- function(data) {
  main_score <- c(
    "microbleed",
    "siderose",
    "lacunes",
    "wmh",
    "atrophy"
  )
  sub_score <- c(
    "microbleed_location___1",
    "microbleed_location___2",
    "microbleed_location___3"
  )
  score_items <- c(main_score, sub_score)

  svd <- data |>
    purrr::pluck("svd_score") |>
    (\(.x){
      split(.x, .x$record_id)
    })() |>
    purrr::map(function(.x) {
      .x |>
        dplyr::select(c(
          "record_id",
          "svd_user",
          paste0("svd_", score_items)
        )) |>
        dplyr::slice(c(1:2))
    })

  svd_scored <- svd |>
    (function(.x) {
      # Creating filter to subset only where 1 assessor has filled, as this indicates
      # that the scoring has been performed at all
      filt <- .x |> purrr::map_vec(function(.y) {
        !all(apply(dplyr::select(.y[1, ], paste0("svd_", main_score)), 1, is.na))
      })
      .x[filt]
    })()

  con_scored <- data |>
    purrr::pluck("consensus_score") |>
    dplyr::tibble() |>
    dplyr::select(c(
      "record_id",
      "consensus_user",
      "consensus_microbleed",
      "consensus_microbleed_location___1",
      "consensus_microbleed_location___2",
      "consensus_microbleed_location___3",
      "consensus_siderose",
      "consensus_lacunes",
      "consensus_wmh",
      "consensus_atrophy"
    )) |>
    remove_prefix("consensus_")

  svd_consensus_merge <- svd_scored |>
    purrr::imap(function(.x, .i) {
      ## Splits steps for easier reading
      ## First replaces the few missings (they should not be there)
      ## This step can hopefully be deleted later
      ds_test <- .x |>
        dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(.x, "missing")))

      # Filtering
      complete <- dplyr::bind_rows(
        remove_prefix(ds_test, prefix = "svd_"),
        remove_prefix(dplyr::filter(con_scored, record_id == .i), prefix = "consensus_")
      )
    })


  .x <- svd_consensus_merge |> purrr::pluck("svd_1")

  out <- svd_consensus_merge |>
    purrr::map(function(.x) {
      agree_index <- .x[-nrow(.x), ] |>
        # dplyr::select(score_items) |>
        # Keeping all variables and set user and locations true to ignore in the correction
        # Handling disagreement in location assessment is minor priority
        dplyr::mutate(
          user = TRUE,
          dplyr::across(dplyr::contains("location"), ~TRUE)
        ) |>
        test_agreement()

      agree <- agree_index |>
        all()

      if (nrow(.x) > 2) {
        if (agree) {
          tibble::tibble(.x[1, ]) |>
            dplyr::mutate(consensus = "agreed")
        } else {
          # Consensus user is NA if not completed yet
          if (is.na(.x[3, 2])) {
            # Consensus needed, but not yet completed
            tibble::tibble(.x[1, ]) |>
              dplyr::mutate(consensus = "disagree")
          } else {
            # Consensus merge
            ifelse(agree_index, .x[1, ], .x[3, ]) |>
              dplyr::bind_cols() |>
              dplyr::mutate(consensus = "consensus")
          }
        }
      } else {
        # If only assessed once, that assessment is kept for now
        remove_prefix(.x[1, ], prefix = "svd_") |>
          dplyr::mutate(consensus = "nonassessed")
      }
    }) |>
    dplyr::bind_rows()

  ## Cleaned and consensus substituted scores
  out
}


#' Testing agreement across tibble/df/matrix
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
test_agreement <- function(data, MARGIN = 2) {
  apply(data, MARGIN = MARGIN, function(.x) {
    length(unique(.x)) == 1
  })
}

#' Remove column name prefix if present
#'
#' @param data
#' @param prefix
#'
#' @return
#' @export
#'
#' @examples
remove_prefix <- function(data, prefix) {
  setNames(data, gsub(prefix, "", names(data)))
}


#' Create flowchart from data based on `consort`
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
      svd_missing_reason = dplyr::case_match(
        svd_missing_reason,
        "Manglende billeder" ~ "No acute MR performed",
        "Ingen adgang til billeder" ~ "MR not available",
        "Andet" ~ "Other reason"
      ),
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

#' Export dot file for most flexible plotting (over png)
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
