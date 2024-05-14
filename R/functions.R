# Functions will live here

# Functions from the SVD-burden annotation project
# source("https://raw.githubusercontent.com/agdamsbo/cSVD-burden-annotation-and-analysis-project/main/R/functions.R")

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
      paste0("consensus_", score_items)
    )) |>
    remove_prefix("consensus_")

  svd_consensus_merge <- svd_scored |>
    purrr::imap(function(.x, .i) {
      # Binding rows
      dplyr::bind_rows(
        remove_prefix(.x, prefix = "svd_"),
        remove_prefix(dplyr::filter(con_scored, record_id == .i), prefix = "consensus_")
      )
    })


  # .x <- svd_consensus_merge |> purrr::pluck("svd_92")

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
            # Consider a coded handling of all these disagreements, as the work
            # to make consensus is huge
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
          dplyr::mutate(consensus = "single")
      }
    }) |>
    dplyr::bind_rows()

  bonus <- data |>
    purrr::pluck("svd_score") |>
    (\(.x){
      split(.x, .x$record_id)
    })() |>
    purrr::map(function(.x) {
      .x |>
        dplyr::select(c(
          "record_id",
          "svd_missing_reason"
        )) |>
        dplyr::slice(1)
    }) |>
    dplyr::bind_rows()


  ## Cleaned and consensus substituted scores
  dplyr::full_join(out, bonus)
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
        "Andet" ~ "Other reasons"
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

################################################################################
########
########      Clinical data collection
########
################################################################################

svd_merging_old <- function(token) {
  REDCapR::redcap_read(
    redcap_uri = "https://redcap.au.dk/api/",
    token = token,
    fields = c("record_id", "trial_id", "trial_name")
  )$data
}

#' Imports specified TALOS data
#'
#' @param token api token
#' @param vars variable names to retrieve (fields in REDCap lingo)
#'
#' @return tibble of REDCap exported data
#' @examples
#' data_svd <- svd_merging_old(keyring::key_get("SVD_REDCAP_API"))
#' data <- data_svd |>
#'   dplyr::filter(trial_name == "TALOS") |>
#'   dplyr::pull(trial_id) |>
#'   clinical_talos()
#' skimr::skim(talos_clin)
clinical_talos <- function(records) {
  REDCapR::redcap_read(
    redcap_uri = "https://redcap.au.dk/api/",
    token = keyring::key_get("TALOS_REDCAP_API"),
    records = records,
    fields = c(
      "record_id",
      # mRS
      "talos_mrs00_4",
      "talos_mrs01_4",
      "rdate",
      "cpr",
      "talos_nihss16_0",
      "rtreat"
    ),
    forms = c("reg", "pase_0")
  )$data
}

#' Title
#'
#' @param data
#' @param nas
#'
#' @return
#' @export
#'
#' @examples
#' data |> talos_data_enrich()
talos_data_enrich <- function(data,
                              nas = c("9.", "", "9. NA", "Not available", "Not relevant (filter/fold)")) {
  data_nas <- data |> dplyr::mutate(dplyr::across(
    dplyr::starts_with("talos"),
    ~ dplyr::if_else(.x %in% nas, NA, .x)
  ))

  out <- data_nas |>
    (\(.x){
      .x |>
        dplyr::select(dplyr::starts_with("talos_pase")) |>
        (\(.y){
          .y |> dplyr::select(-grep("(00|12|x)_0$", names(.y)))
        })() |>
        stRoke::pase_calc() |>
        (\(.y){
          tibble::tibble(.x, .y)
        })()
    })() |>
    dplyr::mutate(
      dob = as.Date(stRoke::cpr_dob(cpr), format = "%d-%m-%Y"),
      age = stRoke::age_calc(dob, as.Date(rdate)),
      female_sex = stRoke::cpr_female(cpr)
    ) |>
    dplyr::select(-cpr) |>
    dplyr::select(-dob)

  out
}


#' Formatting and subsetting clinical data from TALOS
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' data |>
#'   talos_data_enrich() |>
#'   skimr::skim()
#' data |>
#'   talos_data_enrich() |>
#'   format_talos() |>
#'   skimr::skim()
format_talos <- function(data) {
  # Replaces all NAs in factors, and rejoins data frame
  # Factors keeps old levels; doesn't matter when exported as .csv and reimported...
  out <- data |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::starts_with("reg_"),
        ~ as.character(.x)
      ),
      smoker = dplyr::case_match(
        reg_rygning, "1" ~ TRUE,
        c("2", "3", "4") ~ FALSE,
        "9" ~ NA
      ),
      # Living alone defined as not together with somebody
      reg_alone = dplyr::case_match(
        reg_civil, "1" ~ FALSE,
        c("2", "3") ~ TRUE,
        "9" ~ NA
      ),
      reg_more_alc = dplyr::case_match(
        reg_alkohol, "1" ~ FALSE,
        "2" ~ TRUE,
        "9" ~ NA
      ),
      dplyr::across(
        .cols = c(
          "reg_hyperten",
          "reg_diabetes",
          "reg_atriefli",
          "reg_perifer_arteriel",
          "reg_tidl_tci",
          "reg_ami"
        ),
        ~ dplyr::case_match(
          .x, "1" ~ TRUE,
          "2" ~ FALSE,
          "9" ~ NA
        )
      ),
      dplyr::across(
        .cols = c(
          "reg_trombolyse",
          "reg_trombektomi"
        ),
        ~ dplyr::case_match(
          .x, "1" ~ TRUE,
          c("3", "4") ~ FALSE,
          "9" ~ NA
        )
      ),
      any_perf = reg_trombolyse | reg_trombektomi
    ) |>
    # Selecting and renaming in one step
    dplyr::transmute(
      trial_id = as.character(record_id),
      active_treatment = rtreat == "Aktiv (Citalopram)",
      trial = "TALOS",
      age,
      female_sex,
      inclusion_date = as.Date(rdate),
      pase_0 = score_sum,
      smoker,
      alc_more = reg_more_alc,
      alone = reg_alone,
      diabetes = reg_diabetes,
      hyperten = reg_hyperten,
      pad = reg_perifer_arteriel,
      afib = reg_atriefli,
      ami = reg_ami,
      ais = FALSE, # All are first time stroke patients
      tci = reg_tidl_tci,
      tpa = reg_trombolyse,
      evt = reg_trombektomi,
      any_perf,
      nihss = talos_nihss16_0,
      mrs_eos = substr(talos_mrs01_4, 1, 1)
    )

  out
}

resist_pase_items <- function() {
  m <- REDCapR::redcap_metadata_read(
    redcap_uri = "https://redcap.au.dk/api/",
    token = keyring::key_get(
      service = "redcapAPI",
      username = "RESIST-MIGRATION",
      keyring = "REDCAP_APIs"
    )
  )$data

  m |>
    dplyr::filter(grepl("^v\\d*?._", field_name)) |>
    dplyr::pull(field_name)
}



#' Imports specified data
#'
#' @param token api token
#' @param vars variable names to retrieve (fields in REDCap lingo)
#'
#' @return tibble of REDCap exported data
#' @examples
#' data_svd <- svd_merging_old(keyring::key_get("SVD_REDCAP_API"))
#' data <- data_svd |>
#'   dplyr::filter(trial_name == "RESIST") |>
#'   dplyr::pull(trial_id) |>
#'   clinical_resist()
#' skimr::skim(data)
clinical_resist <- function(records) {
  REDCapR::redcap_read(
    redcap_uri = "https://redcap.au.dk/api/",
    token = keyring::key_get(
      service = "redcapAPI",
      username = "RESIST-MIGRATION",
      keyring = "REDCAP_APIs"
    ),
    records = records,
    fields = c(
      "resistid",
      "treatment_arm",
      "civil_status",
      "bolig",
      "smoking",
      "alcohol",
      "ht_art_dicho",
      "diabetes_dicho",
      "afib_dicho",
      "ais_dicho",
      "tia_dicho",
      "ami_dicho_prior",
      "pad_dicho",
      "nihss",
      "mrs3_final",
      "pase_total_score",
      "tpa",
      "evt",
      "cpr",
      "rand_dato_tid",
      resist_pase_items()
    )
  )$data
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' data |> resist_data_enrich()
resist_data_enrich <- function(data) {
  data_nas <- data |>
    # Simple sub-setting by index N
    dplyr::mutate(dplyr::across(
      c(
        "civil_status",
        "bolig",
        "smoking",
        "alcohol",
        "ht_art_dicho",
        "diabetes_dicho",
        "afib_dicho",
        "ais_dicho",
        "tia_dicho",
        "mrs3_final",
        "ami_dicho_prior",
        "pad_dicho",
        "tpa",
        "evt"
      ),
      ~ dplyr::if_else(!.x %in% 0:9, NA, .x)
    ))

  out <- data_nas |>
    dplyr::mutate(
      dob = as.Date(stRoke::cpr_dob(cpr), format = "%d-%m-%Y"),
      age = stRoke::age_calc(dob, as.Date(rand_dato_tid)),
      female_sex = stRoke::cpr_female(cpr)
    ) |>
    dplyr::select(-cpr) |>
    dplyr::select(-dob)

  out
}


#' Formatting and subsetting clinical data from TALOS
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' data |>
#'   resist_data_enrich() |>
#'   skimr::skim()
#' data |>
#'   resist_data_enrich() |>
#'   format_resist() |>
#'   skimr::skim()
format_resist <- function(data) {
  # Replaces all NAs in factors, and rejoins data frame
  # Factors keeps old levels; doesn't matter when exported as .csv and reimported...
  out <- data |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::everything(),
        ~ as.character(.x)
      ),
      smoker = dplyr::case_match(
        smoking, "1" ~ TRUE,
        c("2", "3", "4") ~ FALSE,
        "9" ~ NA
      ),
      # Living alone defined as not together with somebody
      alone = dplyr::case_match(
        civil_status, "1" ~ FALSE,
        c("2", "3") ~ TRUE,
        "9" ~ NA
      ),
      alc_more = dplyr::case_match(
        alcohol, "1" ~ FALSE,
        "2" ~ TRUE,
        "9" ~ NA
      ),
      dplyr::across(
        .cols = c(
          "ht_art_dicho",
          "diabetes_dicho",
          "afib_dicho",
          "ais_dicho",
          "tia_dicho",
          "ami_dicho_prior",
          "pad_dicho",
          "tpa",
          "evt"
        ),
        ~ dplyr::case_match(
          .x, "1" ~ TRUE,
          "0" ~ FALSE,
          "9" ~ NA
        )
      ),
      any_perf = tpa | evt
    ) |>
    # Selecting and renaming in one step
    dplyr::transmute(
      trial_id = as.character(resistid),
      active_treatment = treatment_arm == 1,
      trial = "RESIST",
      age,
      female_sex,
      inclusion_date = as.Date(rand_dato_tid),
      pase_0 = pase_total_score,
      smoker,
      alc_more,
      alone,
      hyperten = ht_art_dicho,
      diabetes = diabetes_dicho,
      afib = afib_dicho,
      tci = tia_dicho,
      ais = ais_dicho,
      ami = ami_dicho_prior,
      pad = pad_dicho,
      nihss,
      tpa,
      evt,
      any_perf,
      mrs_eos = mrs3_final
    )

  out
}

clin_data_all <- function() {
  data_svd <- svd_merging_old(keyring::key_get("SVD_REDCAP_API")) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.character(.x)))

  ls <- list(
    talos = data_svd |>
      dplyr::filter(trial_name == "TALOS") |>
      dplyr::pull(trial_id) |>
      clinical_talos() |>
      talos_data_enrich() |>
      format_talos(),
    resist = data_svd |>
      dplyr::filter(trial_name == "RESIST") |>
      dplyr::pull(trial_id) |>
      clinical_resist() |>
      resist_data_enrich() |>
      format_resist()
  )

  ls |>
    purrr::map(function(.y) {
      .y |> dplyr::mutate(
        dplyr::across(dplyr::all_of(c("age", "pase_0", "nihss")), ~ as.numeric(.x)),
        dplyr::across(dplyr::all_of(c("female_sex")), ~ as.logical(.x))
      )
    }) |>
    purrr::map(dplyr::left_join, dplyr::select(data_svd, dplyr::all_of(c("record_id", "trial_id")))) |>
    dplyr::bind_rows() |> 
    dplyr::mutate(
      pase_0_high = pase_0 > median(as.numeric(pase_0), na.rm = TRUE)
    )
}



#' Simplified SVD score 0-4
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' data <- targets::tar_read("ls_data") |>
#'   purrr::pluck("svd_score") |>
#'   (\(.x)split(.x, .x$svd_user))() |>
#'   purrr::pluck(1)
#' data |>
#'   cut_scores() |>
#'   dplyr::select(record_id, dplyr::starts_with("simple_"))
cut_scores <- function(data,
                       named.cuts = setNames(
                         list(1, 1, 2, 2),
                         c("microbleed", "lacunes", "wmh", "atrophy")
                       ),
                       prefix = "svd_") {
  ds_num <- data |>
    tibble::tibble() |>
    dplyr::select(dplyr::all_of(paste0(prefix, names(named.cuts)))) |>
    setNames(names(named.cuts)) |>
    # na.omit() |>
    # Subset with 9 is dirty trick to order factors correctly to translate to nums
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ gsub(">", "9", .x) |>
      as.factor() |>
      as.numeric() - 1))

  ds_num |>
    purrr::imap(function(.x, .i) {
      dplyr::if_else(.x < named.cuts[[.i]], 0, 1)
    }) |>
    dplyr::bind_cols() |>
    dplyr::rowwise() |>
    dplyr::mutate(score = sum(dplyr::c_across(dplyr::everything()))) |>
    dplyr::ungroup() |>
    (\(.x) setNames(.x, paste0("simple_", names(.x))))() |>
    (\(.x){
      dplyr::bind_cols(data, .x)
    })()
}

################################################################################
########
########      Helper functions
########
################################################################################


dsub <- function(data, pattern, replacement) {
  gsub(pattern = pattern, replacement = replacement, x = data, )
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
#' data <- targets::tar_read("df_complete")
#' data |> clinical_flowchart()
clinical_flowchart <- function(data, export.path = NULL) {
  ds <- data |>
    dplyr::mutate(
      svd_missing_reason = dplyr::case_match(
        svd_missing_reason,
        "Manglende billeder" ~ "No acute MR performed",
        "Ingen adgang til billeder" ~ "MR not available",
        "Andet" ~ "Other reasons"
      ),
      assessed = dplyr::if_else(!is.na(simple_score), "Assessed", "Not assessed"),
      include = dplyr::if_else(is.na(svd_missing_reason), assessed, svd_missing_reason),
      include = dplyr::if_else(include == "Assessed", NA, include),
      pase_available = dplyr::if_else(is.na(pase_0), "Missing", NA)
    )

  out <- ds |> consort::consort_plot(
    orders = c(
      record_id = "All subjects",
      include = "No MR",
      record_id = "Patients considered",
      pase_available = "No PASE available",
      record_id = "Patients included"
    ),
    side_box = c("include", "pase_available"),
    labels = c(
      "1" = "Pooled subjects with AIS",
      "2" = "Subjects available for initial analysis",
      "3" = "Subjects available for final analysis"
    )
  )

  if (!is.null(export.path)) {
    out |> export_consort_dot(path = export.path)
  } else {
    plot(out)
  }
}

#' Store of variable names for data sub-setting
#'
#' @param set character vector of var group. Any of c("pre","clin","post")
#'
#' @return list
#' @examples
#' define_variables()
define_variables <- function(set) {
  ls <- list(
    poster = c(
      "simple_score",
      "age",
      "female_sex",
      "pase_0_high",
      "alone",
      "smoker",
      "hyperten",
      "diabetes",
      "ais",
      "afib"
    ),
    pre = c(
      "simple_score",
      "age",
      "female_sex",
      "pase_0",
      "alone",
      "smoker",
      "alc_more",
      "hyperten",
      "diabetes",
      "tci",
      "ais",
      "ami",
      "afib",
      "pad"
    ),
    clin = c(
      "active_treatment",
      "trial",
      "nihss",
      "tpa",
      "evt"
    ),
    post = c("mrs_eos")
  )
}

#' Get var names in vector from group names. Possibility to keep all vars for as log as possible. Can be supplied to `gtsummary` functions
#'
#' @param groups vector of group names. See names(define_variables()) for options
#'
#' @return
#' @export
#'
#' @examples
get_var_vec <- function(v.groups) {
  define_variables()[{{ v.groups }}] |> purrr::list_c()
}

#' Subsets dataset based on variable group names as defined
#'
#' @param data
#' @param vars.groups character vector of var group. Any of c("pre","clin","post")
#'
#' @return character vector
#'
#' @examples
#' targets::tar_read(df_all_data_formatted) |> get_vars(c("universal", "events"))
get_vars <- function(data, vars.groups) {
  data |> dplyr::select(tidyselect::all_of(get_var_vec(vars.groups)))
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' data <- targets::tar_read(tbl_preds_lin_imp_reg)
#' data |> fix_labels()
fix_labels <- function(data) {
  cls <- class(data)
  data[[1]][["variable"]][data[[1]][["row_type"]] == "label"] |>
    subset_named_labels(var_labels()) |>
    unname() -> data[[1]][["label"]][data[[1]][["row_type"]] == "label"]
  class(data) <- cls
  data
}

################################################################################
########
########       Inter rater reliability calculations
########
################################################################################

#' Simplified SVD score 0-4
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
#' data <- read_instrument() |> inter_rater_data()
#' data |> simple_score()
simple_score <- function(data) {
  data |>
    dplyr::transmute(record_id, redcap_repeat_instance,
      microbleed = dplyr::if_else(svd_microbleed < 1, 0, 1),
      lacunes = dplyr::if_else(svd_lacunes < 1, 0, 1),
      wmh = dplyr::if_else(svd_wmh < 2, 0, 1),
      atrophy = dplyr::if_else(svd_atrophy < 2, 0, 1),
      score = microbleed + lacunes + wmh + atrophy
    )
}


#' Inter rater reliability calculations
#'
#' @param data data
#'
#' @return
#' @export tibble
#'
#' @examples
#' irr_sample |> inter_rater_calc()
#' data <- read_instrument() |>
#'   inter_rater_data()
#' data |> inter_rater_calc()
inter_rater_calc <- function(data, coder_var = redcap_repeat_instance) {
  data |> tidycomm::test_icr(
    unit_var = record_id,
    coder_var = {{ coder_var }},
    holsti = FALSE,
    fleiss_kappa = TRUE,
    brennan_prediger = TRUE,
    na.omit = TRUE
  )
}

#' ICC calculations
#'
#' @param data minimal dataset with only relevant variables
#' @param unit_var subject var
#' @param coder_var rater var
#'
#' @return
#' @export
#'
#' @examples
#' # ds_simple |> icc_multi(unit_var=record_id, coder_var=svd_user)
icc_multi <- function(data, unit_var = record_id, coder_var = redcap_repeat_instance) {
  # The function to calculate ICC
  icc_calc <- function(data) {
    irr::icc(data, model = "twoway", type = "agreement", unit = "single") |>
      purrr::pluck("value")
  }

  # Names of provided variables
  suppressWarnings(nms <- data |>
    dplyr::select(-{{ unit_var }}, -{{ coder_var }}) |>
    names())

  # ICC calculation for each variable
  nms |>
    lapply(function(.x) {
      tidycomm:::unit_coder_matrix(data,
        unit_var = {{ unit_var }},
        coder_var = {{ coder_var }},
        test_var = .x
      )
    }) |>
    purrr::map(icc_calc) |>
    purrr::list_c() |>
    (\(.y){
      tibble::tibble(
        Variable = nms,
        IntraclCorrCoef = .y
      )
    })()
}


#' Join IRR and ICC calculations
#'
#' @param data data
#'
#' @return tibble
#' @export
irr_icc_calc <- function(data) {
  dplyr::left_join(inter_rater_calc(data), icc_multi(data)) |>
    dplyr::select(-tidyselect::all_of(c("n_Coders", "n_Categories", "Level", "n_Units")))
}


#' Functionalised character vector of all labels
#'
#' @return
#' @export
#'
#' @examples
var_labels <- function() {
  c(
    age = "Age",
    female_sex = "Female sex",
    simple_score = "SVD score",
    bmi = "Body mass index",
    smoker = "Current smoker",
    alone = "Living alone",
    alc_more = "High alcohol consumption",
    hyperten = "Hypertension",
    diabetes = "Diabetes",
    afib = "Atrial fibrillation",
    pad = "Peripheral arterial disease",
    tci = "Previous TIA",
    ami = "Previous MI",
    tpa = "Treated with tPA",
    ais = "Previous AIS",
    evt = "Treated with EVT",
    # reg_any_perf,
    # rtreat = "Study group allocation",
    rtreat_placebo = "Placebo trial treatment",
    pase_0 = "Pre-stroke PASE score",
    pase_0_high = "High pre-stroke PA",
    pase_4 = "6 months post-stroke PASE score",
    # pase_change,
    nihss = "Admission NIHSS",
    # soc_status,
    soc_status_work = "Employed",
    soc_status_nowork = "Not employed",
    fam_indk = "Family income group",
    fam_indk_hl = "Lower family income",
    fam_indk_high = "Higher family income",
    fam_indk_low = "Lower family income",
    edu_level = "Educational level group",
    edu_level_hl = "Lower educational level",
    edu_high = "Higher educational level",
    edu_low = "Low educational level",
    who_4 = "WHO-5 score 6 months post-stroke",
    mdi_4 = "MDI score 6 months post-stroke",
    mrs_4_above1 = "mRS > 1 at 6 months post-stroke",
    mfi_gen_4 = "General fatigue (MFI domain) 6 months post-stroke",
    time = "Time",
    status = "Status",
    event.include = "Include event",
    who_0 = "Pre-stroke WHO-5 score",
    mrs_0_above0 = "Pre-stroke mRS > 0",
    mrs_eos = "mRS at End of Study"
  )
}

#' Subset labels
#'
#' @param data
#' @param labels.raw
#'
#' @return character vector
#' @export
#'
subset_named_labels <- function(data, labels.raw) {
  labels.raw[match(data, names(labels.raw))]
}

#' Assign labels to data.frame or tibble
#'
#' @param data
#' @param labels
#'
#' @return
#' @export
#'
#' @examples
assign_labels <- function(data, labels) {
  # data |> labelled::set_variable_labels(labels)

  labelled::var_label(data) <- labels

  data
}

#' Flexible labelling using labelled for nicer tables
#'
#' @param data data set
#'
#' @return
#' @export labelled data.frame/tibble
#'
#' @examples
#' data <- targets::tar_read(df_pred_data)
#' data <- data |> dplyr::mutate(test = "test")
#' data |>
#'   labelling_data() |>
#'   labelled::var_label()
labelling_data <- function(data, label.list = var_labels()) {
  labs <- subset_named_labels(names(data), label.list)
  labs[is.na(labs)] <- names(data)[is.na(labs)]

  data |> assign_labels(labels = labs)
}


#' Main statistical analysis
#'
#' @param include.pase
#' @param data
#'
#' @return
#' @export
#'
main_analysis <- function(include.pase = FALSE, data) {
  if (!include.pase) data <- dplyr::select(data, -pase_0_high)

  data |>
    # (\(.x){
    #   lm(simple_score ~ .,
    #     data = .x
    #   )
    # })() |>
    (\(.x){
      MASS::polr(simple_score ~ .,
        data = .x, Hess = TRUE,
        method = "logistic"
      )
    })() |>
    gtsummary::tbl_regression(
      show_single_row = dplyr::where(is.logical),
      exponentiate = TRUE
    )
}

coef_forrest_plot <- function(data,
                              cols = c(
                                "#CE0045",
                                "#66c1a3"
                              ),
                              x.tics = create_log_tics(c(.2, .4, .6, .8, 1)),
                              legend.title = "") {
  p1 <- data |>
    # dplyr::filter(name=="decrease") |>
    ggplot2::ggplot(ggplot2::aes(
      x = log(estimate),
      y = label,
      color = model,
      fill = model
    )) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linewidth = .5, linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(
      x = log(estimate), xmin = log(conf.low),
      xmax = log(conf.high)
    ), 
    position = ggplot2::position_dodge(width = 0.5), 
    width = .2,
    linewidth=2,
    lineend = "round")+
    ggplot2::geom_point( # ggplot2::aes(shape = model),
      size = 8,
      shape = 18,
      position = ggplot2::position_dodge(width = 0.5)
    )  +
    ggplot2::scale_x_continuous(
      breaks = log(x.tics),
      labels = x.tics,
      limits = log(range(x.tics))
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    # ggplot2::scale_shape_manual(values = c(1, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      # legend.title = ggplot2::element_text(""),
      legend.position = "bottom"
    ) +
    ggplot2::ylab("") +
    ggplot2::xlab("Odds ratio (log)") +
    ggplot2::labs(
      # shape = legend.title,
      color = legend.title,
      fill = legend.title
    )
  p1
}

poster_coef_print <- function(plot, path) {
  x.tics <- create_log_tics(c(.2, .6, 1))

  ggplot2::ggsave(path,
    plot +
      ggplot2::scale_x_continuous(
        breaks = log(x.tics),
        labels = x.tics,
        limits = log(range(x.tics))
      ) +
      ggplot2::facet_wrap(facets = ggplot2::vars(model), ncol = 2) +
      theme_poster(),
    units = "cm",
    width = 25,
    height = 18
  )
}


#' Theme for minimal plotting for posters and other prints
#'
#' @param plot ggplot object
#' @param text on plot text
#' @param bg plot background
#'
#' @return ggplot object list
#' @export
theme_poster <- function(plot, text=ggplot2::element_text(size = 25) ,bg=ggplot2::element_rect(colour = "white")){
  ggplot2::theme(
    legend.position = "none",
    # panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    text = text,
    plot.title = ggplot2::element_text(),
    panel.background = bg,
    panel.border = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    strip.text.x = ggplot2::element_blank()
  )
}



create_log_tics <- function(data) {
  sort(round(unique(c(1 / data, data)), 2))
}

get_label <- function(vars = "pase_0") {
  subset_named_labels(vars, var_labels())
}

get_set_label <- function(set) {
  get_label(get_var_vec(set))
}
