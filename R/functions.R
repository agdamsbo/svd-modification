# Functions will live here

# Functions from the SVD-burden annotation project
# source("https://raw.githubusercontent.com/agdamsbo/cSVD-burden-annotation-and-analysis-project/main/R/functions.R")

#' Title
#'
#' @param key token
#' @param ... passed on to REDCapCAST::read_redcap_tables()
#'
#' @return
#' @export
#'
#' @examples
#' full_dataset(key = "SVD_REDCAP_API")
full_dataset <- function(key, ...) {
  REDCapCAST::read_redcap_tables(
    uri = "https://redcap.au.dk/api/",
    token = keyring::key_get(key),
    fields = c(
      "record_id",
      "trial_name",
      "trial_id"
    ), forms = c("svd_score", "consensus_score"), ...
  )
}


# REDCapR::redcap_metadata_read(redcap_uri = "https://redcap.au.dk/api/",
#   token = keyring::key_get("SVD_REDCAP_API"))$data

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
#'   View()
#' (\(.x) summary(factor(.x$consensus)))()
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
      "talos_mrs01_0",
      "talos_mrs01_4",
      "rdate",
      "cpr",
      "talos_nihss16_0",
      "rtreat",
      "talos_who07_4",
      "toastclassification"
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
  # browser()
  data_nas <- data |> dplyr::mutate(dplyr::across(
    dplyr::starts_with("talos"),
    ~ dplyr::if_else(.x %in% nas, NA, .x)
  ))

  out <- data_nas |>
    (\(.x){
      # browser()
      .x |>
        dplyr::select(dplyr::starts_with("talos_pase")) |>
        (\(.y){
          .y |> dplyr::select(-grep("(00|12|x)_0$", names(.y)))
        })() |>
        (\(.y){
          ns <- names(.y) |>
            strsplit("_") |>
            lapply(\(.z).z[[2]]) |>
            unlist()
          setNames(.y, ns)
        })() |>
        stRoke::pase_calc() |>
        (\(.y){
          tibble::tibble(.x, .y)
        })()
    })() |>
    dplyr::mutate(
      dob = as.Date(stRoke::cpr_dob(cpr), format = "%d-%m-%Y"),
      age = stRoke::age_calc(dob, as.Date(rdate)),
      inc_date = rdate,
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
        reg_rygning,
        c("1", "2") ~ "current",
        c("3") ~ "prior",
        c("4") ~ "never",
        "9" ~ NA
      ),
      # Living alone defined as not together with somebody
      reg_alone = dplyr::case_match(
        reg_civil, "1" ~ FALSE,
        c("2", "3") ~ TRUE,
        "9" ~ NA
      ),
      toast = forcats::fct_relevel(factor(toastclassification,labels = c("Large artery disease",
                                                                         "Cardioembolic",
                                                                         "Small vessel disease",
                                                                         "Other",
                                                                         "Unknown")),
                                   c("Large artery disease",
                                     "Small vessel disease",
                                     "Cardioembolic",
                                     "Other",
                                     "Unknown")),
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
      pase_0 = pase_score_sum,
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
      inc_date = as.Date(inc_date),
      any_perf,
      nihss = talos_nihss16_0,
      mrs_pre = substr(talos_mrs01_0, 1, 1),
      mrs_eos = substr(talos_mrs01_4, 1, 1),
      who5_eos = as.numeric(talos_who07_4),
      toast
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
      "pre_mrs3",
      "mrs3_final",
      "who5_final",
      "pase_total_score",
      "tpa",
      "evt",
      "cpr",
      "rand_dato_tid",
      "ais_etiology_details",
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
        "pre_mrs3",
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
      inc_date = rand_dato_tid,
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
      # smoker = dplyr::case_match(
      #   smoking, c("1","2", "3") ~ TRUE,
      #   c("4") ~ FALSE,
      #   "9" ~ NA
      # ),
      smoker = dplyr::case_match(
        smoking,
        c("1") ~ "current",
        c("2") ~ "prior",
        c("3") ~ "never",
        "4" ~ NA
      ),
      toast = dplyr::case_match(
        ais_etiology_details,
        c("0") ~ "Large artery disease",
        c("1") ~ "Small vessel disease",
        c("2") ~ "Cardioembolic",
        c("3") ~ "Other",
        c("4") ~ "Unknown",
        .default = NA
      ) |> factor(levels=c("Large artery disease",
                           "Small vessel disease",
                           "Cardioembolic",
                           "Other","Unknown")),
      # Living alone defined as not together with somebody
      alone = dplyr::case_match(
        civil_status, "1" ~ FALSE,
        c("2", "98") ~ TRUE,
        "99" ~ NA
      ),
      alc_more = dplyr::case_match(
        alcohol, "1" ~ FALSE,
        "2" ~ TRUE,
        "99" ~ NA
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
      toast,
      nihss,
      inc_date = as.Date(inc_date),
      tpa,
      evt,
      any_perf,
      mrs_pre = pre_mrs3,
      mrs_eos = mrs3_final,
      who5_eos = as.numeric(who5_final)
    )

  out
}

#' Collect all clinical data
#'
#' @return
#' @export
#'
#' @examples
#' ds <- clin_data_all()
#' ds |> View()
#'   dplyr::mutate(mrs_pre = factor(mrs_pre)) |>
#'   skimr::skim()
clin_data_all <- function() {
  # browser()
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
      pase_0_high = pase_0 > median(as.numeric(pase_0), na.rm = TRUE),
      year = format(inc_date, "%Y"),
      month = format(inc_date, "%m"),
      winter = as.numeric(month) %in% c(1, 2, 12)
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
  gsub(pattern = pattern, replacement = replacement, x = data)
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
      "pase_0",
      "alone",
      "smoker",
      "hyperten",
      "diabetes",
      "ais",
      "afib"
    ),
    out = c(
      "simple_score"
    ),
    preds = c(
      "age",
      "female_sex",
      "pase_0",
      # "alone",
      "smoker",
      # "alc_more",
      "hyperten",
      "diabetes",
      # "tci",
      # "ais"#,
      # "ami",
      # "afib",
      # "pad" # ,
      # "year",
      # "winter"
      "ais_tci"
    ),
    pre = c(
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
      "pad" # ,
      # "year",
      # "winter"
    ),
    clin = c(
      "active_treatment",
      "trial",
      "nihss",
      "tpa",
      "evt"
    ),
    post = c(
      "mrs_eos",
      "who5_eos"
    )
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

fix_labels_raw <- function(data, keep = FALSE) {
  out <- data |>
    subset_named_labels(var_labels()) |>
    unname()

  if (keep) {
    ifelse(is.na(out), data, out)
  } else {
    out
  }
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
  old.labels <- data[[1]][["variable"]][data[[1]][["row_type"]] == "label"]
  new.labels <- old.labels |>
    purrr::map_chr(\(.x){
      .x |>
        strsplit(split = ":") |>
        unlist() |>
        fix_labels_raw() |>
        paste(collapse = " * ")
    }) # The above handles interaction term labels
  data[[1]][["label"]][data[[1]][["row_type"]] == "label"] <- ifelse(is.na(new.labels),
    old.labels,
    new.labels
  )
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

#' SImplifies full score based on supplied cut values
#'
#' @param data
#' @param name.score
#' @param lab.mic
#' @param lab.lac
#' @param lab.wmh
#' @param lab.atr
#' @param append
#'
#' @return
#' @export
#'
#' @examples
simplified_score <- function(data,
                             name.score = "new.simple.score",
                             lab.mic = c(0, 1, 1, 1, 1),
                             lab.lac = c(0, 0, 0, 1, 1),
                             lab.wmh = 0:3,
                             lab.atr = 0:3,
                             append = TRUE) {
  out <- as.numeric(factor(data$microbleed_f, labels = lab.mic)) - 1 +
    as.numeric(factor(data$lacunes_f, labels = lab.lac)) - 1 +
    as.numeric(factor(data$wmh_f, labels = lab.wmh)) - 1 +
    as.numeric(factor(data$atrophy_f, labels = lab.atr)) - 1

  if (append) {
    structure(
      cbind(
        data,
        out
      ),
      class = class(data),
      names = c(colnames(data), name.score)
    )
  } else {
    out
  }
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
    kripp_alpha = FALSE,
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
    ## All possible levels across columns
    lvls <- data |>
      lapply(factor) |>
      lapply(levels) |>
      (\(.x){
        unique(do.call(c, .x))
      })()

    ## Factorised and converted to numeric
    data_fmt <- data |>
      lapply(\(.x){
        factor(.x, levels = lvls) |>
          as.numeric()
      }) |>
      (\(.x){
        do.call(cbind, .x)
      })()

    ## ICC calc
    irr::icc(data_fmt, model = "twoway", type = "agreement", unit = "single") |>
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
      ) |> as.data.frame()
    }) |>
    purrr::map(\(.x){
      icc_calc(na.omit(.x))
    }) |>
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
    smoker = "Smoking",
    never = "Smoking (never)",
    current = "Smoking (current)",
    prior = "Smoking (prior)",
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
    pase_0_q = "Pre-stroke PA quartile",
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
    mrs_eos = "mRS at End of Study",
    mrs_pre = "Pre-stroke mRS",
    active_treatment = "Active treatment",
    trial = "Clinical trial",
    year = "Enrollment year",
    month = "Month of enrollment",
    winter = "Winther enrollment",
    simple_score_f_n = "Simplified score",
    simple_score_f = "Simplified score",
    simple_score = "Simplified score",
    microbleed_f = "Microbleeds subscore",
    lacunes_f = "Lacunes subscore",
    atrophy_f = "Atrophy subscore",
    wmh_f = "WMH subscore",
    full_score_f_n = "Complete score",
    full_score_f = "Complete score",
    full_score = "Complete score",
    full_score_svd = "Complete score",
    ais_tci = "Previous ischemic event",
    olama_amended_ext_svd = "extended amended Olama et al 2020",
    olama_amended_svd = "amended Olama et al 2020",
    olama_simple_ext_svd = "extended simple Olama et al 2020",
    olama_simple_svd = "simple Olama et al 2020",
    huijts_simple_ext_svd = "modified Huijts et al 2013",
    lab.mic = "Microbleeds subscore",
    lab.lac = "Lacunes subscore",
    lab.wmh = "WMH subscore",
    lab.atr = "Atrophy subscore",
    by_individual_best_fit = "Best fit subscores",
    clinical_assumption = "Clinical assumptions",
    svd_score_04 = "SVD score"
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
subset_named_labels <- function(data, labels.raw = var_labels()) {
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
  if (!include.pase) data <- dplyr::select(data, -dplyr::starts_with("pase_0"))

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

# coef_forrest_plot <- function(data,
#                               cols = c(
#                                 "#CE0045",
#                                 "#66c1a3"
#                               ),
#                               x.tics = create_log_tics(c(.2, .4, .6, .8, 1)),
#                               legend.title = "") {
#   p1 <- data |>
#     # dplyr::filter(name=="decrease") |>
#     ggplot2::ggplot(ggplot2::aes(
#       x = log(estimate),
#       y = label,
#       color = model,
#       fill = model
#     )) +
#     ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linewidth = .5, linetype = "dashed") +
#     ggplot2::geom_errorbar(
#       ggplot2::aes(
#         x = log(estimate), xmin = log(conf.low),
#         xmax = log(conf.high)
#       ),
#       position = ggplot2::position_dodge(width = 0.5),
#       width = .2,
#       linewidth = 2,
#       lineend = "round"
#     ) +
#     ggplot2::geom_point( # ggplot2::aes(shape = model),
#       size = 8,
#       shape = 18,
#       position = ggplot2::position_dodge(width = 0.5)
#     ) +
#     ggplot2::scale_x_continuous(
#       breaks = log(x.tics),
#       labels = x.tics,
#       limits = log(range(x.tics))
#     ) +
#     ggplot2::scale_color_manual(values = cols) +
#     ggplot2::scale_fill_manual(values = cols) +
#     # ggplot2::scale_shape_manual(values = c(1, 1)) +
#     ggplot2::theme_bw() +
#     ggplot2::theme(
#       panel.grid.minor = ggplot2::element_blank(),
#       # legend.title = ggplot2::element_text(""),
#       legend.position = "bottom"
#     ) +
#     ggplot2::ylab("") +
#     ggplot2::xlab("Odds ratio (log)") +
#     ggplot2::labs(
#       # shape = legend.title,
#       color = legend.title,
#       fill = legend.title
#     )
#   p1
# }

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
theme_poster <- function(plot, text = ggplot2::element_text(size = 25), bg = ggplot2::element_rect(colour = "white")) {
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

#' Forrest plot based on gtsummary tables
#'
#' @param data data table from gtsummary.
#' @param cols Colors. Will be used if column "model" is present.
#' @param x.tics numeric vector of tics. Use create_log_tics().
#' @param legend.title Title for the legend
#'
#' @return ggplot list object
#' @export
#'
#' @examples
#' stRoke::talos |>
#'   dplyr::mutate(mrs_6 = factor(mrs_6)) |>
#'   dplyr::select(-mrs_1) |>
#'   (\(.x){
#'     MASS::polr(mrs_6 ~ .,
#'       data = .x, Hess = TRUE,
#'       method = "logistic"
#'     )
#'   })() |>
#'   gtsummary::tbl_regression(
#'     exponentiate = TRUE
#'   ) |>
#'   # purrr::pluck("table_body") |>
#'   coef_forrest_plot(x.tics = create_log_tics(c(.2, 1)))
coef_forrest_plot <- function(data,
                              cols = c(
                                "#CE0045",
                                "#66c1a3"
                              ),
                              x.tics = create_log_tics(c(.2, .4, .6, .8, 1)),
                              legend.title = "") {
  if ("tbl_merge" %in% class(data)) {
    data <- pivot_gtsummary_long(data)
  }

  if ("gtsummary" %in% class(data)) {
    data <- purrr::pluck(data, "table_body")
  }

  assertthat::assert_that("tbl_df" %in% class(data))

  ds <- data |>
    (\(.x){
      if (!"model" %in% names(.x)) {
        .x$model <- 1
      }
      split(.x, .x$model)
    })() |>
    purrr::map(function(.x) {
      .x |>
        dplyr::mutate(
          variable = factor(variable, levels = unique(variable))
        ) |>
        (\(.x){
          split(.x, .x$variable)
        })() |>
        purrr::map(function(.x) {
          if (nrow(.x) > 1) {
            .x[-1, ] |>
              dplyr::mutate(label = paste(.x[["label"]][1], label))
          } else {
            .x
          }
        }) |>
        dplyr::bind_rows() |>
        dplyr::mutate(
          estimate = dplyr::if_else(reference_row, 1, estimate, missing = estimate),
          label = factor(label, levels = rev(label))
        )
    }) |>
    dplyr::bind_rows()


  if ("model" %in% names(data)) {
    p1 <- ds |>
      dplyr::mutate(model = factor(model, levels = rev(unique(model)))) |>
      # dplyr::filter(name=="decrease") |>
      ggplot2::ggplot(ggplot2::aes(
        x = log(estimate),
        y = label,
        color = model,
        fill = model
      ))
  } else {
    p1 <- ds |>
      # dplyr::filter(name=="decrease") |>
      ggplot2::ggplot(ggplot2::aes(
        x = log(estimate),
        y = label
      ))
  }

  p1 <- p1 +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linewidth = .5, linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        x = log(estimate), xmin = log(conf.low),
        xmax = log(conf.high)
      ),
      position = ggplot2::position_dodge(width = 0.5),
      width = .2,
      linewidth = 2,
      lineend = "round"
    ) +
    ggplot2::geom_point( # ggplot2::aes(shape = model),
      size = 8,
      shape = 18,
      position = ggplot2::position_dodge(width = 0.5)
    ) +
    ggplot2::scale_x_continuous(
      breaks = log(x.tics),
      labels = x.tics,
      limits = log(range(x.tics))
    ) +
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

  if (is.null(cols)) {
    p1
  } else {
    p1 +
      ggplot2::scale_color_manual(values = cols) +
      ggplot2::scale_fill_manual(values = cols)
  }
}

olr_forrest_plot <- function(data,
                             x.tics = create_log_tics(c(.2, .5, .8, 1)),
                             legend.title = "",
                             cols = NULL,
                             p.size = 8,
                             b.width = .2,
                             l.width = 2,
                             d.with = .5) {
  out <- data |> ggplot2::ggplot(ggplot2::aes(
    x = log(estimate),
    y = label,
    color = model,
    fill = model
  )) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0), linewidth = .5, linetype = "dashed") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        x = log(estimate), xmin = log(conf.low),
        xmax = log(conf.high)
      ),
      position = ggplot2::position_dodge(width = d.with),
      width = b.width,
      linewidth = l.width,
      lineend = "round"
    ) +
    ggplot2::geom_point( # ggplot2::aes(shape = model),
      size = p.size,
      shape = 18,
      position = ggplot2::position_dodge(width = d.with)
    ) +
    ggplot2::scale_x_continuous(
      breaks = log(x.tics),
      labels = x.tics,
      limits = log(range(x.tics))
    ) +
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

  return(out)
}


#' Pivot merged gtsummary table for plotting coefficients
#'
#' @param data list of merged gtsummary tables
#'
#' @return tibble
#' @export
#'
pivot_gtsummary_long <- function(data) {
  assertthat::assert_that("tbl_merge" %in% class(data))

  data |>
    purrr::pluck("table_body") |>
    (\(.x){
      dplyr::select(.x, grep("_\\d+$", names(.x)))
    })() |>
    (\(.x){
      f <- REDCapCAST::strsplitx(names(.x), "\\d+$", "before") |>
        purrr::map(function(.y) .y[2]) |>
        purrr::list_c() |>
        factor()
      split.default(.x, f) |>
        purrr::map(function(.y) setNames(.y, strsplit(names(.y), "_\\d+$"))) |>
        setNames(names(purrr::pluck(data, "tbls"))) |>
        purrr::imap(function(.z, .i) dplyr::mutate(.z, model = .i)) |>
        purrr::map(function(.n) {
          dplyr::bind_cols(
            data |>
              purrr::pluck("table_body") |>
              (\(.m){
                dplyr::select(.m, !grep("_\\d+$", names(.m)))
              })(),
            .n
          )
        }) |>
        dplyr::bind_rows()
    })()
}

#' Functional run for multiple OLR models
#'
#' @param data
#' @param out
#' @param gen_exp
#' @param adjust
#' @param main_exp
#' @param include.adjust
#'
#' @return list
#' @export

multi_olr <- function(data, out, gen_exp, adjust, main_exp, include.adjust = FALSE) {
  if (include.adjust) {
    all.exp <- c(main_exp, adjust, gen_exp)
  } else {
    all.exp <- c(main_exp, gen_exp)
  }

  reg.formula <- glue::glue("{out}~{paste0(unique(all.exp), collapse = '+')}")

  data |>
    (\(.x){
      MASS::polr(as.formula(reg.formula),
        data = .x,
        Hess = TRUE,
        method = "logistic"
      )
    })()
}

multi_olr_glue <- function(data,
                           out,
                           gen_exp = NULL,
                           pase,
                           main_exp,
                           formula.glue = "{out}~{main_exp}*{pase}{ifelse(!is.null(gen_exp),paste('+',paste0(unique(gen_exp), collapse = '+')),'')}") {
  reg.formula <- glue::glue(formula.glue)

  data |>
    (\(.x){
      MASS::polr(as.formula(reg.formula),
        data = .x,
        Hess = TRUE,
        method = "logistic"
      )
    })()
}

#' Standard regression table function
#'
#' @param data model
#'
#' @return gtsummary list
#' @export
#'
reg_table <- function(data) {
  data |>
    gtsummary::tbl_regression(
      show_single_row = dplyr::where(is.logical),
      exponentiate = TRUE
    )
}

#' Adding interaction terms results in gtsummary not being able to only show single row for logicals.
#' This fundtion handles the functionality in a crude but efficient way.
#'
#' @param data gtsummary table list
#'
#' @return
#' @export
#'
# out <-
#' @examples
compress_logical_vars <- function(data) {
  data[[1]] <- split(data[[1]], data[[1]][["variable"]]) |>
    purrr::map(\(.d){
      # browser(skipCalls = 2)
      if (!is.na(unique(.d$var_class)) & unique(.d$var_class) == "logical") {
        n <- match("n_obs", names(.d)) # Used si split
        out <- dplyr::bind_cols(
          dplyr::select(.d[1, ], 1:n - 1),
          dplyr::select(.d[nrow(.d), ], n:ncol(.d))
        )
      } else {
        out <- .d
      }
      out
    }) |>
    dplyr::bind_rows()

  data
}


save_table <- function(data, filename) {
  data |>
    gtsummary::as_gt() |>
    gt::gtsave(filename = filename)
}


#' Formatting main data set for analyses
#'
#' @param data
#'
format_df <- function(data) {
  # browser()
  data |>
    dplyr::mutate(
      ais_tci = ais | tci,
      smoker = factor(smoker, levels = c(
        "never",
        "current",
        "prior"
      )),
      microbleed_f = factor(microbleed, levels = c("0", "1", "2-4", "5-10", ">10")),
      lacunes_f = factor(lacunes, levels = c("0", "1", "2", "3-5", ">5")),
      dplyr::across(tidyselect::all_of(c("atrophy", "wmh")),
        ~ as.factor(.x),
        .names = "{.col}_f"
      )
    ) |>
    dplyr::mutate(dplyr::across(c(microbleed_f, lacunes_f, atrophy_f, wmh_f),
      ~ as.numeric(.x) - 1,
      .names = "{.col}_n"
    )) |>
    dplyr::mutate(full_score = dplyr::pick(
      microbleed_f_n,
      lacunes_f_n,
      atrophy_f_n,
      wmh_f_n
    ) |>
      rowSums()) |>
    dplyr::mutate(
      dplyr::across(c(full_score, simple_score),
        ~ as.factor(.x),
        .names = "{.col}_f"
      )
    ) |>
    dplyr::mutate(dplyr::across(c(simple_score_f, full_score_f),
      ~ as.numeric(.x) - 1,
      .names = "{.col}_n"
    )) |>
    dplyr::mutate(
      pase_0_q = stRoke::quantile_cut(pase_0, groups = 4, group.names = paste0("Q", 1:4)),
      pase_0_q = forcats::fct_rev(pase_0_q)
    )
}

svd_scoring_systems <- function() {
  list(
    "olama_amended_ext_svd" =
      list(
        lab.mic = c(0, 1, 1, 1, 1),
        lab.lac = c(0, 0, 1, 2, 3),
        lab.wmh = 0:3,
        lab.atr = 0:3
      ),
    "olama_amended_svd" = list(
      lab.mic = c(0, 1, 1, 1, 1),
      lab.lac = c(0, 0, 1, 2, 3),
      lab.wmh = 0:3,
      lab.atr = rep(0, 4)
    ),
    "olama_simple_ext_svd" = list(
      lab.mic = c(0, 1, 1, 1, 1),
      lab.lac = c(0, 0, 0, 1, 1),
      lab.wmh = c(0, 0, 1, 1),
      lab.atr = c(0, 0, 1, 1)
    ),
    "olama_simple_svd" = list(
      lab.mic = c(0, 1, 1, 1, 1),
      lab.lac = c(0, 0, 0, 1, 1),
      lab.wmh = c(0, 0, 1, 1),
      lab.atr = rep(0, 4)
    ),
    "huijts_simple_ext_svd" = list(
      lab.mic = c(0, 1, 1, 1, 1),
      lab.lac = c(0, 1, 1, 1, 1),
      lab.wmh = c(0, 0, 1, 1),
      lab.atr = c(0, 0, 1, 1)
    ),
    "full_score_svd" = list(
      lab.mic = 0:4,
      lab.lac = 0:4,
      lab.wmh = 0:3,
      lab.atr = 0:3
    ),
    "by_individual_best_fit" = list(
      lab.mic = c(0, 0, 0, 0, 1),
      lab.lac = c(0, 0, 0, 0, 1),
      lab.wmh = c(0, 1, 1, 1),
      lab.atr = c(0, 0, 1, 1)
    ),
    "clinical_assumption" = list(
      lab.mic = c(0, 0, 1, 1, 2),
      lab.lac = c(0, 1, 1, 1, 2),
      lab.wmh = c(0, 0, 1, 1),
      lab.atr = c(0, 0, 1, 1)
    )
  )
}

microbleed_cuts <- function() {
  list(
    "mic_0_some" =
      list(
        lab.mic = c(0, 1, 1, 1, 1),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "mic_0_1_more_many" =
      list(
        lab.mic = c(0, 1, 2, 2, 3),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "mic_0_1_more" =
      list(
        lab.mic = c(0, 1, 2, 2, 2),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "mic_few_more" =
      list(
        lab.mic = c(0, 0, 1, 1, 1),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "mic_some_lots" =
      list(
        lab.mic = c(0, 0, 0, 0, 1),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      )
  )
}

lacunes_cuts <- function() {
  list(
    "lac_0_some" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = c(0, 1, 1, 1, 1),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "lac_0_1_more_many" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = c(0, 1, 2, 2, 3),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "lac_0_1_more" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = c(0, 1, 2, 2, 2),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "lac_few_more" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = c(0, 0, 1, 1, 1),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      ),
    "lac_some_lots" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = c(0, 0, 0, 0, 1),
        lab.wmh = rep(0, 4),
        lab.atr = rep(0, 4)
      )
  )
}

wmh_cuts <- function() {
  list(
    "wmh_0_some" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = c(0, 1, 1, 1),
        lab.atr = rep(0, 4)
      ),
    "wmh_little_some" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = c(0, 0, 1, 1),
        lab.atr = rep(0, 4)
      ),
    "wmh_all" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = c(0, 1, 2, 3),
        lab.atr = rep(0, 4)
      ),
    "wmh_some_much" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = c(0, 1, 1, 2),
        lab.atr = rep(0, 4)
      )
  )
}

atr_cuts <- function() {
  list(
    "atr_0_some" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = c(0, 1, 1, 1)
      ),
    "atr_little_some" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = c(0, 0, 1, 1)
      ),
    "atr_all" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = c(0, 1, 2, 3)
      ),
    "atr_some_much" =
      list(
        lab.mic = rep(0, 5),
        lab.lac = rep(0, 5),
        lab.wmh = rep(0, 4),
        lab.atr = c(0, 1, 1, 2)
      )
  )
}

svdscore_analysis_compare <- function(data, list) {
  list |> purrr::imap(\(.z, .i){
    df_subs <- data |> svd_system_calc(list.system = .z)
    preds <- c("pase_0", purrr::pluck(model_input(), "main_exp"))
    # browser()
    purrr::map(preds, \(.y){
      vars <- names(df_subs) |>
        tail(length(.z))

      vars |>
        purrr::map(\(.x){
          .df <- df_subs
          .df[.x] <- factor(.df[[.x]])
          reg.formula <- paste0(.x, "~", paste0(c(.y), collapse = "+"))

          if (length(levels(.df[[.x]])) == 2) {
            out <- glm(as.formula(reg.formula),
              family = "binomial",
              data = .df
            )
          } else {
            out <- MASS::polr(as.formula(reg.formula),
              data = .df,
              Hess = TRUE,
              method = "logistic"
            )
          }
          out
        }) |>
        setNames(vars)
    }) |>
      setNames(preds) |>
      purrr::map(performance::compare_performance, rank = TRUE)
  })
}

svd_system_table <- function(data,
                             system = "new_score",
                             lab.mic = c(0, 1, 1, 1, 1),
                             lab.lac = c(0, 0, 0, 1, 1),
                             lab.wmh = 0:3,
                             lab.atr = 0:3) {
  list(
    microbleed_f = tibble::tibble(levels(data$microbleed_f), lab.mic),
    lacunes_f = tibble::tibble(levels(data$lacunes_f), lab.lac),
    wmh_f = tibble::tibble(levels(data$wmh_f), lab.wmh),
    atrophy_f = tibble::tibble(levels(data$atrophy_f), lab.atr)
  ) |>
    purrr::map(setNames, c("Annotation", system)) |>
    dplyr::bind_rows()
}

# |>
#   gt::gt()
# |>
#   purrr::map(gt::gt) |>
#   purrr::imap(\(.x,.i){
#     # browser()
#     .x |> gt::tab_row_group(
#       label = fix_labels_raw(.i),
#       rows = seq_len(nrow(.x$`_data`))
#     )|>
#       purrr::pluck(1)
#   })



svd_systems_overview <- function(data, list.system = svd_scoring_systems(), subset = seq_along(svd_scoring_systems())) {
  # browser()
  out <- purrr::imap(list.system[subset], \(.x, .i){
    params <- list(.x, "data" = data, "system" = .i) |>
      purrr::list_flatten()
    do.call(svd_system_table, params)
  }) |>
    purrr::reduce(function(d1, d2) {
      dplyr::bind_cols(d1, d2[, !names(d2) %in% names(d1)])
    })

  names(out) <- fix_labels_raw(names(out), keep = TRUE)

  list <- list.system |>
    purrr::pluck(1) |>
    setNames(nm = NULL)
  list_seqs <- list |> purrr::imap(\(.x, .i){
    seq_along(.x) + sum(lengths(list[seq_len(.i - 1)]))
  })

  out <- out |>
    gt::gt()

  for (i in rev(seq_len(length(list.system |> purrr::pluck(1))))) {
    out <- out |>
      gt::tab_row_group(
        label = fix_labels_raw(names(purrr::pluck(list.system, 1))[i]),
        rows = list_seqs[[i]]
      )
  }

  out
}

#' Appends SVD scoring systems with predefined arguments from list of lists
#'
#' @param data data containing relevant variables
#'
#' @return
#' @export
#'
#' @examples
svd_system_calc <- function(data, list.system = svd_scoring_systems()) {
  message(paste("Added scores:", paste(names(list.system), collapse = ", ")))

  dplyr::bind_cols(
    data,
    purrr::map(list.system, \(.x){
      params <- list(.x, "data" = data, "append" = FALSE) |>
        purrr::list_flatten()
      do.call(simplified_score, params)
    }) |> dplyr::bind_cols()
  )
}

#' Vertical stacked bar plot wrapper
#'
#' @param data
#' @param score
#' @param group
#' @param strata
#' @param t.size
#'
#' @return
#' @export
#'
#' @examples
vertical_stacked_bars <- function(data,
                                  score = "full_score",
                                  group = "pase_0_q",
                                  strata = NULL,
                                  t.size = 10,
                                  l.color = "black",
                                  l.size = .5,
                                  draw.lines = FALSE) {
  df.table <- data |>
    dplyr::select(
      {{ score }},
      {{ group }},
      {{ strata }}
    ) |>
    table()

  p <- df.table |>
    rankinPlot::grottaBar(
      scoreName = score,
      groupName = group,
      textColor = c("black", "white"),
      strataName = strata,
      textCut = 6,
      textSize = 20,
      printNumbers = "none",
      lineColor = l.color,
      lineSize = l.size,
      drawLines = draw.lines,
      returnData = TRUE
    )

  colors <- viridisLite::viridis(nrow(df.table))
  contrast_cut <-
    sum(project.aid::contrast_text(colors, threshold = .3) == "white")

  p |>
    (\(.x){
      .x$plot +
        ggplot2::geom_text(
          data = .x$rectData[which(.x$rectData$n >
            0), ],
          size = t.size,
          fontface = "plain",
          ggplot2::aes(
            x = group,
            y = p_prev + 0.49 * p,
            color = as.numeric(score) > contrast_cut,
            # label = paste0(sprintf("%2.0f", 100 * p),"%"),
            label = sprintf("%2.0f", 100 * p)
          )
        ) +
        ggplot2::labs(fill = "Score") +
        ggplot2::scale_fill_manual(values = rev(colors)) +
        ggplot2::theme(
          legend.position = "bottom",
          axis.title = ggplot2::element_text(),
        ) + ggplot2::xlab("Physical activity score quartile") + ggplot2::ylab(NULL)
      # viridis::scale_fill_viridis(discrete = TRUE, direction = -1, option = "D")
    })()
}

multivariable_olr <- function(data, out) {
  MASS::polr(as.formula(glue::glue("{out}~.")),
    data = data, Hess = TRUE,
    method = "logistic"
  )
}

regression_all_olr <- function(data, out) {
  list(
    "Univariable" = data |>
      gtsummary::tbl_uvregression(
        method = MASS::polr,
        method.args = list(
          Hess = TRUE,
          method = "logistic"
        ),
        y = !!dplyr::sym(out),
        show_single_row = dplyr::where(is.logical),
        exponentiate = TRUE # ,
        # pvalue_fun = ~ style_pvalue(.x, digits = 2)
      ),
    "Multivariable" = multivariable_olr(
      data = data,
      out = out
    ) |>
      gtsummary::tbl_regression(
        show_single_row = dplyr::where(is.logical), exponentiate = TRUE
      )
  )
}


models_uni_mini_multi_olr <- function(data,
                                      out,
                                      mini = c(
                                        "age",
                                        # "mrs_pre",
                                        "female_sex"
                                      )) {
  list(
    "Univariable" = list(data = dplyr::select(
      data,
      tidyselect::any_of(
        c(
          out,
          "pase_0_q"
        )
      )
    )),
    "Minimal" = list(data = dplyr::select(
      data,
      tidyselect::any_of(
        c(
          out,
          "pase_0_q",
          mini
        )
      )
    )),
    "Multivariable" = list(data = data)
  ) |>
    lapply(\(.x){
      do.call(
        multivariable_olr,
        c(
          .x,
          list(out = out)
        )
      )
    })
}

regression_uni_mini_multi_olr <- function(...) {
  models_uni_mini_multi_olr(...) |>
    lapply(\(.x){
      .x |>
        gtsummary::tbl_regression(
          show_single_row = dplyr::where(is.logical),
          exponentiate = TRUE
        ) |>
        gtsummary::add_glance_table()
    })
}


svd_dist <- function(data) {
  dplyr::as_tibble(x = data) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram()
}


model_input <- function() {
  list(
    # gen_exp <- c("age", "female_sex"),
    # gen_exp <- c("age"),
    gen_exp = NULL, # Considered outcome variables/subscores
    outcomes = c("simple_score_f", "full_score_f", "microbleed_f", "lacunes_f", "wmh_f", "atrophy_f"),
    # Considered main exposures
    main_exp = c(
      "age",
      "female_sex",
      "alone",
      "smoker",
      "alc_more",
      "hyperten",
      "diabetes",
      "ais_tci",
      "afib",
      "ami",
      "mrs_pre"
    ),
    # Adjustment variable name
    adjust = "pase_0_q"
  )
}



#' Title
#'
#' @param data
#' @param outcomes
#' @param main_exp
#' @param gen_exp
#' @param adjust
#'
#' @return
#' @export
#'
#' @examples
#'
#' planned_multi_olr
#'
planned_multi_olr <- function(data,
                              outcomes = purrr::pluck(model_input(), "outcomes"),
                              main_exp = purrr::pluck(model_input(), "main_exp"),
                              gen_exp = purrr::pluck(model_input(), "gen_exp"),
                              adjust = purrr::pluck(model_input(), "adjust")) {
  complete_model_list <- outcomes |>
    # purrr::set_names() |>
    purrr::map(\(.o){
      main_exp |>
        # purrr::set_names() |>
        purrr::map(\(.z){
          c("No PA adjustment" = FALSE, "With PA adjustment" = TRUE) |>
            (\(.y){
              .y |>
                # purrr::set_names() |>
                purrr::map(\(.x){
                  multi_olr(data,
                    out = .o,
                    gen_exp = gen_exp,
                    adjust = adjust,
                    main_exp = .z,
                    include.adjust = .x
                  )
                })
            })()
        })
    })

  complete_list <- complete_model_list |>
    purrr::map(\(.x){
      purrr::map(.x, \(.y){
        .y |>
          purrr::map(reg_table) |>
          gtsummary::tbl_merge(tab_spanner = names(.y)) |>
          gtsummary::remove_row_type(
            variables = tidyselect::any_of(c(adjust, gen_exp)),
            type = "all"
          )
      }) |> gtsummary::tbl_stack()
    }) |>
    setNames(fix_labels_raw(outcomes))

  regression_stack <- complete_list |>
    (\(.x){
      .x |> gtsummary::tbl_stack(group_header = names(.x))
    })() |>
    fix_labels()

  list(
    "complete_model_list" = complete_model_list,
    "complete_list" = complete_list,
    "regression_stack" = regression_stack
  )
}

#' Title
#'
#' @param data
#'
#' @returns
#' @export
#'
#' @examples
#' formatted_df(targets::tar_read(df_complete))
#'
formatted_df <- function(data) {
  data |>
    dplyr::filter(include) |>
    format_df() |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c("afib", "ami", "ais_tci", "hyperten", "diabetes")), \(.x){
      # Thsi only concerns one patient. Assumption is that if nothing is known, it is not present
      .x[is.na(.x)] <- FALSE
      .x
    })) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c("mrs_pre")), \(.x){
      # This concerns one patient only, meidan is used to not affect estimates
      .x[is.na(.x)] <- median(.x, na.rm = TRUE)
      .x
    }))
}

stack2long <- function(data) {
  purrr::pluck(data, "regression_stack") |>
    purrr::pluck("table_body") |>
    (\(.x){
      purrr::map(c(1, 2), \(.y){
        # browser(skipCalls = 1)
        .x[c(1:7, grep(pattern = paste0("_", .y, "$"), names(.x)))] |>
          dplyr::mutate(model = .y)
      })
    })() |>
    purrr::map(\(.y){
      names(.y) <- gsub("_[12]$", "", names(.y))
      .y
    }) |>
    dplyr::bind_rows()
}

format_long <- function(data,
                        model.labels = c(
                          "No PA adjustment",
                          "With PA adjustment"
                        )) {
  data |>
    dplyr::mutate(
      groupname_col = factor(groupname_col,
        levels = unique(groupname_col)
      ),
      model = factor(model,
        labels = model.labels
      ),
      ## Fixing label labels by using previous function to name levels
      label = ifelse(row_type == "level",
        glue::glue("{fix_labels_raw(var_label)} ({label})"),
        fix_labels_raw(var_label)
      ),
      label = forcats::as_factor(label)
    ) |>
    ## Filters out reference and header rows
    dplyr::filter(
      reference_row == FALSE &
        header_row == FALSE | is.na(header_row)
    )
}

#' Title
#'
#' @param data
#' @param model.labels
#'
#' @returns
#' @export
#'
#' @examples
#' targets::tar_read(lst_multi_olr) |> multi_coef_plot()
#'
multi_coef_plot <- function(data,
                            model.labels = c(
                              "No PA adjustment",
                              "With PA adjustment"
                            )) {
  complete_long <- data |> stack2long()

  complete_long_format <- complete_long |> format_long(model.labels)

  p1 <- complete_long_format |>
    dplyr::filter(groupname_col %in% c("Simplified score", "Complete score")) |>
    olr_forrest_plot() +
    ggplot2::facet_wrap(ggplot2::vars(groupname_col)) +
    ggplot2::xlab("Odds ratio of higher SVD score (log)") +
    ggplot2::ylab("Possible effect modifier")



  # ggplot2::ggsave(
  #   filename = "olr_plot_simple_score.png",
  #   plot = p1, width = 140, height = 140, units = "mm", dpi = 600
  # )

  p2 <- complete_long_format |>
    dplyr::filter(!groupname_col %in% c("Simplified score", "Complete score")) |>
    olr_forrest_plot() +
    ggplot2::facet_wrap(ggplot2::vars(groupname_col)) +
    ggplot2::xlab("Odds ratio of higher SVD sub-score (log)") +
    ggplot2::ylab("Possible effect modifier")

  list(p1, p2)
}

write_collinearity_test_rds <- function(data, file) {
  # check_orig <- targets::tar_read(lst_multi_olr_orig)[[1]][[1]] |> purrr::map(\(.x){
  #   performance::check_collinearity(.x[[2]], check = "vif")
  # })

  check_orig <- data[[1]][[1]] |> purrr::map(\(.x){
    performance::check_collinearity(.x[[2]], check = "vif")
  })
  # here::here("check_orig.rds")
  saveRDS(check_orig, file = file)


  # check_reduced <- targets::tar_read(lst_multi_olr)[[1]][[1]] |> purrr::map(\(.x){
  #   performance::check_collinearity(.x[[2]], check = "vif")
  # })

  # check_reduced <- ds2[[1]][[1]] |> purrr::map(\(.x){
  #   performance::check_collinearity(.x[[2]], check = "vif")
  # })
  #
  # saveRDS(check_reduced,file = here::here("check_reduced.rds"))

  # check_orig |>
  #   purrr::map(plot) |>
  #   patchwork::wrap_plots(ncol = 3)
  # check_reduced |>
  #   purrr::map(plot) |>
  #   patchwork::wrap_plots(ncol = 3)
}

tbl_merge_named <- function(data) {
  gtsummary::tbl_merge(tbls = data, tab_spanner = names(data))
}

gt_theme <- function(data, ...) {
  data |>
    gt::tab_options(
      row.striping.background_color = "#fafafa",
      table_body.hlines.color = "#f6f7f7",
      source_notes.font.size = 12,
      table.font.size = 12,
      table.width = gt::px(700),
      heading.align = "left",
      heading.title.font.size = 24,
      table.border.top.color = "transparent",
      table.border.top.width = gt::px(3),
      data_row.padding = gt::px(7),
      # row_group.font.size = 16,
      # row_group.background.color = "grey"
      ...
    ) |>
    gt::opt_row_striping() |>
    gt::opt_all_caps() |>
    gt::opt_table_font(
      font = list(
        gt::google_font("opensans"),
        gt::default_fonts()
      )
    )
}


planned_multi_olr <- function(data,
                              outcomes = purrr::pluck(model_input(), "outcomes"),
                              main_exp = purrr::pluck(model_input(), "main_exp"),
                              gen_exp = purrr::pluck(model_input(), "gen_exp"),
                              adjust = purrr::pluck(model_input(), "adjust")) {
  complete_model_list <- outcomes |>
    purrr::map(\(.o){
      main_exp |>
        purrr::map(\(.z){
          c("No PA adjustment" = FALSE, "With PA adjustment" = TRUE) |>
            (\(.y){
              .y |>
                purrr::map(\(.x){
                  multi_olr(data,
                    out = .o,
                    gen_exp = gen_exp,
                    adjust = adjust,
                    main_exp = .z,
                    include.adjust = .x
                  )
                })
            })()
        })
    })

  complete_list <- complete_model_list |>
    purrr::map(\(.x){
      purrr::map(.x, \(.y){
        .y |>
          purrr::map(reg_table) |>
          gtsummary::tbl_merge(tab_spanner = names(.y)) |>
          gtsummary::remove_row_type(
            variables = tidyselect::any_of(c(adjust, gen_exp)),
            type = "all"
          )
      }) |> gtsummary::tbl_stack()
    }) |>
    setNames(fix_labels_raw(outcomes))

  regression_stack <- complete_list |>
    (\(.x){
      .x |> gtsummary::tbl_stack(group_header = names(.x))
    })() |>
    fix_labels()

  list(
    "complete_model_list" = complete_model_list,
    "complete_list" = complete_list,
    "regression_stack" = regression_stack
  )
}


#' Title
#'
#' @param data
#'
#' @returns
#' @export
#'
#' @examples
#' targets::tar_read("df_formatted") |> mulitvar_olr()
#' targets::tar_read("df_formatted") |> mulitvar_olr(reg.fun = regression_uni_mini_multi_olr)
mulitvar_olr <- function(data, reg.fun = "regression_all_olr") {
  # regression_uni_mini_multi_olr
  ds <- data |>
    svd_system_calc() |>
    dplyr::select(tidyselect::all_of(c(
      purrr::pluck(model_input(), "gen_exp"),
      purrr::pluck(model_input(), "main_exp"),
      "pase_0_q",
      "huijts_simple_ext_svd"
    ))) |>
    dplyr::mutate(
      huijts_simple_ext_svd = factor(huijts_simple_ext_svd) # ,
      # pase_0_q = forcats::fct_rev(pase_0_q)
    )

  do.call(
    what = reg.fun,
    args = list(
      data = ds,
      out = "huijts_simple_ext_svd"
    )
  ) |>
    purrr::map(fix_labels) |>
    tbl_merge_named()

  # regression_all_olr(out = "huijts_simple_ext_svd") |>
  # # purrr::pluck("Univariable") |>
  # purrr::map(fix_labels) |>
  # tbl_merge_named()
}



mulitvar_olr_interact <- function(data) {
  data |>
    svd_system_calc() |>
    dplyr::mutate(
      huijts_simple_ext_svd = factor(huijts_simple_ext_svd) # ,
      # pase_0_q = forcats::fct_rev(pase_0_q)
    ) |>
    (\(.x){
      vars <- purrr::pluck(model_input(), "main_exp")
      vars |>
        purrr::map(\(.y){
          multi_olr_glue(
            data = .x,
            out = "huijts_simple_ext_svd",
            pase = "pase_0_q",
            gen_exp = c("age", "female_sex"),
            main_exp = .y
          ) |>
            # reg_table() |>
            gtsummary::tbl_regression(exponentiate = TRUE)
        }) |>
        setNames(vars)
    })()
}

#' Multi list model merger
#'
#' @param data1
#' @param data2
#'
#' @return
#' @export
#'
#' @examples
#' multi_var_interaction_merge(data1 = targets::tar_read(lst_multi_olr_orig), data2 = targets::tar_read(list_multi_var_olr_interact))
multi_var_interaction_merge <- function(data1, data2) {
  # browser()
  require(gtsummary)
  d1_format <- data1[[1]] |>
    purrr::pluck(match("simple_score_f", purrr::pluck(model_input(), "outcomes"))) #|>

  d1_format |>
    purrr::imap(\(.y, .i){
      .y |>
        purrr::map(reg_table)
    })

  list_merge <- d1_format |>
    purrr::imap(\(.x, .i){
      list(.x,
        "With PA interaction" = data2[[.i]] |>
          compress_logical_vars()
      ) |>
        purrr::list_flatten()
    }) |>
    setNames(purrr::pluck(model_input(), "main_exp"))

  gt_merge <- list_merge |>
    purrr::imap(\(.x, .i){
      .x |>
        tbl_merge_named() |>
        gtsummary::remove_row_type(
          variables = -dplyr::one_of(.i),
          type = "all"
        )
    })

  gt_merge |>
    gtsummary::tbl_stack() |>
    fix_labels()
}


# |>
#   gtsummary::tbl_regression(
#     show_single_row = dplyr::where(is.logical),
#     exponentiate = TRUE
#   )

#' Title
#'
#' @param data
#' @param out
#' @param gen_exp
#' @param pase
#' @param main_exp
#' @param formula.strings
#'
#' @return
#' @export
#'
#' @examples
#' multi_string_olr(data = df_format, out = "simple_score_f", pase = "pase_0_q", gen_exp = c("age", "female_sex"), main_exp = purrr::pluck(model_input(), "main_exp"))
multi_string_olr <- function(data,
                             out,
                             gen_exp = NULL,
                             pase,
                             main_exp,
                             formula.strings = c(
                               "No PA adjustment" = "{out}~{.exp}",
                               "With PA adjustment" = "{out}~{.exp}+{pase}",
                               "With PA interaction" = "{out}~{.exp}*{pase}"
                             )) {
  gen_adj <- ifelse(!is.null(gen_exp), paste0("+", paste(unique(gen_exp), collapse = "+")), "")
  main_exp <- main_exp[!main_exp %in% gen_exp]

  main_exp |>
    purrr::set_names() |>
    purrr::map(\(.exp){
      formula.strings |>
        purrr::map(\(.s){
          polr_map(data,
            formula = formula(glue::glue(.s, "{gen_adj}")),
            Hess = TRUE,
            method = "logistic"
          )
        }) |>
        setNames(names(formula.strings))
    })
}


#' Convenience function to insert correct formula in model output for mapping
#'
#' @param data data
#' @param formula formula
#' @param ... passed on to MASS::polr
#'
#' @return
#' @export
#'
#' @examples
polr_map <- function(data,
                     formula,
                     ...) {
  m <- MASS::polr(
    formula = formula,
    data = data,
    ...
  )
  m$call$formula <- formula

  m
}


mini_multi_merge <- function(data) {
  data |>
    purrr::imap(\(.x, .i) {
      purrr::map(.x, \(.y) {
        gtsummary::tbl_regression(.y, exponentiate = TRUE) |>
          compress_logical_vars()
      }) |>
        tbl_merge_named() |>
        gtsummary::remove_row_type(
          variables = -dplyr::one_of(.i),
          type = "all"
        )
    }) |>
    gtsummary::tbl_stack() |>
    fix_labels()
}


######
###### DAGs
######

dags_list <- function() {
  list(
    mediator = ggdag::dagify(
      SVD ~ PA + CVRF,
      PA ~ AGE_SEX + CVRF,
      CVRF ~ AGE_SEX,
      outcome = "SVD",
      exposure = c("CVRF", "PA"),
      coords = dag_coords()
    ),
    confounder = ggdag::dagify(
      SVD ~ PA + CVRF,
      PA ~ AGE_SEX,
      CVRF ~ AGE_SEX + PA,
      outcome = "SVD",
      exposure = c("CVRF", "PA"),
      coords = dag_coords()
    ),
    interaction = ggdag::dagify(
      SVD ~ PA + CVRF + PA_CVRF,
      PA ~ AGE_SEX,
      CVRF ~ AGE_SEX,
      PA_CVRF ~ PA + CVRF,
      outcome = "SVD",
      exposure = c("CVRF", "PA"),
      coords = dag_coords()
    ),
    modification = ggdag::dagify(
      SVD ~ CVRF + PA_CVRF,
      PA ~ AGE_SEX,
      CVRF ~ AGE_SEX,
      PA_CVRF ~ PA + CVRF,
      outcome = "SVD",
      exposure = c("CVRF", "PA", "PA_CVRF"),
      coords = dag_coords()
    ),
    independent = ggdag::dagify(
      SVD ~ PA + CVRF,
      PA ~ AGE_SEX,
      CVRF ~ AGE_SEX,
      outcome = "SVD",
      exposure = c("CVRF", "PA"),
      coords = dag_coords()
    ),
    intermediate_pa = ggdag::dagify(
      SVD ~ CVRF,
      PA ~ AGE_SEX,
      CVRF ~ PA + AGE_SEX,
      outcome = "SVD",
      exposure = c("PA"),
      coords = dag_coords()
    ),
    intermediate_cvrf = ggdag::dagify(
      SVD ~ PA,
      PA ~ AGE_SEX + CVRF,
      CVRF ~ AGE_SEX,
      outcome = "SVD",
      exposure = c("CVRF"),
      coords = dag_coords()
    )
  )
}


dag_labels <- function() {
  c(
    AGE_SEX = "Age +\n sex",
    "AGE+SEX" = "Age and sex",
    PA = "Physical\nactivity",
    CVRFs = "CVRF",
    CVRF = "Risk\nfactor",
    "SVD burden" = "SVD burden",
    SVD = "SVD\nburden",
    "PA*CVRFs" = "Interaction effect",
    "PA_CVRF" = "Interaction effect"
  )
}

keyed_dict_lookup <- function(data, dict = dag_labels()) {
  unname(dict[match(data, names(dict))])
}

add_dag_labels <- function(data) {
  data$data$label <- keyed_dict_lookup(data$data$name, dict = dag_labels())

  data
}

dag_colors <- function() {
  c(
    AGE_SEX = "grey65",
    PA = "#3B528BFF",
    CVRF = "#5DC863FF",
    SVD = "#440154FF",
    PA_CVRF = "#21908CFF"
  )
}

dag_coords <- function() {
  list(
    x = c(SVD = 4, CVRF = 2, PA = 3, AGE_SEX = 1, PA_CVRF = 3),
    y = c(SVD = 1, CVRF = 1, PA = sqrt(3) + 1, AGE_SEX = sqrt(3) + 1, PA_CVRF = sqrt(3) / 3 + 1)
  )
}

dags_names <- function() {
  c(
    mediator = "Mediator",
    confounder = "Confounder",
    interaction = "Effect modification",
    modification = "Effect modification",
    independent = "Independent",
    intermediate_pa = "PA is intermediate",
    intermediate_cvrf = "CVRF is intermediate"
  )
}

dags_plot <- function(data, plot.name = NULL) {
  if (!is.null(plot.name)) {
    plot.name <- keyed_dict_lookup(plot.name, dags_names())
  }
  data |>
    ggdag::tidy_dagitty() |>
    add_dag_labels() |>
    dplyr::mutate(
      name = factor(name),
      adjusts = factor(ifelse(name == "AGE_SEX", "adjusted", "unadjusted"))
    ) |>
    # node_dconnected(from = c("PA", "CVRF"), controlling_for = "AGE_SEX") |>
    ggplot2::ggplot(ggplot2::aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      shape = adjusts,
      col = name
    )) +
    ggdag::geom_dag_edges(end_cap = ggraph::circle(7, "mm")) +
    # geom_dag_collider_edges() +
    # ggdag::geom_dag_point(size = 8) +
    # ggdag::geom_dag_text(col = "white") +
    ggdag::geom_dag_label(
      ggplot2::aes(
        label = label,
        fill = name
      ),
      col = "white",
      size = 3
      # label.padding = .15
    ) +
    ggdag::theme_dag() +
    ggdag::scale_adjusted() +
    ggdag::expand_plot(expand_y = ggplot2::expansion(c(0.2, 0.2))) +
    # ggplot2::ggtitle(plot.name)+
    ggplot2::scale_color_manual(
      name = "name",
      values = dag_colors()
    ) +
    ggplot2::scale_fill_manual(
      name = "name",
      values = dag_colors()
    ) +
    ggplot2::theme(legend.position = "none")
}


#' Title
#'
#' @param data
#'
#' @returns
#' @export
#'
#' @examples
#' targets::tar_read(lst_olr_interact_merged) |> merged_interact_table()
merged_interact_table <- function(data) {
  data |>
    purrr::imap(\(.x, .i) {
      purrr::map(.x, \(.y) {
        gtsummary::tbl_regression(.y, exponentiate = TRUE) |>
          compress_logical_vars()
      }) |>
        tbl_merge_named() |>
        gtsummary::remove_row_type(
          variables = -dplyr::one_of(.i),
          type = "all"
        )
    }) |>
    gtsummary::tbl_stack() |>
    fix_labels()
}


default_baseline <- function(data, by, ...) {
  out <- data |>
    gtsummary::tbl_summary(
      by = by,
      type = list(svd_score_04 = "continuous"),
      label = list(
        simple_microbleed = "Microbleeds (≥1)",
        simple_lacunes = "Lacunes (≥1)",
        simple_wmh = "Begininng-large confluenting areas of WMH (Fazekas 2-3)",
        simple_atrophy = "Moderate-severe global atrophy (GCA 2-3)"
      ),
      digits = list(age ~ 0, svd_score_04 ~ 0)
    )

  if (!is.null(by)) {
    out <- out |>
      gtsummary::add_p() |>
      gtsummary::add_overall()
  }

  out |>
    bstfun::add_variable_grouping(
      "SVD score features" = c("simple_microbleed", "simple_lacunes", "simple_wmh", "simple_atrophy")
    )
}

default_regression_table <- function(data, remove_all_but=TRUE) {
  out <- data |>
    gtsummary::tbl_regression(
      show_single_row = dplyr::where(is.logical), exponentiate = TRUE, p.values = FALSE
    ) |>
    fix_labels()
  
  if (remove_all_but){
    out <- out |> gtsummary::remove_row_type(variables = -pase_0_q, type = "all")
  }
  out
}

## https://www.danieldsjoberg.com/gtsummary/articles/themes.html?q=theme#writing-themes
local_gtsummary_theme <-
  list(
    # round large p-values to two places
    `pkgwide-fn:pvalue_fun` = gtsummary::label_style_pvalue(digits = 2),
    `pkgwide-fn:prependpvalue_fun` = gtsummary::label_style_pvalue(digits = 2, prepend_p = TRUE),
    # report median (Q1 - Q2) and n (percent) as default stats in `tbl_summary()`
    `tbl_summary-arg:statistic` = list(
      gtsummary::all_continuous() ~ "{median} ({p25}-{p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    `tbl_summary-arg:missing` = "ifany",
    `tbl_summary-arg:missing_text` = "Missing",
    `add_difference-fn:addnl-fn-to-run` = function(x) {
      tryCatch(
        {
          new_header_text <- paste0(
            dplyr::pull(dplyr::filter(
              x$table_styling$header,
              .data$column == "estimate"
            ), "label"), " **(**",
            dplyr::pull(dplyr::filter(
              x$table_styling$header,
              .data$column == "conf.low"
            ), "label"),
            "**)**"
          )
          estimate_footnote <- paste(c(
            dplyr::pull(dplyr::filter(dplyr::filter(
              x$table_styling$footnote_abbrev,
              .data$column %in% "estimate"
            ), dplyr::row_number() ==
              dplyr::n(), !is.na(.data$footnote)), "footnote"),
            "CI = Confidence Interval"
          ), collapse = ", ")
          modify_footnote(
            modify_header(
              x %>% modify_column_merge(
                rows = !!expr(.data$variable %in%
                  !!x$table_body$variable & !is.na(.data$estimate)),
                pattern = "{estimate} ({conf.low} to {conf.high})"
              ),
              estimate = new_header_text
            ),
            estimate = estimate_footnote,
            abbreviation = TRUE
          )
        },
        error = function(e) x
      )
    }, `tbl_regression-fn:addnl-fn-to-run` = function(x) {
      tryCatch(
        {
          new_header_text <- paste0(
            dplyr::pull(dplyr::filter(
              x$table_styling$header,
              .data$column == "estimate"
            ), "label"), " **(",
            style_number(x$inputs$conf.level, scale = 100),
            "% CI)**"
          )
          estimate_footnote <- paste(c(
            dplyr::pull(dplyr::filter(dplyr::filter(
              x$table_styling$footnote_abbrev,
              .data$column %in% "estimate"
            ), dplyr::row_number() ==
              dplyr::n(), !is.na(.data$footnote)), "footnote"),
            "CI = Confidence Interval"
          ), collapse = ", ")
          modify_footnote(
            modify_header(
              x %>% modify_column_merge(rows = !!expr(.data$variable %in%
                !!x$table_body$variable & !is.na(.data$estimate) &
                !.data$reference_row %in% TRUE), pattern = "{estimate} ({conf.low} to {conf.high})"),
              estimate = new_header_text
            ),
            estimate = estimate_footnote,
            abbreviation = TRUE
          )
        },
        error = function(e) x
      )
    }
  )

# gtsummary::check_gtsummary_theme(local_gtsummary_theme)

gtsummary::set_gtsummary_theme(local_gtsummary_theme)



strat_cut <- function(data, strat = NULL, ...) {
  data <- tibble::tibble(data, id = dplyr::row_number())

  if (is.null(strat)) {
    ls <- list(data)
  } else {
    ls <- split(data, strat)
  }

  out <- ls |>
    lapply(\(.x){
      .x$data_cut <- project.aid::quantile_cut(x = .x$data, ...)
      .x
    }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(id)

  out$data_cut
}


gg_export <- function(
    plot,
    filename,
    height = 50,
    ...) {
  ggplot2::ggsave(
    filename = here::here(filename),
    plot = plot,
    width = 140,
    height = height,
    units = "mm",
    dpi = 600,
    scale = 2.4, ...
  )
}

#' Create regression model based on new outcome variable
#'
#' @param data
#' @param vars
#'
#' @returns
#' @export
#'
#' @examples
#' ds <- targets::tar_read(df_formatted)
#' create_n_models(ds, paste0("simple_", c("atrophy", "lacunes", "microbleed", "wmh")), covars = c("pase_0_q", "age", "female_sex"))

create_n_models <- function(data, vars, covars, prefix="simple_") {
  score <- rowSums(data[vars])
  type <- freesearcheR::possible_functions(score)[[1]]
  
  ds <- data.frame(score=score,data[covars])
  
  ls <- do.call(
    freesearcheR::regression_model_list,
      list(data = ds,
           outcome.str = "score",
           fun.descr = type)
  )
  
  m_name <- gsub(prefix,"",vars) |> paste(collapse = "_")
  
  list(ls$model) |> setNames(m_name)
}

# pak::pak("agdamsbo/project.aid")
