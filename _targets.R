# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("tibble") # Packages that your targets need for their tasks.
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(
    name = ls_data,
    command = full_dataset(key = "SVD_REDCAP_API"),
    # cue = targets::tar_cue(mode = "always"),
    description = "Main dataset acquisition"
    # format = "qs" # Efficient storage for general data objects.
  ),
  tar_target(
    name = flowchart_dot,
    command = create_flowchart(ls_data, export.path = here::here("images/flow.dot")),
    description = "Flowchart dot file export for better plotting"
    # format = "qs" # Efficient storage for general data objects.
  ),
  tar_target(
    name = flowchart_clin,
    command = clinical_flowchart(df_complete, export.path = here::here("images/flow_main.dot")),
    description = "Flowchart dot file export for better plotting"
    # format = "qs" # Efficient storage for general data objects.
  ),
  tar_target(
    name = complete_scores,
    command = final_score(ls_data),
    description = "Tibble of merged SVD scores with consensus assessments"
  ),
  tar_target(
    name = clin_data,
    command = clin_data_all(),
    description = "Clinical data acquisition"
  ),
  tar_target(
    name = df_complete,
    command = dplyr::left_join(clin_data, 
                               cut_scores(data = complete_scores, prefix = "")) |>
      dplyr::mutate(include = !is.na(pase_0) & !is.na(simple_score)),
    description = "Joined data set of SVD scores and clinical data"
  ),
  tar_target(
    name = df_formatted,
    command = formatted_df(df_complete),
    description = "Formatted complete dataset"
  ),
  tar_target(
    name = list_multi_var_olr,
    command = mulitvar_olr(df_formatted),
    description = "Merged uni and mulitvariable olr analyses"
  ),
  tar_target(
    name = list_multi_var_olr_steps,
    command = mulitvar_olr(df_formatted,reg.fun = regression_uni_mini_multi_olr),
    description = "Merged univariable, minimally adjusted and mulitvariable olr analyses"
  ),
  tar_target(
    name = list_multi_var_olr_interact,
    command = mulitvar_olr_interact(df_formatted),
    description = "List of minimally adjust interaction analyses"
  ),
  tar_target(
    name = lst_multi_olr,
    command = planned_multi_olr(df_formatted),
    description = "Combined list of planned OLR models and formatted tables"
  ),
  tar_target(
    name = lst_multi_olr_orig,
    command = planned_multi_olr(df_formatted,gen_exp = c("age", "female_sex")),
    description = "Combined list of originally planned OLR models and formatted tables"
  ),
  tar_target(
    name = lst_multi_olr_plot,
    command = multi_coef_plot(lst_multi_olr),
    description = "Plots planned OLR model coefs"
  ),
  tar_target(
    name = rds_collinearity_test_orig,
    command = write_collinearity_test_rds(data = lst_multi_olr_orig,file = here::here("check_orig.rds")),
    description = "Save orig collinearity RDS"
  ),
  tar_target(
    name = rds_collinearity_test_mod,
    command = write_collinearity_test_rds(data = lst_multi_olr,file = here::here("check_reduced.rds")),
    description = "Save modified collinearity RDS"
  ),
  tar_target(
    name = lst_olr_interact_merged,
    command = multi_string_olr(data = df_formatted, out = "simple_score_f", pase = "pase_0_q", gen_exp = c("age", "female_sex"), main_exp = purrr::pluck(model_input(), "main_exp")),
    description = "Merged list of all (extended) planned anayses"
  # ),
  # tar_target(
  #   name = tbl_olr_interact_merged,
  #   command = merged_interact_table(lst_olr_interact_merged),
  #   description = "Merged table of all (extended) planned anayses"
  )#,
  # tar_target(
  #   name = gt_mini_multi_merge,
  #   command = mini_multi_merge(lst_olr_interact_merged),
  #   description = "Merged table of all (extended) planned anayses"
  # )
)
