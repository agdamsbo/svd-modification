
## Get completed scores
# df <- targets::tar_read(complete_scores)

df <- full_dataset(key = "SVD_REDCAP_API",raw_or_label = "raw") |> 
  final_score()

## Get instrument metadata to reuse
meta <- REDCapR::redcap_metadata_read(redcap_uri = "https://redcap.au.dk/api/",
                              token = keyring::key_get("SVD_REDCAP_API"))

## Renaming to avoid confusion later on
names(df) <- gsub("svd_","",names(df))

## Selecting fields to use for new instrument
## and removing all "svd_" references
final_meta <- meta$data |> 
  dplyr::filter(gsub("svd_","",field_name) %in% gsub("*___[0-9]$","",names(df))) |> 
  dplyr::mutate(dplyr::across(dplyr::everything(),\(.x) gsub("svd_","final_",.x)),
                branching_logic="")

merged_meta <- df |> dplyr::transmute(final_consensus=factor(consensus)) |> 
  REDCapCAST::ds2dd_detailed() |> 
  purrr::pluck("meta") |> 
  (\(.x){
    dplyr::bind_rows(
      final_meta,
      .x
    )
  })() |> 
  dplyr::mutate(form_name="final_score")

## Export new instrument
merged_meta[-c(2:3),] |> REDCapCAST::create_instrument_meta(record.id = TRUE)

## Renaming data
names(df)[-1] <- paste0("final_",names(df)[-1])

df |> dplyr::select(-final_user,-final_missing_reason) |> 
  dplyr::mutate(final_consensus=as.numeric(factor(final_consensus))) |> 
  REDCapR::redcap_write(redcap_uri = "https://redcap.au.dk/api/",
                        token = keyring::key_get("SVD_REDCAP_API"))