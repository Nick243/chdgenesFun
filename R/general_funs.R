read_chdgenes <- function(my_path){
  df <- vroom::vroom(my_path,
                     delim = ";", guess_max = 10000) |>
    dplyr::select(1:57) |>
    dplyr::rename_all(tolower) |>
    dplyr::mutate(collection_date = lubridate::mdy(collection_date),
                  date_shipped = lubridate::mdy(date_shipped),
                  create_date = lubridate::mdy(create_date),
                  disposal_date = lubridate::mdy(disposal_date),
                  dna_qcdate = lubridate::mdy(dna_qcdate),
                  date_ordered = lubridate::mdy(date_ordered))
  return(df)
}

remove_relatives <- function(my_df){
  df <- my_df |>
    tidyr::separate(blind_id, into = c("drop", "family_id", "member"), sep = "-", remove = F) |>
    dplyr::mutate(member = ifelse(is.na(member), "00", member)) |>
    dplyr::filter(member %in% c("00", "01", "02"))
  return(df)
}

remove_probands <- function(my_df){
  df <- my_df |>
    dplyr::filter(nchar(blind_id) != 7)
  return(df)
}

remove_decommissioned <- function(my_df){
  df <- my_df |>
    dplyr::filter(nchar(blind_id) != 14) |>
    dplyr::filter(nchar(blind_id) != 16)
  return(df)
}




