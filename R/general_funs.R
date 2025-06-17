#' read_chdgenes
#'
#' @param my_path The full CHDGENES sample data inventory file as a semicolon delimited text file.
#'
#' @returns The full CHDGENES sample data inventory file as a R data.frame
#' @export
#'
#' @examples
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

#' remove_relatives
#'
#' @param my_df A data.frame containing the a list of blind IDs. Often this will be the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with specimens from only probands, mom (-01), and dad (-02). All relative samples will be removed.
#' @export
#'
#' @examples
remove_relatives <- function(my_df){
  df <- my_df |>
    tidyr::separate(blind_id, into = c("drop", "family_id", "member"), sep = "-", remove = F) |>
    dplyr::mutate(member = ifelse(is.na(member), "00", member)) |>
    dplyr::filter(member %in% c("00", "01", "02"))
  return(df)
}

#' remove_probands
#'
#' @param my_df A data.frame containing the a list of blind IDs. Often this will be the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with all proband data removed.
#' @export
#'
#' @examples
remove_probands <- function(my_df){
  df <- my_df |>
    dplyr::filter(nchar(blind_id) != 7)
  return(df)
}

#' remove_decommissioned
#'
#' @param my_df A data.frame containing the a list of blind IDs. Often this will be the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with all decommissioned samples removed.
#' @export
#'
#' @examples
remove_decommissioned <- function(my_df){
  df <- my_df |>
    dplyr::filter(nchar(blind_id) != 14) |>
    dplyr::filter(nchar(blind_id) != 16)
  return(df)
}




