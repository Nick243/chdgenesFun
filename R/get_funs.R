
#' get_shipped_wgs
#'
#' Generates a data.frame of blinded IDs shipped for short read whole genome sequencing
#'
#' @param my_df A data.frame containing the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_shipped_wgs <- function(my_df){
  df <- my_df |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose == "WGS") |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_wgs = "Y") |>
    dplyr::rename("wgs_ship_date" = "date_shipped",
                  "wgs_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_wgs, wgs_ship_date, wgs_order_id)
  return(df)
}

#' get_shipped_wes
#'
#' Generates a data.frame of blinded IDs shipped for exome sequencing
#'
#' @param my_df A data.frame containing the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_shipped_wes <- function(my_df){
  df <- my_df |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose == "WES") |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_wes = "Y") |>
    dplyr::rename("wes_ship_date" = "date_shipped",
                  "wes_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_wes, wes_ship_date, wes_order_id)
  return(df)
}

#' get_shipped_array
#'
#'  Generates a data.frame of blinded IDs shipped for genotying array
#'
#' @param my_df A data.frame containing the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_shipped_array<- function(my_df){
  df <- my_df |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose == "Genotyping") |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_array = "Y") |>
    dplyr::rename("array_ship_date" = "date_shipped",
                  "array_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_array, array_ship_date, array_order_id)
  return(df)
}

#' get_shipped_mips
#'
#'  Generates a data.frame of blinded IDs shipped for MIPs
#'
#' @param my_df A data.frame containing the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_shipped_mips<- function(my_df){
  df <- my_df |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose %in% c("MIPs", "MIPs targeted sequencing")) |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_mips = "Y") |>
    dplyr::rename("mips_ship_date" = "date_shipped",
                  "mips_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_mips, mips_ship_date, mips_order_id)
  return(df)
}

#' get_shipped_lrwgs
#'
#' Generates a data.frame of blinded IDs shipped for long read whole genome sequencing (PacBio or Oxford Nanopore)
#'
#' @param my_df A data.frame containing the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_shipped_lrwgs<- function(my_df){
  df <- my_df |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose %in% c("Long Read", "Oxford Nanopore", "PacBio  WGS")) |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_lrwgs = "Y") |>
    dplyr::rename("lrwgs_ship_date" = "date_shipped",
                  "lrwgs_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_lrwgs, lrwgs_ship_date, lrwgs_order_id)
  return(df)
}

#' get_total_avail_dna
#'
#' Generates a data.frame of blinded IDs and the total available (in circulation or reserve) QC pass DNA in ug for each.
#'
#' @param my_df A data.frame containing the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_total_avail_dna <- function(my_df){
  df <- my_df |>
    dplyr::filter(dna_qcstatus == "QC completed") |>
    dplyr::filter(status %in% c("In Circulation", "Reserved")) |>
    dplyr::group_by(blind_id) |>
    dplyr::summarise(total_avail_qcpass_dna = sum(nanodrop_mass_ug, na.rm = T)) |>
    dplyr::ungroup()
  return(df)
}

#' get_source_type
#'
#'  Generates a data.frame of blinded IDs and the source type of the submitted specimens. If a blinded ID has multiple submission types they will be concatenated (e.g., wb|sal).
#'
#' @param my_df A data.frame containing the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_source_type <- function(my_df){
  df <- my_df |>
    dplyr::distinct(blind_id, source_type) |>
    dplyr::group_by(blind_id) |>
    dplyr::summarise(dna_source = paste(source_type, collapse = "|"))
  return(df)
}

#' get_probands_only
#'
#' Removes non-proband blinded IDs from a data.frame
#'
#' @param my_df A data.frame containing the a list of blind IDs. Often this will be the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_probands_only <- function(my_df){
  df <- my_df |>
    dplyr::filter(nchar(blind_id) == 7)
  return(df)
}

#' get_in_trio
#'
#' Generates a data.frame of family IDs and a column denoting whether the ID is part of a complete trio (proband, mom, dad)
#'
#' @param my_df  A data.frame containing the a list of blind IDs. Often this will be the full CHDGENES sample inventory (e.g., aliquot level data).
#'
#' @returns A data.frame with columns as denoted by the function title.
#' @export
#'
#' @examples
get_in_trio <- function(my_df){
  df <- remove_relatives(my_df) |>
    add_family_id() |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::group_by(family_id) |>
    dplyr::count(family_id) |>
    dplyr::rename("in_biobank_trio" = "n") |>
    dplyr::mutate(in_biobank_trio = ifelse(in_biobank_trio == 3, "Y", "N"))
  return(df)
}
