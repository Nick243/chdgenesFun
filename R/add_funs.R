#' add_shipped_wgs
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns denoting if the blinded ID has been shipped for short read whole genome sequencing, the date of shipment, and the external order ID number.
#' @export
#'
#' @examples
add_shipped_wgs <- function(my_df, full_chdgenes){
  wgs_ship_df <- full_chdgenes |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose == "WGS") |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_wgs = "Y") |>
    dplyr::rename("wgs_ship_date" = "date_shipped",
                  "wgs_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_wgs, wgs_ship_date, wgs_order_id)

  working_df <- dplyr::left_join(my_df, wgs_ship_df, by = "blind_id")
  return(working_df)
}

#' add_shipped_wes
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns denoting if the blinded ID has been shipped for exome sequencing, the date of shipment, and the external order ID number.
#' @export
#'
#' @examples
add_shipped_wes <- function(my_df, full_chdgenes){
  wes_ship_df <- full_chdgenes |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose == "WES") |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_wes = "Y") |>
    dplyr::rename("wes_ship_date" = "date_shipped",
                  "wes_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_wes, wes_ship_date, wes_order_id)

  working_df <- dplyr::left_join(my_df, wes_ship_df, by = "blind_id")
  return(working_df)
}

#' add_shipped_array
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns denoting if the blinded ID has been shipped for genotyping array, the date of shipment, and the external order ID number.
#' @export
#'
#' @examples
add_shipped_array <- function(my_df, full_chdgenes){
  array_ship_df <- full_chdgenes |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose == "Genotyping") |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_array = "Y") |>
    dplyr::rename("array_ship_date" = "date_shipped",
                  "array_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_array, array_ship_date, array_order_id)

  working_df <- dplyr::left_join(my_df, array_ship_df, by = "blind_id")
  return(working_df)
}

#' add_shipped_mips
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns denoting if the blinded ID has been shipped for MIPs, the date of shipment, and the external order ID number.
#' @export
#'
#' @examples
add_shipped_mips <- function(my_df, full_chdgenes){
  mips_ship_df <- full_chdgenes |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose %in% c("MIPs", "MIPs targeted sequencing")) |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_mips = "Y") |>
    dplyr::rename("mips_ship_date" = "date_shipped",
                  "mips_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_mips, mips_ship_date, mips_order_id)

  working_df <- dplyr::left_join(my_df, mips_ship_df, by = "blind_id")
  return(working_df)
}

#' add_shipped_lrwgs
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns denoting if the blinded ID has been shipped for lrWGS, the date of shipment, and the external order ID number.
#' @export
#'
#' @examples
add_shipped_lrwgs <- function(my_df, full_chdgenes){
  lr_ship_df <- full_chdgenes |>
    dplyr::arrange(blind_id, date_shipped) |>
    dplyr::filter(shipment_purpose %in% c("Long Read", "Oxford Nanopore", "PacBio  WGS")) |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::mutate(shipped_lrwgs = "Y") |>
    dplyr::rename("lrwgs_ship_date" = "date_shipped",
                  "lrwgs_order_id" = "external_order_id") |>
    dplyr::select(blind_id, shipped_lrwgs, lrwgs_ship_date, lrwgs_order_id)

  working_df <- dplyr::left_join(my_df, lr_ship_df, by = "blind_id")
  return(working_df)
}

#' add_family_id
#'
#' Adds the five digit family ID to an existing data.frame
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#'
#' @returns Additional column with the 5 digit family ID
#' @export
#'
#' @examples
add_family_id <- function(my_df){
  df <- my_df |>
    tidyr::separate(blind_id, into = c("drop", "family_id", "member"), sep = "-", remove = F) |>
    dplyr::select(-drop, -member)
}

#' add_mom_wgs
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_mom_wgs <- function(my_df, full_chdgenes){
  mom_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "01") |>
    get_shipped_wgs() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_wgs, wgs_ship_date, wgs_order_id) |>
    dplyr::rename("mom_shipped_wgs" = "shipped_wgs",
                  "mom_wgs_ship_date" = "wgs_ship_date",
                  "mom_wgs_order_id" = "wgs_order_id")

  working_df <- dplyr::left_join(my_df, mom_df, by = "family_id")
  return(working_df)
}

#' add_mom_wes
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_mom_wes <- function(my_df, full_chdgenes){
  mom_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "01") |>
    get_shipped_wes() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_wes, wes_ship_date, wes_order_id) |>
    dplyr::rename("mom_shipped_wes" = "shipped_wes",
                  "mom_wes_ship_date" = "wes_ship_date",
                  "mom_wes_order_id" = "wes_order_id")

  working_df <- dplyr::left_join(my_df, mom_df, by = "family_id")
  return(working_df)
}

#' add_dad_wgs
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_dad_wgs <- function(my_df, full_chdgenes){
  dad_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "02") |>
    get_shipped_wgs() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_wgs, wgs_ship_date, wgs_order_id) |>
    dplyr::rename("dad_shipped_wgs" = "shipped_wgs",
                  "dad_wgs_ship_date" = "wgs_ship_date",
                  "dad_wgs_order_id" = "wgs_order_id")

  working_df <- dplyr::left_join(my_df, dad_df, by = "family_id")
  return(working_df)
}

#' add_dad_wes
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_dad_wes <- function(my_df, full_chdgenes){
  dad_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "02") |>
    get_shipped_wes() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_wes, wes_ship_date, wes_order_id) |>
    dplyr::rename("dad_shipped_wes" = "shipped_wes",
                  "dad_wes_ship_date" = "wes_ship_date",
                  "dad_wes_order_id" = "wes_order_id")

  working_df <- dplyr::left_join(my_df, dad_df, by = "family_id")
  return(working_df)
}

#' add_mom_array
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_mom_array <- function(my_df, full_chdgenes){
  mom_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "01") |>
    get_shipped_array() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_array, array_ship_date, array_order_id) |>
    dplyr::rename("mom_shipped_array" = "shipped_array",
                  "mom_array_ship_date" = "array_ship_date",
                  "mom_array_order_id" = "array_order_id")

  working_df <- dplyr::left_join(my_df, mom_df, by = "family_id")
  return(working_df)
}

#' add_mom_mips
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_mom_mips <- function(my_df, full_chdgenes){
  mom_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "01") |>
    get_shipped_mips() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_mips, mips_ship_date, mips_order_id) |>
    dplyr::rename("mom_shipped_mips" = "shipped_mips",
                  "mom_mips_ship_date" = "mips_ship_date",
                  "mom_mips_order_id" = "mips_order_id")

  working_df <- dplyr::left_join(my_df, mom_df, by = "family_id")
  return(working_df)
}

#' add_mom_lrwgs
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_mom_lrwgs <- function(my_df, full_chdgenes){
  mom_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "01") |>
    get_shipped_lrwgs() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_lrwgs, lrwgs_ship_date, lrwgs_order_id) |>
    dplyr::rename("mom_shipped_lrwgs" = "shipped_lrwgs",
                  "mom_lrwgs_ship_date" = "lrwgs_ship_date",
                  "mom_lrwgs_order_id" = "lrwgs_order_id")

  working_df <- dplyr::left_join(my_df, mom_df, by = "family_id")
  return(working_df)
}

#' add_dad_array
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_dad_array <- function(my_df, full_chdgenes){
  dad_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "02") |>
    get_shipped_array() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_array, array_ship_date, array_order_id) |>
    dplyr::rename("dad_shipped_array" = "shipped_array",
                  "dad_array_ship_date" = "array_ship_date",
                  "dad_array_order_id" = "array_order_id")

  working_df <- dplyr::left_join(my_df, dad_df, by = "family_id")
  return(working_df)
}

#' add_dad_mips
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_dad_mips <- function(my_df, full_chdgenes){
  dad_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "02") |>
    get_shipped_mips() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_mips, mips_ship_date, mips_order_id) |>
    dplyr::rename("dad_shipped_mips" = "shipped_mips",
                  "dad_mips_ship_date" = "mips_ship_date",
                  "dad_mips_order_id" = "mips_order_id")

  working_df <- dplyr::left_join(my_df, dad_df, by = "family_id")
  return(working_df)
}

#' add_dad_lrwgs
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional columns as denoted by the function title.
#' @export
#'
#' @examples
add_dad_lrwgs <- function(my_df, full_chdgenes){
  dad_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "02") |>
    get_shipped_lrwgs() |>
    add_family_id() |>
    dplyr::select(family_id, shipped_lrwgs, lrwgs_ship_date, lrwgs_order_id) |>
    dplyr::rename("dad_shipped_lrwgs" = "shipped_lrwgs",
                  "dad_lrwgs_ship_date" = "lrwgs_ship_date",
                  "dad_lrwgs_order_id" = "lrwgs_order_id")

  working_df <- dplyr::left_join(my_df, dad_df, by = "family_id")
  return(working_df)
}

#' add_mom_total_dna
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional column as denoted by the function title. Note this is only QC pass DNA that is listed as either in circulation or in reserve.
#' @export
#'
#' @examples
add_mom_total_dna <- function(my_df, full_chdgenes){
  mom_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "01") |>
    get_total_avail_dna() |>
    add_family_id() |>
    dplyr::select(family_id, total_avail_qcpass_dna) |>
    dplyr::rename("mom_total_avail_qcpass_dna" = "total_avail_qcpass_dna")

  working_df <- dplyr::left_join(my_df, mom_df, by = "family_id")
  return(working_df)
}

#' add_dad_total_dna
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional column as denoted by the function title. Note this is only QC pass DNA that is listed as either in circulation or in reserve.
#' @export
#'
#' @examples
add_dad_total_dna <- function(my_df, full_chdgenes){
  dad_df <- full_chdgenes |>
    remove_relatives() |>
    dplyr::filter(member == "02") |>
    get_total_avail_dna() |>
    add_family_id() |>
    dplyr::select(family_id, total_avail_qcpass_dna) |>
    dplyr::rename("dad_total_avail_qcpass_dna" = "total_avail_qcpass_dna")

  working_df <- dplyr::left_join(my_df, dad_df, by = "family_id")
  return(working_df)
}

#' add_total_avail_dna
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional column as denoted by the function title. This will add columns for all blind ID (proband or parent) passed to the function. Note this is only QC pass DNA that is listed as either in circulation or in reserve.
#' @export
#'
#' @examples
add_total_avail_dna <- function(my_df, full_chdgenes){
  dna_df <- full_chdgenes |>
    dplyr::filter(dna_qcstatus == "QC completed") |>
    dplyr::filter(status %in% c("In Circulation", "Reserved")) |>
    dplyr::group_by(blind_id) |>
    dplyr::summarise(total_avail_qcpass_dna = sum(nanodrop_mass_ug, na.rm = T)) |>
    dplyr::ungroup()

  working_df <- dplyr::left_join(my_df, dna_df, by = "blind_id")
  return(working_df)
}

#' add_in_trio
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional column denoting if the blinded ID is contained in a biobanking trio (e.g., at least one specimen was submitted for the proband and both parents).
#' @export
#'
#' @examples
add_in_trio <- function(my_df, full_chdgenes){
  trio_df <- remove_relatives(chd_df) |>
    add_family_id() |>
    dplyr::distinct(blind_id, .keep_all = T) |>
    dplyr::group_by(family_id) |>
    dplyr::count(family_id) |>
    dplyr::rename("in_biobank_trio" = "n") |>
    dplyr::ungroup() |>
    dplyr::mutate(in_biobank_trio = ifelse(in_biobank_trio == 3, "Y", "N")) |>
    dplyr::mutate(a = paste0("1-", family_id),
                  b = paste0("1-", family_id, "-01"),
                  c = paste0("1-", family_id, "-02")) |>
    tidyr::pivot_longer(cols = c("a", "b", "c"), names_to = "var", values_to = "blind_id") |>
    dplyr::select(in_biobank_trio, blind_id)

  working_df <- dplyr::left_join(my_df, trio_df, by = "blind_id")
  return(working_df)
}

#' add_specimens_in_biobank
#'
#' @param my_df A data.frame containing the blinded IDs of interest
#' @param full_chdgenes A data.frame containing the full CHDGENES sample inventory
#'
#' @returns Additional column as denoting if at least one specimen was received by the biobank for the provided blinded IDs.
#' @export
#'
#' @examples
add_specimens_in_biobank <- function(my_df, full_chdgenes){
  check_df <- full_chdgenes |>
    dplyr::distinct(blind_id) |>
    dplyr::mutate(specimens_in_biobank = "Y")

  working_df <- dplyr::left_join(my_df, check_df, by = "blind_id")
  return(working_df)
}



