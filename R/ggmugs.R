#' Visualizing Multiple Genetic Studies
#'
#' This function processes and combines summary statistics from multiple genetic studies and
#' creates a visualization for all studies. The genetic loci are colored based on three significance
#' thresholds to facilitate the visualization of highly significant genomic regions.
#'
#' @param study_name A character vector of names for the studies.
#' @param summary_stat A character vector of file paths where each path points to the summary
#'        statistics data file for the corresponding study. Files should be in a tabular format
#'        readable by `fread` from the `data.table` package.
#'        The files should contain 3 fields: `chr` (Chromosome), `pos` (chromosome position), and `p` (association p-value).
#'        The positions of multiple GWAS summary statistics should have consistent genome builds.
#' @param p1 The first significance level threshold for p-values (default is 1e-3).
#' @param p2 The second, more stringent significance level threshold for p-values (default is 5e-5).
#' @param p3 The most stringent significance level threshold for p-values (default is 1e-8).
#' @param color1 The color for points below the first significance level (default is "#FFFFE0").
#' @param color2 The color for points between the first and second significance levels (default is "#FFC300").
#' @param color3 The color for points above the second significance level (default is "#FF5733").
#'
#' @return A `ggplot` object representing the visualization with the specified data.
#' @importFrom ggplot2 ggplot aes geom_point ylab xlab scale_x_continuous theme_classic theme
#' @importFrom dplyr arrange distinct select filter full_join group_by summarize left_join mutate
#' @importFrom data.table fread
#' @importFrom purrr map
#' @importFrom tidyr gather
#' @examples
#' ### NOT RUN
#' # ggmugs(
#' #   study_name = c("study1", "study2", "study3", "study4", "study5"),
#' #   summary_stat = c("https://raw.githubusercontent.com/Broccolito/ggmugs_data/main/sumstat1.txt",
#' #                    "https://raw.githubusercontent.com/Broccolito/ggmugs_data/main/sumstat2.txt",
#' #                    "https://raw.githubusercontent.com/Broccolito/ggmugs_data/main/sumstat3.txt",
#' #                    "https://raw.githubusercontent.com/Broccolito/ggmugs_data/main/sumstat4.txt",
#' #                    "https://raw.githubusercontent.com/Broccolito/ggmugs_data/main/sumstat5.txt"),
#' #   p1 = 1e-4,
#' #   p2 = 1e-6,
#' #   p3 = 1e-8,
#' #   color1 = "#FFFFE0",
#' #   color2 = "#FFC300",
#' #   color3 = "#FF5733"
#' # )
#' @export

ggmugs = function(
    study_name = c("sumstat1", "sumstat2", "sumstat3", "sumstat4", "sumstat5"),
    summary_stat = c("data/sumstat1.txt", "data/sumstat2.txt", "data/sumstat3.txt",
                     "data/sumstat4.txt", "data/sumstat5.txt"),
    p1 = 1e-3, p2 = 5e-5, p3 = 1e-8,
    color1 = "#FFFFE0", color2 = "#FFC300", color3 = "#FF5733"
){

  # Create a data frame to map study names to summary statistic files
  input_files = data.frame(
    study_name = study_name,
    summary_stat = summary_stat
  )
  # Split the data by study names
  input_files = split(input_files, input_files[["study_name"]])

  # Load and process data files for each study
  data = purrr::map(input_files, function(x){
    message(paste0("Loading ", x[["study_name"]], "...\n"))
    d = data.table::fread(x[["summary_stat"]])
    d[["study_name"]] = x[["study_name"]]
    d[["chr_pos"]] = paste0(d[["chr"]], "_", d[["pos"]])
    d[["logp"]] = -log(d[["p"]],10)
    d = dplyr::arrange(d, chr, pos, dplyr::desc(logp))
    d = dplyr::distinct(d, chr_pos, .keep_all = TRUE)
    d = dplyr::select(d, chr, pos, chr_pos, logp, study_name)
    return(d)
  })

  # Combine data from multiple studies
  data1 = data[[1]]
  message(paste0("Combining ", data1[["study_name"]][1], "...\n"))
  names(data1)[names(data1) == "logp"] = data1[["study_name"]][1]
  data1 = dplyr::select(data1, -study_name)
  for(i in 2:length(data)){
    data2 = data[[i]]
    message(paste0("Combining ", data2[["study_name"]][1], "...\n"))
    data_combined = dplyr::full_join(data1, data2, by = "chr_pos", suffix = c("", "_combined"))
    names(data_combined)[names(data_combined) == "logp"] = data2[["study_name"]][1]
    data_combined = dplyr::select(data_combined, -c(chr_combined, pos_combined, study_name))
    data1 = data_combined
  }

  # Final data selection and cleanup
  data = dplyr::select(data_combined, chr, pos, chr_pos, dplyr::everything())
  rm(data_combined, data1, data2)

  # Calculate chromosome lengths and positions
  data_chromosome_position = data |>
    dplyr::group_by(chr) |>
    dplyr::summarise(chromosome_length = max(pos)) |>
    dplyr::mutate(tot = cumsum(as.numeric(chromosome_length)) - chromosome_length) |>
    dplyr::select(-chromosome_length)

  # Merge and arrange data with cumulative base pair positions
  data = dplyr::left_join(data, data_chromosome_position, by = "chr") |>
    dplyr::arrange(chr, pos) |>
    dplyr::mutate(bp_cum = pos + tot) |>
    dplyr::select(chr, pos, chr_pos, tot, bp_cum, dplyr::everything())

  # Prepare the plot data and axis labels
  axisdf = data |>
    dplyr::group_by(chr) |>
    dplyr::summarize(center = (max(bp_cum) + min(bp_cum))/2)

  # Filter data for plotting based on significance thresholds
  data = data |>
    tibble::as_tibble() |>
    tidyr::gather(key = "study_name", value = "logp", -c(chr, pos, chr_pos, tot, bp_cum)) |>
    dplyr::filter(logp >= -log(p1, base = 10))

  data_p1 = dplyr::filter(data, logp < -log(p2, base = 10))
  data_p2 = dplyr::filter(data, logp >= -log(p2, base = 10), logp < -log(p3, base = 10))
  data_p3 = dplyr::filter(data, logp >= -log(p3, base = 10))

  # Create the ggplot object with custom styling and color based on significance
  message("Visualizing summary statistics...\n")
  plt = ggplot2::ggplot(data = data, ggplot2::aes(x = bp_cum, y = study_name)) +
    ggplot2::geom_point(data = data_p1, color = color1, shape = 15) +
    ggplot2::geom_point(data = data_p2, color = color2, shape = 15) +
    ggplot2::geom_point(data = data_p3, color = color3, shape = 15) +
    ggplot2::ylab("") +
    ggplot2::xlab("Chromosome") +
    ggplot2::scale_x_continuous(label = axisdf$chr, breaks = axisdf$center, expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position="none",
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank()); plt

  # Shade alternate chromosomes for visual distinction
  chrs = unique(data[["chr"]])
  shaded_chrs = chrs[-seq(1, length(chrs), 2)]
  for(shaded_chr in shaded_chrs){
    plt = plt + ggplot2::annotate("rect", xmin = min(dplyr::filter(data, chr == shaded_chr)[["bp_cum"]]),
                                  xmax = max(dplyr::filter(data, chr == shaded_chr)[["bp_cum"]]),
                                  ymin = 0, ymax = Inf, alpha = 0.3, fill="gray")
  }

  return(plt)
}
