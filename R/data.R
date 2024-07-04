#' Velmeshev et al. DEGs
#'
#' From https://www.science.org/doi/10.1126/science.aav8130,
#' data supplement S4 (i.e., \code{aav8130_data-s4.xls})
#'
#' @format ## `velmeshev_gene_df`
#' A data frame with 692 rows (genes) and 14 columns:
#' \describe{
#'   \item{Cell type}{Cell type that the DEG was in}
#'   \item{gene ID}{ENSG gene ID}
#'   \item{Gene name}{Gene name}
#'   \item{Gene biotype}{Biotype}
#'   \item{Fold change}{See original paper for details}
#'   \item{Sample fold change}{See original paper for details}
#'   \item{q value}{See original paper for details}
#'   \item{correlation (bulk mRNA/bulkized nuclear RNA)}{See original paper for details}
#'   \item{Epilepsy DEG}{boolean}
#'   \item{gene group}{boolean}
#'   \item{SFARI gene}{boolean}
#'   \item{Satterstrom}{boolean}
#'   \item{Sanders}{boolean}
#'   \item{cell type-specific expression}{boolean}
#' }
#' @source <https://www.science.org/doi/10.1126/science.aav8130>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' tib <- readxl::read_xls("aav8130_data-s4.xls",
#' sheet = "ASD_DEGs")
#' velmeshev_gene_df <- as.data.frame(tib)
#' usethis::use_data(velmeshev_gene_df)
#' }
"velmeshev_gene_df"

#' Hounkpe et al. housekeeping genes
#'
#' From https://academic.oup.com/nar/article/49/D1/D947/5871367#supplementary-data,
#' gkaa609_Supplemental_Files, specifically, \code{Supplementary_Table1.xlsx})
#'
#' @format ## `housekeeping_df`
#' A data frame with 1003 rows (genes) and 1 column:
#' \describe{
#'   \item{Gene}{Gene name}
#' }
#' @source <https://academic.oup.com/nar/article/49/D1/D947/5871367>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' df <- openxlsx::read.xlsx("Supplementary_Table1.xlsx",
#' sheet = "Gene Model")
#' colnames(df) <- df[1,1]
#' housekeeping_df <- df[-1,,drop = FALSE]
#' usethis::use_data(housekeeping_df)
#' }
"housekeeping_df"

#' SFARI genes
#'
#' From https://gene.sfari.org/database/human-gene/, after clicking on "Download this dataset",
#' accessed on June 15, 2024
#'
#' @format ## `sfari_df`
#' A data frame with 1176 rows (genes) and 10 columns:
#' \describe{
#'   \item{status}{See website for details}
#'   \item{gene.symbol}{Gene symbol}
#'   \item{gene.name}{Full gene name}
#'   \item{ensembl.id}{Gene ID}
#'   \item{chromosome}{Chromosome the gene is on}
#'   \item{genetic.category}{See website for details}
#'   \item{gene.score}{See website for details}
#'   \item{syndromic}{See website for details}
#'   \item{eagle}{See website for details}
#'   \item{number.of.reports}{See website for details}
#' }
#' @source <https://www.science.org/doi/10.1126/science.aav8130>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' sfari_df <- read.csv("SFARI-Gene_genes_03-28-2024release_06-16-2024export.csv")
#' usethis::use_data(sfari_df)
#' }
"sfari_df"

#' Gandal et al. results
#'
#' From https://www.nature.com/articles/s41586-022-05377-7,
#' Supplementary Data 3 (i.e., \code{41586_2022_5377_MOESM5_ESM.xlsx})
#'
#' @format ## `gandal_df`
#' A data frame with 24836 rows (genes) and 3 columns:
#' \describe{
#'   \item{external_gene_name}{Gene symbol}
#'   \item{WholeCortex_ASD_logFC}{Estimated log fold change}
#'   \item{WholeCortex_ASD_FDR}{FDR q-value}
#' }
#' @source <https://www.nature.com/articles/s41586-022-05377-7>
#' @examples
#' \dontrun{
#' # How the data was loaded
#' df <- openxlsx::read.xlsx("/Users/kevinlin/Downloads/41586_2022_5377_MOESM5_ESM.xlsx",
#'                           sheet = "DEGene_Statistics")
#' gandal_df <- df[,c("external_gene_name", "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR")]
#' usethis::use_data(gandal_df)
#' }
"gandal_df"


