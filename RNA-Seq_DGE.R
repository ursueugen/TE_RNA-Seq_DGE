
get_tx2gene <- function() {
  df = read.table(TX2GENE_PATH, header=TRUE, sep='\t', stringsAsFactors = FALSE)
  df_copy = copy(df)
  df[, 1] = df_copy[, 2]
  df[, 2] = df_copy[, 1]
  names(df) = c(names(df_copy)[2], names(df_copy)[1])
  return(df)
}

get_comparisons_table <- function() {
  return(data.table::data.table(
    read.table(
      "TE_Datasets_Aging_Comparisons.tsv",
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  ))
}

samples_from_GSE <- function(gse) {}

get_salmon_filepath <- function(gse, gsm) {
  dirs = list.dirs(path = file.path(SAMPLES_BASE_PATH, gse, gsm), recursive = FALSE)
  if ( length(dirs) > 1 ) {stop("Not accounted for: > 1 run per sample.")}
  dir_split = strsplit(dirs[1], "/")[[1]]
  run = dir_split[length(dir_split)]
  
  return(
    file.path(SAMPLES_BASE_PATH, gse, gsm, run, paste0("quant_", gse, "_", gsm, "_", run), "quant.sf")
  )
}


SAMPLES_BASE_PATH = "/home/rstudio/samples"
TX2GENE_PATH = "/home/rstudio/ensembl/101/species/Homo_sapiens/tx2gene_human.tsv"
GSE = "GSE39170"
GSM = "GSM957471"
# Get samples' table
comparisons = get_comparisons_table()
gsms_favorable = strsplit(comparisons[GSE == "GSE39170", GSMs_1], ",")[[1]]
gsms_unfavorable = strsplit(comparisons[GSE == "GSE39170", GSMs_2], ",")[[1]]
samples = data.table::data.table(
  GSM = c(gsms_favorable, gsms_unfavorable),
  condition = c(rep("favorable", length(gsms_favorable)), rep(
    "unfavorable", length(gsms_unfavorable)
  ))
)

# Load the data
files = sapply(samples$GSM, function(gsm) {get_salmon_filepath(gse = GSE, gsm = gsm)} )
tx2gene = get_tx2gene()
txi = tximport::tximport(files, type="salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# DGE
ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)

dds <- DESeq2::DESeq(ddsTxi)
DESeq2::resultsNames(dds)
res = DESeq2::results(dds, contrast = c("condition", "faborable", "unfavorable"))

# res = DESeq2::lfcShrink(dds, coef="favorable_vs_unfavorable", type="apeglm")
