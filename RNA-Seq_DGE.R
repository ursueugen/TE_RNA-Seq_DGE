

get_tx2gene <- function() {
  df = read.table(
    TX2GENE_PATH,
    header = TRUE,
    sep = '\t',
    stringsAsFactors = FALSE
  )
  df_copy = data.frame(df)
  df[, 1] = df_copy[, 2]
  df[, 2] = df_copy[, 1]
  names(df) = c(names(df_copy)[2], names(df_copy)[1])
  return(df)
}

get_comparisons_table <- function() {
  return(
    read.table(
      "TE_Datasets_Aging_Comparisons.tsv",
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
  )
}

get_salmon_filepath <- function(gse, gsm) {
  
  
  if (gse == "GSE106669") {
    gse_adjusted = "GSE106669_GSE106670"
  } else {
    gse_adjusted = gse
  }
  
  dirs = list.dirs(path = file.path(SAMPLES_BASE_PATH, gse_adjusted, gsm),
                   recursive = FALSE)
  if (length(dirs) > 1) {
    stop("Not accounted for: > 1 run per sample.")
  }
  dir_split = strsplit(dirs[1], "/")[[1]]
  run = dir_split[length(dir_split)]
  
  return(file.path(
    SAMPLES_BASE_PATH,
    gse_adjusted,
    gsm,
    run,
    paste0("quant_", gse_adjusted, "_", gsm, "_", run),
    "quant.sf"
  ))
}

SAMPLES_BASE_PATH = "/home/rstudio/samples"
TX2GENE_PATH = "/home/rstudio/ensembl/101/species/Homo_sapiens/tx2gene_human.tsv"

# MAIN
GSE = "GSE106669"  # GSE39170 GSE106669

comparisons = get_comparisons_table()
comparisons = comparisons[comparisons$GSE == GSE,]

for (i in 1:nrow(comparisons)) {

  gsms_favorable = strsplit(comparisons[i, "GSMs_1"], ",")[[1]]
  gsms_unfavorable = strsplit(comparisons[i, "GSMs_2"], ",")[[1]]
  title_favorable = comparisons[i, "Title_1"]
  title_unfavorable = comparisons[i, "Title_2"]
  
  samples = data.table::data.table(
    GSM = c(gsms_favorable, gsms_unfavorable),
    condition = c(rep("favorable", length(gsms_favorable)), rep(
      "unfavorable", length(gsms_unfavorable)
    ))
  )
  # Load the data
  files = sapply(samples$GSM, function(gsm) {
    get_salmon_filepath(gse = GSE, gsm = gsm)
  })
  tx2gene = get_tx2gene()
  txi = tximport::tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE
  )
  
  # DGE
  ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi,
                                             colData = samples,
                                             design = ~ condition)
  
  dds <- DESeq2::DESeq(ddsTxi)
  #DESeq2::resultsNames(dds)
  res = DESeq2::results(dds, contrast = c("condition", "favorable", "unfavorable"))
  # res = DESeq2::lfcShrink(dds, coef="favorable_vs_unfavorable", type="apeglm")
  
  res = data.frame(res)
  write.table(
    res,
    file = paste0(
      GSE,
      "_favorable-",
      title_favorable,
      "_vs_unfavorable-",
      title_unfavorable,
      ".tsv"
    ),
    sep = "\t"
  )
}
