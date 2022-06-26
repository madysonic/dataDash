# Dashboard Functions

# DGE object from SQL database
getLBankDGE <- function(){
  
  # Subset count data for only good quality Dx samples
  
  sampleAnno <- dbGetQuery(LBank, "SELECT * FROM Sample_Anno WHERE [exclude] != 'Yes'")
  
  # Select counts columns
  countNames <- dbGetQuery(LBank, "PRAGMA table_info(Counts);")
  
  sampleAnno <- sampleAnno[sampleAnno$LB_ID %in% countNames$name,]
  
  LBIDs <- countNames[countNames$name %in% sampleAnno$LB_ID, "name"]
  LBIDs <- gsub("^", "[", LBIDs)
  LBIDs <- gsub("$", "]", LBIDs)
  LBIDs <- paste(c(LBIDs), collapse = ", ")
  
  
  counts <- dbGetQuery(LBank, paste0("SELECT ", LBIDs, ", [gene_id] FROM Counts;"))
  
  # Make the row names the gene ID
  row.names(counts) <- counts$gene_id
  geneIDcol <- which(colnames(counts) == "gene_id")
  counts <- counts[ ,-geneIDcol]
  
  # # Remove rows containing na.values
  counts <- counts[complete.cases(counts), ]
  
  
  genes <- dbGetQuery(LBank, "SELECT * FROM Gene_Anno;")
  
  # Remove duplicates
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
  genes <- genes[!duplicated(genes$external_gene_name),]
  
  # Reorder genes to match the counts data
  gene_order <- data.frame(ensembl_gene_id = rownames(counts))
  genes <- gene_order %>% left_join(genes, by="ensembl_gene_id")
  
  # Drop genes without a gene symbol from both genes and counts datasets
  
  drop.genes <- is.na(genes$external_gene_name)
  genes <- genes[!drop.genes,]
  counts <- counts[!drop.genes,]
  
  # Set the ensembl gene ID as the rownames
  rownames(genes) <- genes$ensembl_gene_id
  
  # Remove genes with 0 count across all samples
  drop.genes <- rowSums(counts==0)==ncol(counts)
  counts <- counts[!drop.genes,]
  genes <- genes[!drop.genes,]
  
  
  # Normalise
  dge <- DGEList(counts = counts,
                 genes = genes,
                 samples = sampleAnno) %>% 
    calcNormFactors(method = "TMM")
  return(dge)
}

# DGE object from local files
getDGE <- function(){
  # Make the row names the gene ID
  row.names(counts) <- counts$gene_id
  geneIDcol <- which(colnames(counts) == "gene_id")
  counts <- counts[ ,-geneIDcol]
  
  # Remove rows containing na.values
  counts <- counts[complete.cases(counts), ]
  
  # Subset count data for only good quality Dx samples
  sampleAnno <- sample_anno %>% filter(sample %in% c("Dx"), 
                                       exclude != "Yes")
  counts <- counts[ ,match(intersect(sampleAnno$LB_ID, colnames(counts)), colnames(counts))]
  sampleAnno <- sampleAnno[which(sampleAnno$LB_ID %in% colnames(counts)),]
  
  ### Turn into a DGEList
  y <- DGEList(counts = counts)
  
  # Add annotation
  y$samples <- sampleAnno
  
  # Remove duplicates
  genes <- genes[!duplicated(genes$ensembl_gene_id),]
  genes <- genes[!duplicated(genes$external_gene_name),]
  
  # Reorder genes to match the counts data
  gene_order <- data.frame(ensembl_gene_id = rownames(counts))
  genes <- gene_order %>% left_join(genes, by="ensembl_gene_id")
  
  # Drop genes without a gene symbol from both genes and counts datasets
  
  drop.genes <- is.na(genes$external_gene_name)
  genes <- genes[!drop.genes,]
  counts <- counts[!drop.genes,]
  
  # Set the ensembl gene ID as the rownames
  rownames(genes) <- genes$ensembl_gene_id
  # Add to DGE list
  y$genes <- genes
  
  # Remove genes with 0 count across all samples
  x <- ncol(counts)
  drop.genes <- rowSums(y$counts==0)==x
  y <- y[!drop.genes,,keep.lib.sizes=FALSE]
  
  # Remove genes lacking entrez gene ID
  drop.genes <- is.na(y$genes$entrezgene)
  y <- y[!drop.genes,,keep.lib.sizes=FALSE]  
  
  # Normalise
  dge <- calcNormFactors(y, method = "TMM")
  return(dge)
}

plotTSNE <- function(){
  dge2 <- dge[,dge$samples$sample == "Dx", keep.lib.sizes = FALSE]
  # Get normalised counts without batch effects
  logCPM <- cpm(dge, log = TRUE) 
  batch <- as.factor(c(dge$samples$ref))
  batch2 <- as.factor(c(dge$samples$sex))
  logcounts <- removeBatchEffect(logCPM, batch=batch, batch2 = batch2)
  
  # Select most highly variable genes (using logcounts)
  var_genes <- apply(logcounts, 1, var)
  top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:2000]  
  
  # Subset logcounts matrix with the selected genes
  top_var_lcpm <- logcounts[top_var_genes, ]
  
  # Transpose for plotting
  t_highly_variable_lcpm <- data.frame(t(top_var_lcpm))
  
  # Generate dataframe for plotting
  tsne_anno <- as.data.frame(factor(rownames(t_highly_variable_lcpm)))
  
  colnames(tsne_anno) <- "Sample"
  
  # Ensure data remains in same order as columns in logcounts
  tsne_anno$sex <- dge$samples$sex[match(tsne_anno$Sample, dge$samples$LB_ID)]
  tsne_anno$key_alt <- dge$samples$key_alt[match(tsne_anno$Sample, dge$samples$LB_ID)]
  tsne_anno$col <- dge$samples$col[match(tsne_anno$Sample, dge$samples$LB_ID)]
  tsne_anno$patientID <- dge$samples$patientID[match(tsne_anno$Sample, dge$samples$LB_ID)]
  
 
  # Run tSNE 
  set.seed(3) # for reproducibility with MTB tSNE
  tsne_out <- Rtsne(t_highly_variable_lcpm, 
                    pca = FALSE, 
                    perplexity = 30, 
                    theta = 0.0, 
                    check_duplicates = FALSE)
  
  # Plot
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], 
                          subtype = tsne_anno$col, 
                          patientID = tsne_anno$patientID,
                          LB_ID = tsne_anno$Sample,
                          key_alt = tsne_anno$key_alt)
  
  
  
return(tsne_plot)
  
}

plotGenes <- function(genes_to_plot){
  dge <- dge
  logCPM <- cpm(dge, log = TRUE) 
  batch <- as.factor(c(dge$samples$ref))
  batch2 <- as.factor(c(dge$samples$sex))
  logcounts <- removeBatchEffect(logCPM, batch=batch, batch2 = batch2)
  
  
  # Extract count data
  expression <- logcounts %>%
    as.data.frame %>%
    rownames_to_column("external_gene_name") %>%
    dplyr::filter(external_gene_name %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("external_gene_name", "patientID", "logCPM")) %>%
    left_join(dge$samples, by = "patientID")  %>%
    left_join(dge$genes, by = "external_gene_name")
 
  
  # Plot expression
  
  expressionPlot <- ggplot(expression, aes(x=primary.subtype, y=logCPM, fill=primary.subtype)) +  
    geom_point(aes(color=primary.subtype), position = position_jitterdodge()) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=9))+
    labs(x="Subtype") +
    #facet_wrap(.~external_gene_name, scales='free_y') +
    ggtitle("Expression of selected genes across subtypes")
  
  return(expressionPlot)
}


