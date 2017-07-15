ensgGeneMap_latest_Hum_Mou <-
function (infile, ensAttributes, org) 
{
    if (missing(infile)) {
        stop("ensmble ID/Probe ID data required.\n")
    }
    if (missing(ensAttributes)) {
        stop("Mention a measure of central tendceny parameter mean, median or sd.\n ")
    }
    if (missing(org)) {
        stop("Give organism name in ensmeble biomaRt format.\n")
    }
    cat("You have given all parameters as below \n\n")
    Sys.sleep(3)
    print(infile)
    print(ensAttributes)
    print(org)
    wdDir <- getwd()
    out2Dir <- paste0(wdDir, "/mappingOutput/")
    dir.create(out2Dir, showWarnings = FALSE)
    cat("Your output will be directed into following directory.\n\n ")
    print(out2Dir)
    if (org == "human" | org == "Human") {
        cat("Reading your dataset for gene ID mapping \n\n")
        file = basename(infile)
        filename = unlist(strsplit(file, "\\."))[[1]]
        data = read.delim(infile, header = TRUE, row.names = 1)
        dim(data)
        dim(data)
        mart.hs <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
            dataset = "hsapiens_gene_ensembl")
        genes <- getBM(attributes = c("hgnc_symbol", ensAttributes), 
            values = c(rownames(data)), filters = ensAttributes, 
            mart = mart.hs)
        gene_name = as.data.frame(cbind(genes$ensembl_gene_id, 
            genes$hgnc_symbol))
        gene_name$V2[gene_name$V2 == ""] <- NA
        gene_name = na.omit(gene_name)
        dim(gene_name)
        colnames(gene_name) = c("ID", "Gene_Name")
        write.table(gene_name, paste0(out2Dir, filename, "_ID_gene_name.txt", 
            sep = "\t"), quote = FALSE, row.names = FALSE)
        merge_probe_gene_name <- data.frame(Gene = gene_name$Gene_Name[match(rownames(data), 
            gene_name$ID)], data)
        dim(merge_probe_gene_name)
        final_gene_map = merge_probe_gene_name[, 1:ncol(merge_probe_gene_name)]
        dim(final_gene_map)
        with_gene_nameFile = paste0(out2Dir, "/", Sys.Date(), 
            "_", filename, "_with", "_gene_name_data.txt")
        write.table(final_gene_map, with_gene_nameFile, sep = "\t", 
            quote = FALSE, row.names = FALSE)
    }
    else if (org == "mouse" | org == "Mouse") {
        cat("Reading your dataset for gene ID mapping \n\n")
        file = basename(infile)
        filename = unlist(strsplit(file, "\\."))[[1]]
        data = read.delim(infile, header = TRUE, row.names = 1)
        dim(data)
        dim(data)
        mart.hs <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
            dataset = "mmusculus_gene_ensembl")
        genes <- getBM(attributes = c("mgi_symbol", ensAttributes), 
            values = c(rownames(data)), filters = ensAttributes, 
            mart = mart.hs)
        gene_name = as.data.frame(cbind(genes$ensembl_gene_id, 
            genes$mgi_symbol))
        gene_name$V2[gene_name$V2 == ""] <- NA
        gene_name = na.omit(gene_name)
        dim(gene_name)
        colnames(gene_name) = c("ID", "Gene_Name")
        write.table(gene_name, paste0(out2Dir, filename, "_ID_gene_name.txt", 
            sep = "\t"), quote = FALSE, row.names = FALSE)
        merge_probe_gene_name <- data.frame(Gene = gene_name$Gene_Name[match(rownames(data), 
            gene_name$ID)], data)
        dim(merge_probe_gene_name)
        final_gene_map = merge_probe_gene_name[, 1:ncol(merge_probe_gene_name)]
        dim(final_gene_map)
        with_gene_nameFile = paste0(out2Dir, "/", Sys.Date(), 
            "_", filename, "_with", "_gene_name_data.txt")
        write.table(final_gene_map, with_gene_nameFile, sep = "\t", 
            quote = FALSE, row.names = FALSE)
    }
}
