library(dplyr)
library(data.table)
library(GenomicRanges)

read_region_annot <- function(region_annot_file, format_file, type_filtering =  'none')
{
    df_region_annot = read.table(region_annot_file, sep ='\t')
    df_format = read.table(format_file, sep ='\t')
    colnames(df_region_annot) = as.character(df_format$V2)
    head(df_region_annot)
    dim(df_region_annot)
    #df_region_annot %>% distinct(region_name, .keep_all = TRUE)

    df_region_annot <- dplyr::distinct(df_region_annot, region_name, .keep_all = TRUE)
    dim(df_region_annot)
    head(df_region_annot)
    length(unique(df_region_annot$region_name))


    if(type_filtering != 'none')
    {
      df_region_annot  = df_region_annot[df_region_annot$region_type == type_filtering, ]
    }

    rownames(df_region_annot) =  df_region_annot$region_name

    head(df_region_annot)

    return(df_region_annot)
}

get_promoters <- function(df_gene_annot)
{

    head(df_gene_annot)
    df_promoters = df_gene_annot
    head(df_promoters)
    attach(df_gene_annot)
    df_promoters$start[strand == '+'] = df_gene_annot$start[strand == '+'] - 2000
    df_promoters$end[strand == '+'] = df_gene_annot$start[strand == '+'] + 500

    df_promoters$end[strand == '-'] = df_gene_annot$end[strand == '-'] + 2000
    df_promoters$start[strand == '-'] = df_gene_annot$end[strand == '-'] - 500
    detach(df_gene_annot)
    head(df_promoters)

    df_promoters$start[df_promoters$start < 0] = 0

    return(df_promoters)
}

convert_to_granges <- function(df_region)
{

    gr_region <- with(df_region, GRanges(chrom,
                                           IRanges(start+1, end),
                                           strand = '*',
                                           region_name,
                                           region_type))

    return(gr_region)
}


compute_region_met_matrix <- function(df_region, methylation_calls_dir,
                                       methylation_type = 'CpG',
                                       min_call_count_threshold = 10)
{

  gr_region = convert_to_granges(df_region)
  setwd(methylation_calls_dir)
  pattern  = paste0(methylation_type, '_calls.*organism.cov.*')
  cov_files = list.files(methylation_calls_dir, pattern)

  result_list = list()


  cl <- makeCluster(num_cores, outfile="", type = 'SOCK')
  registerDoSNOW(cl)

  #for(i in  1:length(cov_files))
  #result_list <- foreach(i=1:length(cov_files)) %dopar%
  result_list <- foreach(i=1:length(cov_files)) %dopar%
  {
    library(data.table)
    library(GenomicRanges)

    print(i)
    cov_file = cov_files[i]
    cell_id = gsub('.organism.cov.gz', '', cov_file)
    cell_id = gsub(paste0(methylation_type,'_calls.'), '', cell_id)

    dt_cov = fread(paste0(methylation_calls_dir, cov_file) )
    colnames(dt_cov) = c('chrom', 'start', 'end', 'met_rate', 'met', 'demet')
    #dt_cov$chr = paste0('chr', dt_cov$chr)
    head(dt_cov)
    unique(dt_cov$chr)

    gr_cov <- with(dt_cov, GRanges(chrom, IRanges(start+1, end), strand = '*', met_rate, met, demet)  )
    gr_cov

    #dt_inter = data.table(intersect_bed(gr_region, gr_cov))
    #dim(dt_inter)
    #head(dt_inter)

    hits_obj <- findOverlaps(gr_region, gr_cov)
    class(hits_obj)
    da = as.data.frame(gr_region[queryHits(hits_obj)])
    db = as.data.frame(gr_cov[subjectHits(hits_obj)])
    dt_inter <- data.table(cbind(da, db))


    quant_cols = c('met', 'demet')

    dt_aggr <- dt_inter[, lapply(.SD, sum), by = .(region_name), .SDcols = quant_cols  ]
    length(unique(dt_aggr$region_name))
    head(dt_aggr)
    dim(dt_aggr)

    dt_aggr$call_count = dt_aggr$met + dt_aggr$demet

    dt_aggr$met_rate = dt_aggr$met / dt_aggr$call_count
    head(dt_aggr)

    dt_aggr$met_rate[dt_aggr$call_count < min_call_count_threshold] = NA

    df_aggr = data.frame(dt_aggr)
    head(df_aggr)

    df_region = data.frame(gr_region)
    head(df_region)
    length(df_region)

    df_merged = merge(df_region, df_aggr, by = "region_name",  all.x = T)
    head(df_merged)
    dim(df_merged)

    met_rate_vector = df_merged$met_rate
    names(met_rate_vector) = df_merged$region_name
    head(met_rate_vector)

    #result_list[[cell_id]] = met_rate_vector

    met_rate_vector
  }


  cell_ids = cov_files
  cell_ids = gsub('.organism.cov.gz', '', cell_ids)
  cell_ids = gsub(paste0(methylation_type,'_calls.'), '', cell_ids)
  names(result_list) = cell_ids

  met_rate_matrix = do.call("cbind", result_list)
  dim(met_rate_matrix)
  met_rate_matrix[1:5, 1:5]

  #stopCluster(cl)

  return(met_rate_matrix)

}



