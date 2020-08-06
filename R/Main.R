library(doSNOW)


test <- function()
{
  cat('SINBAD installation is ok.\n')
}
#

read_configs <- function(config_dir)
{
  print('Reading config.general.R for program paths...')
  source(paste0(config_dir, 'config.general.R') )
  print('Reading config.genome.R for genome paths...')
  source(paste0(config_dir, 'config.genome.R') )
  print('Reading config.project.R for prroject settings...')

  source(paste0(config_dir, 'config.project.R') )
  #system('echo $PATH')

}



#sample_name = 'Test'
#readRenviron(path = 'variables.env')
#system('source ~/.bashrc')

#image_file = paste0(working_dir, '/Image_2020_06_03.img')
#image_file = paste0(working_dir, '/Image_2020_06_11.img')
#save.image(image_file)

construct_sinbad_object <- function(raw_fastq_dir ,
                                    demux_index_file ,
                                    working_dir ,
                                    sample_name)
{


  main_log_dir <<- paste0(working_dir, '/logs/')
  demux_fastq_dir <<- paste0(working_dir, '/demux_fastq/')
  trimmed_fastq_dir <<- paste0(working_dir, '/trimmed_fastq/')
  alignment_dir <<- paste0(working_dir, '/alignments/')
  methylation_calls_dir <<- paste0(working_dir, '/methylation_calls/')
  summary_dir <<- paste0(working_dir, '/summary/')
  plot_dir <<- paste0(working_dir, '/plots/')

  dir.create(main_log_dir, showWarnings = F, recursive = T)
  dir.create(trimmed_fastq_dir, showWarnings = F, recursive = T)
  dir.create(alignment_dir, showWarnings = F, recursive = T)
  dir.create(methylation_calls_dir, showWarnings = F, recursive = T)
  dir.create(summary_dir, showWarnings = F, recursive = T)
  dir.create(plot_dir, showWarnings = F, recursive = T)
  dir.create(demux_fastq_dir, showWarnings = F, recursive = T)

  sinbad_object = list('main_log_dir' = main_log_dir,
                       'raw_fastq_dir' = raw_fastq_dir,
                       'demux_index_file' = demux_index_file,
                       'demux_fastq_dir' = demux_fastq_dir,
                       'trimmed_fastq_dir' =  trimmed_fastq_dir,
                       'alignment_dir' = alignment_dir,
                       'methylation_calls_dir'  = methylation_calls_dir,
                       'summary_dir' = summary_dir,
                       'plot_dir' = plot_dir,
                       'sample_name' = sample_name
                       )

  class(sinbad_object) = 'Sinbad'

  return(sinbad_object)


}


wrap_demux_fastq_files <- function(sinbad_object)
{
  #Demultiplex fastq files
  #TODO: Check the input fastq dir and demux file exists, give error and stop otherwise

  demux_fastq_files(sinbad_object$raw_fastq_dir,
                    sinbad_object$demux_index_file,
                    demux_index_length,
                    sinbad_object$demux_fastq_dir,
                    sinbad_object$main_log_dir)

  sinbad_object$df_demux_reports = read_demux_logs(main_log_dir)

  demux_summary_file = paste0(summary_dir, '/Demux_statistics.tsv')
  write.table(sinbad_object$df_demux_reports, file = demux_summary_file, sep = '\t', quote = F, row.names = F, col.names = T)

  #Count demuxd reads
  sinbad_object$demux_read_counts =  count_fastq_reads(demux_fastq_dir)

  return(sinbad_object)

}

wrap_trim_fastq_files <- function(sinbad_object)#Trim adapters
{
  trim_fastq_files(sinbad_object$demux_fastq_dir,
                   sinbad_object$trimmed_fastq_dir,
                   sinbad_object$main_log_dir)

  sinbad_object$trimmed_read_counts = count_fastq_reads(sinbad_object$trimmed_fastq_dir)

  plot_file = paste0(plot_dir, '/Preprocessing_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = T, title = sample_name)
  plot_preprocessing_results(sample_name = sinbad_object$sample_name,
                             demux_reports = sinbad_object$df_demux_reports,
                             demux_read_counts = sinbad_object$demux_read_counts,
                             trimmed_read_counts = sinbad_object$trimmed_read_counts)
  dev.off()

  return(sinbad_object)
}

wrap_align_sample <- function(sinbad_object)
{

  #Run aligner
  align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
               genomic_sequence_path,
               sinbad_object$alignment_dir, num_cores, mapq_threshold, main_log_dir)

  df_alignment_reports = process_bismark_alignment_reports(sinbad_object$alignment_dir)
  df_bam_read_counts = count_bam_files(alignment_dir)
  dim(df_alignment_reports)
  dim(df_bam_read_counts)


  df_alignment_stats = base::merge(df_alignment_reports, df_bam_read_counts, by = 0)
  dim(df_alignment_stats)
  head(df_alignment_stats)

  coverage_rates = compute_coverage_rates(alignment_dir)
  df_alignment_stats$coverage_rate = coverage_rates[as.character(df_alignment_stats$Cell_ID)]

  df_alignment_stats$Row.names = NULL
  rownames(df_alignment_stats) = df_alignment_stats$Cell_ID
  alignment_summary_file = paste0(summary_dir, '/Alignment_statistics.tsv')
  write.table(df_alignment_stats, file = alignment_summary_file, sep = '\t', quote = F, row.names = F, col.names = T)


  plot_file = paste0(plot_dir, '/Alignment_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = F, title = sample_name)
  plot_alignment_stats(sample_name, df_alignment_stats, coverage_rates)
  dev.off()

  df_org_split_reports = process_bismark_split_reports(methylation_calls_dir, genome_type = 'organism')
  df_lambda_split_reports = process_bismark_split_reports(methylation_calls_dir, genome_type = 'lambda')

  list_org_bias_reports = process_bismark_bias_reports(methylation_calls_dir, genome_type = 'organism')
  list_lambda_bias_reports = process_bismark_bias_reports(methylation_calls_dir, genome_type = 'lambda')

  plot_file = paste0(plot_dir, '/Met_call_statistics.eps')
  postscript(plot_file, paper = 'a4', horizontal = F, title = sample_name)
  plot_split_reports(df_org_split_reports, df_lambda_split_reports, list_org_bias_reports)
  dev.off()

  head(df_org_split_reports)

}

process_sample_wrapper <- function(sinbad_object)
{


  par(mfrow = c(2,2))

  sinbad_object = wrap_demux_fastq_files(sinbad_object)

  sinbad_object = wrap_trim_fastq_files(sinbad_object)




  #Call Methylation Sites
  call_methylation_sites_for_sample(alignment_dir, methylation_calls_dir, main_log_dir)

  #Read regions
  df_gene_annot = read_region_annot(gene_annot_file, format_file)
  head(df_gene_annot)
  df_promoters = get_promoters(df_gene_annot = df_gene_annot)
  df_100k_bins = read_region_annot(bins_100k_file, format_file)
  head(df_100k_bins)
  df_10k_bins = read_region_annot(bins_10k_file, format_file)
  head(df_10k_bins)


  #Quantify methylation in regions
  df_for_dim_red = df_10k_bins
  name_for_dim_red = '10k_bins'
  df_for_features = df_promoters
  name_for_features = 'Promoters'

  met_mat_for_dim_red = compute_region_met_matrix(df_for_dim_red, methylation_calls_dir,
                                                  methylation_type = methylation_type,
                                                  min_call_count_threshold = 10)
  met_mat_for_features = compute_region_met_matrix(df_for_features, methylation_calls_dir,
                                                   methylation_type = methylation_type,
                                                   min_call_count_threshold = 10)

  #Reduce dimensionality
  marker_genes = get_marker_genes('leuk')

  dim_red_list = reduce_dims_for_sample(met_mat_for_dim_red = met_mat_for_dim_red,
                         met_mat_for_features = met_mat_for_features,
                         name_for_dim_red = name_for_dim_red,
                         name_for_features = name_for_features,
                         methylation_calls_dir, plot_dir,
                         methylation_type = 'CpG',
                         min_call_count_threshold = 10,
                         features = marker_genes)

  #DM Analysis

  dm_stat_list_for_clusters = dm_stat_test_for_clusters(dim_red_list)
  dm_result_file = paste0(summary_dir, '/DM_Analysis.xlsx')
  write.xlsx(dm_stat_list_for_clusters$dm_result_list_with_summary, file = dm_result_file)

  all_dms = do.call("rbind", dm_stat_list_for_clusters$dm_result_list_with_summary[2:15] )
  all_dms_sorted = all_dms[order(all_dms$p.value), ]
  head(all_dms_sorted)

  all_dms_sig = all_dms_sorted[all_dms_sorted$adjusted.p.value < dmr_adj_p_value_cutoff, ]
  heatmap_dm_regions = as.character(all_dms_sorted$region)

  if(nrow(all_dms_sig) > dm_num_heatmap_regions)
  {
     heatmap_dm_regions = heatmap_dm_regions[1:dm_num_heatmap_regions]
  }
  #sig_regions = dm_stat_list_for_clusters$sig_ids
  sorted_clusters = sort(dim_red_list$clusters)

  met_mat = replace_nas_by_column_mean(dim_red_list$met_mat_for_features)
  sig_mat = as.matrix(met_mat[heatmap_dm_regions, names(sorted_clusters) ])
  df_annot = data.frame(cluster = sorted_clusters)

  head(df_annot)

  ph <- pheatmap::pheatmap(sig_mat, cluster_cols = F,
                           cluster_rows = T,
                           show_colnames = F, fontsize = 13,
                           annotation_col = df_annot,
                           color = viridis(100),
                           #annotation_colors = ann_colors,
                           fontsize_row = 9)

  plot_file = paste0(plot_dir, '/Heatmap.DMR_',dim_red_list$name_for_features,'.clusters.eps')
  ggsave(ph, filename = plot_file, device = 'eps', width = 20, height = 20, units = 'cm')

}


