#library(devtools)
#install_github("yasin-uzun/SINBAD.1.0")

library(SINBAD)
SINBAD::test()

config_dir = '/mnt/isilon/tan_lab/uzuny/projects/jamboree/config_files/'
read_configs(config_dir)

example_dir = '/mnt/isilon/tan_lab/uzuny/projects/jamboree/data/example/'
raw_fastq_dir = paste0(example_dir, '/reads/')
demux_index_file = paste0(example_dir, '/demux_index.txt')
#working_dir = paste0(example_dir, '/working_dir/')

wenbao_project_dir = '/mnt/isilon/tan_lab/yuw1/projects/sinbad/data/example/'
working_dir = paste0(wenbao_project_dir, '/working_dir/')
dir.create(working_dir, recursive = T)


library(doSNOW)


wrap_align_sample <- function(sinbad_object)
{

  #Run aligner
  align_sample(read_dir = sinbad_object$trimmed_fastq_dir,
               genomic_sequence_path = genomic_sequence_path,
               alignment_dir = sinbad_object$alignment_dir,
               aligner = aligner,
               num_cores= num_cores,
               mapq_threshold =mapq_threshold,
               main_log_dir = sinbad_object$main_log_dir)


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


sample_name <- 'Test'
sinbad_object = construct_sinbad_object(working_dir = working_dir,
                                        raw_fastq_dir = raw_fastq_dir,
                                        demux_index_file = demux_index_file,
                                        sample_name = sample_name)

#process_sample_wrapper(sinbad_object)

#Demux
sinbad_object = wrap_demux_fastq_files(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.01.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print('Step 1')
print(sinbad_object)

#Trim
sinbad_object = wrap_trim_fastq_files(sinbad_object)

sinbad_object_file = paste0(working_dir, '/sinbad_object.02.rds')
saveRDS(sinbad_object,   file = sinbad_object_file)
sinbad_object = readRDS(sinbad_object_file)
print('Step 2')
print(sinbad_object)




