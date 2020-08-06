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


sample_name <- 'Test'
sinbad_object = construct_sinbad_object(working_dir = working_dir,
                                        raw_fastq_dir = raw_fastq_dir,
                                        demux_index_file = demux_index_file,
                                        sample_name = sample_name)

#process_sample_wrapper(sinbad_object)

sinbad_object = wrap_demux_fastq_files(sinbad_object)

sinbad_object = wrap_trim_fastq_files(sinbad_object)






