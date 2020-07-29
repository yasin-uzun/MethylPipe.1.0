library(SINBAD)
SINBAD::test()

config_dir = '/mnt/isilon/tan_lab/uzuny/projects/jamboree/config_files/'
read_configs(config_dir)

example_dir = '/mnt/isilon/tan_lab/uzuny/projects/jamboree/data/example/'
raw_fastq_dir = paste0(example_dir, '/reads/')
demux_index_file = paste0(example_dir, '/demux_index/')
working_dir = paste0(example_dir, '/working_dir/')

process_sample_wrapper(raw_fastq_dir = raw_fastq_dir,
                       demux_index_file = demux_index_file,
                       working_dir = working_dir,
                       sample_name = 'Test')


