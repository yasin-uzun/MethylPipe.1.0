#qsub -pe smp 10 -l h_vmem=16G -e /mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/Test.respublica.log -o /mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/Test.respublica.log  /mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/Test.respublica.sh 

module load R/3.6.2

Rscript /mnt/isilon/tan_lab/uzuny/projects/jamboree/package/p99/MethylProc/Test.R

