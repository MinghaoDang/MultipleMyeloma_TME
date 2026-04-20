#BSUB -J SCENIC
#BSUB -W 230:00
#BSUB -o /rsrch3/scratch/genomic_med/mdang1/logs
#BSUB -e /rsrch3/scratch/genomic_med/mdang1/errors
#BSUB -cwd /rsrch3/scratch/genomic_med/mdang1/Orlowski/TME/
#BSUB -q long
#BSUB -n 12
#BSUB -M 200
#BSUB -R rusage[mem=200]
#BSUB -B
#BSUB -N

# run the following code in R
# library('SCopeLoomR')
# library('SCENIC')
# CD8T.obj<-readRDS('MM.TME.235samples.CD8T.filter.scs.rds')
# exprMatrix<-GetAssayData(CD8T.obj,slot = 'counts')
# exprMatrix<-exprMatrix[apply(exprMatrix,1,function(x){sum(x>0)>0.01*ncol(exprMatrix)}),]
# cellInfo <- data.frame(CellType=CD8T.obj$cell.status)
# loom <- build_loom('/rsrch3/scratch/genomic_med/mdang1/Orlowski/TME/process/SCENIC/CD8T.filter.loom', dgem=exprMatrix)
# loom <- add_cell_annotation(loom, cellInfo)
# close_loom(loom)

sampleIndex=$LSB_JOBINDEX-1

celltypes=(CD16_Mono CD4T NK Bcell Macrophage Neutrophil DC CD14_Mono CD8T)
ct=${celltypes[$sampleIndex]}

cd /rsrch3/scratch/genomic_med/mdang1/Orlowski/TME/process/SCENIC

singularity exec -B /rsrch3:/rsrch3 /rsrch3/home/genomic_med/mdang1/software/pySCENIC/aertslab-pyscenic-0.11.0.sif \
arboreto_with_multiprocessing.py \
--num_workers 12 \
--output $ct.adj.tsv \
--method grnboost2 \
$ct.filter.loom \
/rsrch3/home/genomic_med/mdang1/data/SCENIC/hs_hgnc_tfs.txt \
--seed 777

singularity exec -B /rsrch3:/rsrch3 /rsrch3/home/genomic_med/mdang1/software/pySCENIC/aertslab-pyscenic-0.11.0.sif \
pyscenic ctx \
$ct.adj.tsv \
/rsrch3/home/genomic_med/mdang1/data/SCENIC/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname /rsrch3/home/genomic_med/mdang1/data/SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname $ct.filter.loom \
--mode "dask_multiprocessing" \
--output $ct.reg.csv \
--num_workers 12 \
--mask_dropouts

singularity exec -B /rsrch3:/rsrch3 /rsrch3/home/genomic_med/mdang1/software/pySCENIC/aertslab-pyscenic-0.11.0.sif \
pyscenic aucell \
$ct.filter.loom \
$ct.reg.csv \
--output $ct.filter.SCENIC.loom \
--num_workers 12
