
dir=pyscenic/TF 
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl


input_loom=OV-Endo.loom
ls $tfs  $feather  $tbl  


anaconda3/envs/pyscenic_1/bin/pyscenic grn \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
$input_loom  $tfs 


anaconda3/envs/pyscenic_1/bin/pyscenic ctx \
adj.sample.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname $input_loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 10  \
--mask_dropouts


anaconda3/envs/pyscenic_1/bin/pyscenic aucell \
$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 10

