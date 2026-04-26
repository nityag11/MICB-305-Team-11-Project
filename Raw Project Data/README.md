QIIME2 Export Files
Last login: Sat Jan 31 15:21:39 2026 from 128.189.118.161
(qiime2-amplicon-2025.4) root@stu-6277:~# cd /
(qiime2-amplicon-2025.4) root@stu-6277:/# ls
bin   data	dev  home  lib64  mnt  proc  run   srv	tmp  var
boot  datasets	etc  lib   media  opt  root  sbin  sys	usr  work
(qiime2-amplicon-2025.4) root@stu-6277:/# cd work
(qiime2-amplicon-2025.4) root@stu-6277:/work# ls
jimmy_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work# cd jimmy_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work/jimmy_group_assignment# ls
aligned-rep-seqs.qza	     rooted-tree.qza  taxa-bar-plots.qzv
alpha-rarefaction.qzv	     stats.qza	      taxonomy.qza
demux.qza		     stats.qzv	      taxonomy.qzv
masked-aligned-rep-seqs.qza  table.qza	      unrooted-tree.qza
rep-seqs.qza		     table.qzv	      visualization.qzv
(qiime2-amplicon-2025.4) root@stu-6277:/work/jimmy_group_assignment# cd ../../
(qiime2-amplicon-2025.4) root@stu-6277:/# ls
bin   data	dev  home  lib64  mnt  proc  run   srv	tmp  var
boot  datasets	etc  lib   media  opt  root  sbin  sys	usr  work
(qiime2-amplicon-2025.4) root@stu-6277:/# cd project_2
-bash: cd: project_2: No such file or directory
(qiime2-amplicon-2025.4) root@stu-6277:/# cd datasets
(qiime2-amplicon-2025.4) root@stu-6277:/datasets# cd project_2
(qiime2-amplicon-2025.4) root@stu-6277:/datasets/project_2# ls
alcohol     dog		  fish		  infant      prefetch	     vaccine
anemia	    dorms	  fmt		  IVF	      smoking	     voles
colombia    duck_flu	  gastric_cancer  melanoma    soil	     wetlands
COVID	    dysautonomia  hiseas	  MS	      space_station  zoo
depression  evelyn	  hiv		  nasa	      starch
diabetes    farm	  human_ibd	  parkinsons  tanzania
(qiime2-amplicon-2025.4) root@stu-6277:/datasets/project_2# cd MS
(qiime2-amplicon-2025.4) root@stu-6277:/datasets/project_2/MS# ls
corrected_ms_metadata.tsv  ms_manifest.tsv  ms_metadata.tsv  seqs  test
(qiime2-amplicon-2025.4) root@stu-6277:/datasets/project_2/MS# cd ../../..
(qiime2-amplicon-2025.4) root@stu-6277:/# ls
bin   data	dev  home  lib64  mnt  proc  run   srv	tmp  var
boot  datasets	etc  lib   media  opt  root  sbin  sys	usr  work
(qiime2-amplicon-2025.4) root@stu-6277:/# cd work
(qiime2-amplicon-2025.4) root@stu-6277:/work# ls
jimmy_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work# mkdir nitya_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work# ls
jimmy_group_assignment	nitya_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work# cd nitya_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/MS/ms_manifest.tsv \
  --output-path /demux_seqs.qza
Imported /datasets/project_2/MS/ms_manifest.tsv as SingleEndFastqManifestPhred33V2 to /demux_seqs.qza
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# ls
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# cd ../
(qiime2-amplicon-2025.4) root@stu-6277:/work# ls
jimmy_group_assignment	nitya_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work# cd nitya_group_assignment
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# ls
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# mv /demux_seqs.qza .
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# ls
demux_seqs.qza
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# qiime demux summarize \
  --i-data demux_seq.qza \
  --o-visualization demux.qzv
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# qiime demux summarize   --i-data demux_seqs.qza   --o-visualization demux.qzv
Saved Visualization to: demux.qzv
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# ls
demux.qzv  demux_seqs.qza
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 151 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza
Saved FeatureTable[Frequency] to: table.qza
Saved FeatureData[Sequence] to: rep-seqs.qza
Saved SampleData[DADA2Stats] to: stats.qza
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv
Saved Visualization to: stats.qzv
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /datasets/project_2/MS/corrected_ms_metadata.tsv
Saved Visualization to: table.qzv
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
Saved Visualization to: rep-seqs.qzv
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# ls
demux.qzv	rep-seqs.qza  stats.qza  table.qza
demux_seqs.qza	rep-seqs.qzv  stats.qzv  table.qzv
(qiime2-amplicon-2025.4) root@stu-6277:/work/nitya_group_assignment# client_loop: send disconnect: Broken pipe
