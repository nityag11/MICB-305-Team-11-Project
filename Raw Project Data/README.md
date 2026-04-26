
qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path /datasets/project_2/MS/ms_manifest.tsv \
  --output-path /demux_seqs.qza

qiime demux summarize \
  --i-data demux_seq.qza \
  --o-visualization demux.qzv

qiime demux summarize  
 --i-data demux_seqs.qza  
 --o-visualization demux.qzv

 qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 151 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

qiime metadata tabulate \
  --m-input-file stats.qza \
  --o-visualization stats.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /datasets/project_2/MS/corrected_ms_metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

  qiime feature-classifier extract-reads \
  --i-sequences /datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-trunc-len 151 \
  --o-reads ref-seqs-trimmed.qza
