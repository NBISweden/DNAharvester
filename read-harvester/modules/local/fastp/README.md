# Notes on FastP module implementation

Kept option to process single-end data (pairedness is automatically detected from samplesheet)

Reads that do not pass filtering and trimming are discarded

Adapter auto-detection for PE data is enabled:
- For PE data, the adapter sequence auto-detection is disabled by default since the adapters 
  can be trimmed by overlap analysis. However, you can specify --detect_adapter_for_pe to enable it.
- For PE data, fastp will run a little slower if you specify the sequence adapters or enable adapter 
  auto-detection, but usually result in a slightly cleaner output, since the overlap analysis may 
  fail due to sequencing errors or adapter dimers.
- Consider the option to provide an adapter fasta file instead of auto-detection

Enabled base correction in overlapped regions (only for PE data; --correction) with default settings
- --overlap_len_require            the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
- --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
- --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])