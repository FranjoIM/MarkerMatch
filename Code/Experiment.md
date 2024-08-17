### Environments
#### Local
Local analyses were ran on a Windows machine. Specifications:
|Component| Value|
|---|---|
|Processor|	12th Gen Intel(R) Core(TM) i9-12900KS   3.40 GHz|
|Installed RAM|	128 GB (128 GB usable)|
|System type|	64-bit operating system, x64-based processor|
|Edition|	Windows 11 Pro|
|Version|	23H2|
|OS build|	22631.4037|

#### Cloud
Cloud analyses were ran oh HiPerGator, a Linux cloud computing service, using SLURM scheduler to schedule/run jobs.

### Prepare Data

Ran locally. Using Illumina GenomeStudio, prepare data by clustering microarray genotyping samples _de novo_. Once done, exclude non-autosomes using the filter 
function, and export following files:
1. Illumina Final Report (Columns: `SNP Name`, `Sample ID`, `Chr`, `Position`, `B Allele Freq`, `Log R Ratio`),
2. PLINK Input Report (Leaving default options, except `UseForwardStrand` should be set to `True`),
3. Copy-numer metrics (CNM file, export `Name`, `Chr`, `Position` and `Log R Ratio` for all samples columns from the `Full Data Table` tab as tab-delimeted file),
4. Manifest file (MAN file, export `Name`, `Chr`, `Position` and `Minor Freq` columns from the `SNP Table` tab as a tab-delimeted file).

### Complete the Manifest Files

In R, run on a local machine, use thr [ManifSum.R](ManifSum.R) function to suummarize the copy-number metrics (CNM) export into copu-number summary (CNS) file.
