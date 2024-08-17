#### Prepare Data

Using Illumina GenomeStudio, prepare data by clustering microarray genotyping samples _de novo_. Once done, exclude non-autosomes using the filter 
function, and export following files:
1. Illumina Final Report (Columns: `SNP Name`, `Sample ID`, `Chr`, `Position`, `B Allele Freq`, `Log R Ratio`),
2. PLINK Input Report (Leaving default options, except `UseForwardStrand` should be set to `True`),
3. Copy-numer metrics (CNM file, export `Name`, `Chr`, `Position` and `Log R Ratio` for all samples columns from the `Full Data Table` tab as tab-delimeted file),
4. Manifest file (MAN file, export `Name`, `Chr`, `Position` and `Minor Freq` columns from the `SNP Table` tab as a tab-delimeted file).

### Complete the Manifest Files

