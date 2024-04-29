# Mats-Blue-Hole-Maldives

Analysis for the paper:
## Genome-resolved metagenomics revealed novel microbial taxa and distinct metabolism from giant filamentous microbial mats inhabiting anoxic deep reefs of a Maldivian Blue Hole
___Lapo Doni 1,2, Annalisa Azzola 1,2, Caterina Oliveri 1, Emanuele Bosi 1,2, Manon Auguste 1,2, Carla Morri 1,3, Carlo Nike Bianchi 1,3, Monica Montefalcone 1,2, Luigi Vezzulli 1,2*___


_1 Department of Earth, Environmental and Life Sciences (DiSTAV), University of Genoa, Corso Europa 26, 16132 Genoa, Italy.
2 National Biodiversity Future Center, Palermo, Italy
3 Department of Integrative Marine Ecology (EMI), Stazione Zoologica Anton Dohrn - National Institute of Marine Biology, Ecology and Biotechnology, Genoa Marine Centre (GMC), Villa del Principe, Piazza del Principe 4, 16126 Genoa, Italy_

 

## MAGs with [metaWRAP](https://github.com/bxlab/metaWRAP) 
```
conda activate metawrap-env

mkdir READ_QC

metawrap read_qc -1 PELONI_SPSEA-07-22-N1715_S6_L007_R1_001.fastq-005.gz -2 PELONI_SPSEA-07-22-N1715_S6_L007_R2_001.fastq-010.gz -t 24 -o READ_QC/

for i in PELONI_SPSEA-07-22-N1715_S6_L007_R1_001.fastq-005.gz_val_1.fq.gz
do 
	prefix=$(basename $i _SPSEA-07-22-N1715_S6_L007_R1_001.fastq-005.gz_val_1.fq.gz)
mv PELONI_SPSEA-07-22-N1715_S6_L007_R1_001.fastq-005.gz_val_1.fq.gz ../CLEAN_READS/ALL_READS_1.fastq.gz
mv PELONI_SPSEA-07-22-N1715_S6_L007_R2_001.fastq-010.gz_val_2.fq.gz ../CLEAN_READS/ALL_READS_2.fastq.gz
done

metawrap assembly -1 CLEAN_READS/ALL_READS_1.fastq.gz -2 CLEAN_READS/ALL_READS_2.fastq.gz  -t 30 -o ASSEMBLY

metawrap binning -o INITIAL_BINNING -t 10 -a ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/*fastq.gz
metawrap bin_refinement -o BIN_REFINEMENT -t 10 -A INITIAL_BINNING/metabat2_bins/ -B INITIAL_BINNING/maxbin2_bins/ -C INITIAL_BINNING/concoct_bins/ -c 50 -x 10
metawrap quant_bins -b BIN_REFINEMENT/metawrap_50_10_bins -o QUANT_BINS -a ASSEMBLY/final_assembly.fasta CLEAN_READS/*fastq.gz

```




### MAGs taxonomy with  with [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)

```
conda activate gtdbtk-2.1.0
mkdir gtdbtk
gtdbtk classify_wf -x .fa --genome_dir . --out_dir gtdbtk --scratch_dir gtdbtk/scratch_dir --cpus 15
```

### MAGs QC with  with [CheckM](https://github.com/Ecogenomics/CheckM)

```
checkm lineage_wf -t 10 -x fa . ./checkm -f ./checkm/checkm_lineage_wf_overview_qa.txt --tab_table 
```





## Phylogenomic analysis of Dehalococcoidales MAGs
###  order Dehalococcoidales  tax id 1202465

 
### Dereplication with [dRep](https://github.com/MrOlm/drep)
 ```
dRep dereplicate ncbi_derep/ -g *.fna --contamination 10 
 ```
### ANI with [pyani](https://github.com/widdowquinn/pyani)
 ```
average_nucleotide_identity.py -i . -o ANIb_output -m ANIb -g -v 2>&1 > ANIb_output/error.log
 ```
### AAI with [compareM](https://github.com/donovan-h-parks/CompareM)
 ```
comparem aai_wf . AAI  --file_ext fa
 ```
### Tree with [GToTree](https://github.com/AstrobioMike/GToTree) & [iqtree](https://github.com/Cibiv/IQ-TREE)
```
GToTree  -f fasta_files.txt -H Bacteria  -j 30
iqtree -s Aligned_SCGs.faa -B 1000 -alrt 1000
 ```

## Metabolic analysis of all MAGs with [METABOLIC](https://github.com/AnantharamanLab/METABOLIC)
 ```
conda activate METABOLIC_v4.0
perl ~/METABOLIC_running_folder/METABOLIC/METABOLIC-C.pl  -t 30 -in-gn metawrap_bins -r reads.txt -o METABOLIC-C_out -m-cutoff 0.75 -kofam-db full
 ```


## LUCA-likeness with [Melange](https://github.com/sandragodinhosilva/melange)
```
conda activate snakemake
snakemake --use-conda --cores 30
Outuput->>Cog_PA.csv & Cog_description.csv
```

in R

```
library(tidyverse)


cog_luca <- read_lines("cog_luca.txt")  
cog_pa <- read_csv("Cog_PA.csv")  
Cog_description <- read_csv("Cog_description.csv")
 
COG_LUCA_in_bins <- cog_pa %>%
  filter(index %in% cog_luca)
dim(COG_LUCA_in_bins)
COG_LUCA_with_description <- COG_LUCA_in_bins %>%
  left_join(Cog_description, by = "index")
dim(COG_LUCA_with_description)



luca_cog_counts <- cog_pa %>%
  filter(index %in% cog_luca) %>%
  select(-index) %>%
  summarise(across(everything(), sum))
total_cogs_per_genome <- cog_pa %>%
  select(-index) %>%
  summarise(across(everything(), sum))
luca_cog_percentages <- luca_cog_counts / total_cogs_per_genome * 100
print(luca_cog_percentages,digits = 3)
 ```

## Detection of eukaryotes in the Mat with [EukDetect](https://github.com/allind/EukDetect)
```
conda activate eukdetect
snakemake --snakefile rules/eukdetect_eukfrac.rules --configfile PELONIdefault_configfile.yml --cores 20 runall 
```






