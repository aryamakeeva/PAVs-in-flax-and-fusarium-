# 🌿🧬 Presence-absence variations in Fusarium oxysporum and Linum usitatissimum 🧬🌿

> **Arina Makeeva**  
> Research Project, Bioinformatics Institute
> 📬 arina.makeeva@bio.msu.ru  
> 🧪 Supervisor: A. Samsonova, A.Kanapin, V.Stanin

---

## About this project

Plant genomes are dynamic landscapes shaped by genome duplications, gene loss, and the activity of mobile elements. These processes lead to structural variants, including **Presence/Absence Variations (PAVs)** — genomic regions present in some individuals but completely absent in others (Figure 1).

![Alt text](./imgs/Figure1.png)

*Figure 1: Presence/Absence Variations (PAVs) in plant genomes*

In this study, we investigate the role of PAVs in shaping:

- **Virulence** in the pathogenic fungus *Fusarium oxysporum*  
- **Resistance** in two *Linum usitatissimum* cultivars

*Keywords:* PAV, structural variation, Flax, Fusarium

---

## Datasets

We analysed:

- **12 genomes of *F. oxysporum*** (PRJNA630722, PRJNA721899) [1]
- **Two *Linum usitatissimum* cultivars** [2]:  
  - *AT* (resistant)  
  - *LM98* (susceptible) 

---

## PAV Analysis Pipeline 

![](./imgs/Figure2.png)

The required software environment is defined in the `environment.yaml` file.

```bash
conda env create -f environment.yaml
conda activate scanPAV
```

### 1. Scan for PAVs

A custom script based on the [repository example](https://github.com/wtsi-hpag/scanPAV). 
Location: `./code/1_PAV.sh`

---

### 2. Index the Reference Genome

```bash
bwa index reference.fasta
```

---

### 3. Align Genomes to Reference

```bash
bwa mem "$reference_genome" "$fasta_file" > "$output_dir/${base_name}.sam"
```

### 4. Convert SAM → Sorted BAM + Index

```bash
bam_file="./alignment_bam/$(basename "${sam_file%.sam}.bam")"
samtools view -Sb "$sam_file" > "$bam_file"

sorted_bam_file="./alignment_bam/$(basename "${sam_file%.sam}_sorted.bam")"
samtools sort -o "$sorted_bam_file" "$bam_file"

samtools index "$sorted_bam_file"
```
### 5. Annotate Genes

```bash
bedtools intersect -a "$reference_gff" -b "$bam_file" -wa -wb > "$output_dir/${base_name}_annotated_genes.txt"
```

### 6. Generate PAV Matrix

A binary matrix is created from gene lists to summarise gene presence in each pairwise comparison (`./code/2_pav_matrix.py`):

- `1` — gene is unique to genome A
- `0` - gene is shared between genomes A and B

The downstream analysis and data visualisation were performed using the following Python packages: `pandas`, `numpy`, `scipy.stats`, `seaborn`, `matplotlib`.

Example of the output:

```bash
🧬 Matrix saved to gene_presence_matrix.tsv
🔢 Total genes: 2487 | Comparisons: 156
      COMP_1  COMP_2  COMP_3  ...
GeneA       1       0       1
GeneB       1       1       0
...
```
---

## Results

- Among 12 Fusarium genomes we identified 233 highly divergent genes (Figure 2.A), and examined domain differences (Figure 2.B);

![Alt text](./imgs/Figure4.png)

*Figure 2. A - Heatmap of the most variable genes across the analysed Fusarium strains, based on presence/absence variation. B - GO enrichment analysis of protein domains associated with the top 5% of genes ranked by variation score. Bars represent significantly enriched functional categories; values indicate p-values calculated using the hypergeometric test*

- Compared the gene content between a susceptible and a resistant flax variety (Figure 3.A)and performed GO enrichment analysis. Genes specific to the susceptible variety are associated with processes related to cell degradation and viral activity, whereas genes unique to the resistant variety are enriched in cellular defence mechanisms, including mitochondrial regulation and calcium signalling (Figure 3.B).

![Alt text](./imgs/Figure3.png)

*Figure 3. A - Venn diagram showing the number of unique and shared genes between the susceptible and resistant flax varieties. B - Gene Ontology (GO) enrichment analysis of variety-specific genes. Bars represent enriched biological processes; values indicate p-values calculated via hypergeometric test*

## Future plans

We are going to evaluate genomic context of identified PAVs. Check for updates ⭐️

---

## Literature

1. Logachev, A., Kanapin, A., Rozhmina, T., Stanin, V., Bankin, M., Samsonova, A., Orlova, E., & Samsonova, M. (2024). Pangenomics of flax fungal parasite Fusarium oxysporum f. sp. lini. Frontiers in plant science, 15, 1383914. https://doi.org/10.3389/fpls.2024.1383914 
2. Rozhmina, T. A., Kanapin, A. A., Bankin, M. P., & Samsonova, M. G. (2024). Identification of Two QTLs Controlling Flax Resistance to Fusarium Wilt. Biophysics, 69(1), 57–62. https://doi.org/10.1134/S0006350924700076 






