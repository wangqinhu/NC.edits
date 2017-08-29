About
=====
This repository keeps the major scripts developed for analysis of A-to-I RNA editing in _Neurospora_


Notes
=====

- Comparing the RNA secondary structures of editing sites and random editable sites with 30-bp flanking sequence
```bash
# Real editing sites for RNA secondary structure analysis
perl ./cds_red_flank.pl 30 60
# Random editable sites for RNA secondary structure analysis
# Consider position -2 to +3
# Step 1
perl ./cds_rnd_sampos.pl
# Step 2
perl ./cds_rnd_flank.pl
# For RNA secondary structure analysis:
# See https://github.com/wangqinhu/red_loop
```

- Comparing the nonsyn/syn ratio of editing sites and random editable sites
```bash
# Consider position -2 to +3
# Total
perl ./level_group_ns_sam.pl 1
# Five groups
perl ./level_group_ns_sam.pl 5
# Statistic
perl ./level_group_ns_stat.pl
# Plot
R CMD BATCH ./level_group_ns.R
```

- Comparing the amino acid substitution types caused by editing sites and random editable sites
```bash
# Sampling and amnio acid changes calculation
perl ./sample_a2a.pl
# Plot
R CMD BATCH ./sample_a2a.nonsyn.R
```

- Study on the relationship of dN/dS and density of editing sites in the editing sites and random editable sites
```bash
# Convert random sample data to edit data format
perl ./bin/sam2red.pl
# Calculate dN/dS and fn,fs for real and random data (several minutes)
perl ./dnds_vs_fnfs.red_sam.pl
# Plot dN/dS vs. fn,fs for real data
R CMD BATCH ./dnds_vs_fnfs.red.R
# Prepare random data required for plotting
perl ./bin/sam_stat.pl
# Plot dN/dS vs. fn,fs for sample data
R CMD BATCH ./dnds_vs_fnfs.sam.R
```

- Study on the relationship of dN/dS and cumulative editing level
```bash
# Calculate dN/dS and cumulative editing level
perl ./dnds_vs_level.pl
# Plot dN/dS vs. cumulative editing level
R CMD BATCH ./dnds_vs_level.R
```

- Study on the relationship of gene expression and density of editing sites
```bash
# Calculate grouped expression level and fn,fs
perl ./expr_vs_fnfs.pl
# Plot expression (FPKM value) vs. cumulative editing level
R CMD BATCH ./expr_vs_fnfs.R
```

- Comparing the percentage of ancestral As replaced with G/C/T in _Neurospora tetrasperma_
```bash
# Real and random editing data obtained from above
# The conservation of real and random editing sites were calculated by the following scripts:
# https://github.com/wangqinhu/red_ortholog
# After processing the raw data by bin/ANA_extract.pl, two files were obtained:
# real: data/red.ANA.txt
# rand: data/rnd.ANA.txt
# Statistic
sh ./ANA_stat.sh
# Plot
R CMD BATCH ./ANA_plot.R
```

Citation
========
Liu H, Li Y, Chen D, Qi Z, Wang Q, Wang J, Jiang C, Xu JR. (2017) A-to-I RNA editing is developmentally regulated and generally adaptive for sexual reproduction in _Neurospora crassa_. Proc Natl Acad Sci USA, doi: [10.1073/pnas.1702591114][1].

[1]: http://www.pnas.org/cgi/doi/10.1073/pnas.1702591114
