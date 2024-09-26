# Extracting iBAQ values from DIA data from DIA-NN 1.9.1 search results

The following packages are required to run this script:

```
library(arrow)
library(here)
library(tidyverse)
library(DIAgui)
```

### Load and filter the report.parquet file from DIANN search

```
diann_report <- arrow::read_parquet("report.parquet") %>%
    dplyr::filter(PG.MaxLFQ.Quality >= 0.75 & Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
    dplyr::mutate(File.Name = Run)
```

### transforming the diann_report in a format compatible with DIAgui
We get the gene names and peptides from the diann_report.

⚠️ the function `diann_matrix` used here is not the same from `diann` package. It comes from the `DIAgui` package.

```
unique_genes <- DIAgui::diann_matrix(diann_report,
    id.header = "Precursor.Id",
    quantity.header = "Precursor.Normalised",
    proteotypic.only = TRUE,
    pg.q = .01
)
```

### Get the protein sequences and iBAQ values
The fasta file needs to be in the same folder as the script, so we can use the `getallseq` function properly.

```
sequence_list <- DIAgui::getallseq(
    spec = "Homo sapiens",
    pr_id = unique_genes$Protein.Group,
    bank_name = "UP000005640_9606.fasta",
    fasta_file = TRUE
    )

# get the iBAQ values for each protein
iBAQ <- DIAgui::get_iBAQ(
  unique_genes,
  proteinDB = sequence_list,
  id_name = "Protein.Group",
  ecol = 7:21,       # here we need to specify the indexes of the columns containing the intensity values
  peptideLength = c(5, 35),
  nbMiscleavages = 0,
  proteaseRegExp = getProtease("trypsin"),
  log2_transformed = TRUE,
  keep_original = FALSE
)

# save the file in the working directory
write_tsv(iBAQ, "DIANN_iBAQ_values.tsv")
```

## Additional features that can be interesting

### Annotate the protease inhibitors and peptidases found in the dataset.
Use the **IDmapping** in the [**UniProtKB**](https://www.uniprot.org) to recover the protein families information for each ID.
Dowload the query as **.tsv** file with the usefull information.

```
annotated_iBAQ <- iBAQ %>%
  left_join(read_tsv("idmapping_2024_09_19.tsv"), 
                        by = join_by(Protein.Group == From)) %>%
  dplyr::mutate(protease_or_inhibitor = case_when(
    str_detect(`Protein families`, "Protease inhibitor") ~ "Protease inhibitor",
    str_detect(`Protein families`, "Peptidase") ~ "Protease",
    TRUE ~ "General"
  ))

View(annotated_iBAQ)
```

### plot the iBAQ values idstribution higlihting the protease inhibitors and peptidases

```
boxplot_iBAQ <- annotated_iBAQ %>%
  dplyr::select(Protein.Group, Genes, protease_or_inhibitor, starts_with("iBAQ")) %>%
  pivot_longer(-c(Protein.Group, Genes, protease_or_inhibitor), names_to = "sample", values_to = "iBAQ") %>%
  ggplot(aes(x = protease_or_inhibitor, y = iBAQ, fill = protease_or_inhibitor)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("grey60", "firebrick3", "dodgerblue2")) +
  facet_wrap(~sample, scales = "free_y", ncol = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(y = "log2(iBAQ)", x = NULL, fill = NULL)

ggsave("boxplot_iBAQ.png",
    boxplot_iBAQ, width = 10, 
    height = 10, dpi = 300)

annotated_iBAQ %>% 
  dplyr::count(protease_or_inhibitor)
```
