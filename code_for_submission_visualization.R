#######################################################
## UMAP Plot by Cluster or Gene-based Cell Type
#######################################################

colors <- c("#FFE119", "#000075", "#9A6324", "#17BECF", "#F032E6", "#B7E544", "#DDB673", "#A9A9A9", "#FFBB78", "#FABED4",
            "#808000", "#469990", "#911EB4", "#AEC7E8", "#BCBD22", "#FF7F0E", "#4363D8", "#E6194B", "#3CB44B", "#DCBEFF",
            "#800000", "#98DF8A")

DimPlot(obj.srt, reduction = "umap", group.by = "cluster_name", cols = colors)

# Gene-based Neuron classification: Tbr1 (Excitatory), Gad1 (Inhibitory)
genes <- c("Tbr1", "Gad1")
gene_data <- FetchData(obj.srt, vars = genes) %>% mutate(condition = obj.srt$condition)
cells_exn <- gene_data %>% filter(Tbr1 > 0 & Gad1 == 0) %>% rownames()
cells_inn <- gene_data %>% filter(Gad1 > 0 & Tbr1 == 0) %>% rownames()

umap_data <- Embeddings(obj.srt, "umap") %>% as.data.frame() %>% mutate(Neuron = "NA", condition = obj.srt$condition)
umap_data[cells_exn, "Neuron"] <- "ExN: Tbr1 > 0 & Gad1 == 0"
umap_data[cells_inn, "Neuron"] <- "InN: Gad1 > 0 & Tbr1 == 0"

ggplot(umap_data, aes(UMAP_1, UMAP_2, color = Neuron)) +
  geom_point(size = 0.2) +
  facet_wrap(~condition) +
  scale_color_manual(values = c("grey88", "#40B0A6", "#E1BE6A")) +
  theme_classic()

#######################################################
## Violin Plot of Cell-Type Marker Genes
#######################################################

marker_genes <- c("Aspm", "Hes1", "Ednrb", "Eomes", "Fezf2", "Tbr1", "Satb2", "Dlx2", "Mcm2", "Cenpa", "Six3", "Foxp2",
                  "Nkx2-1", "Lhx6", "Adarb2", "Magel2", "Shox2", "Reln", "Cx3cr1", "Emcn")

Idents(obj.srt) <- "RNA_snn_res.0.4"
cluster_levels <- c(17, 14, 13, 18, 3, 1, 12, 10, 8, 9, 7, 16, 5, 2, 4, 0, 6, 11, 15, 19, 20, 21)
obj.srt$RNA_snn_res.0.4 <- factor(obj.srt$RNA_snn_res.0.4, levels = cluster_levels)

VlnPlot(obj.srt, features = marker_genes, group.by = "RNA_snn_res.0.4", idents = cluster_levels, stack = TRUE, flip = TRUE)

#######################################################
## Split Violin Plots for Selected Genes by Condition
#######################################################

genes <- c("Tbr1", "Gad1")  # Selected genes
vln_colors <- c("grey88", "red")

for (gene in genes_subset) {
  p <- VlnPlot(obj.srt, features = gene, group.by = "Cell_Type", split.by = "condition",
               cols = vln_colors, pt.size = 0.1, alpha = 0.3)
  print(p)
}

#######################################################
## Volcano Plot Function (DEG Visualization)
#######################################################

volcanoplot <- function(deg, title, fc = 1.5, pval = 0.05) {
  deg_sig <- deg %>%
    filter(DE %in% c("UP", "DN")) %>%
    filter(abs(avg_log2FC) > log2(fc)) %>%
    filter(p_val_adj < pval)

  genes_to_label <- deg_sig$gene

  p <- ggplot(deg, aes(avg_log2FC, -log10(p_val_adj), color = DE)) +
    geom_point(size = 0.2, alpha = 0.8) +
    scale_color_manual(values = c("red", "blue", "grey90")) +
    geom_vline(xintercept = c(-log2(fc), log2(fc)), color = "grey90") +
    geom_hline(yintercept = -log10(pval), color = "grey90") +
    theme_bw() +
    ggrepel::geom_text_repel(data = filter(deg, gene %in% genes_to_label),
                             aes(label = gene), size = 2, max.overlaps = 10, max.time = 1, max.iter = 1000) +
    ggtitle(title)
  print(p)
}