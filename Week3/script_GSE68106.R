library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ecoli2.db)
library(AnnotationDbi)
library(umap)

#PART C. PENGAMBILAN DATA DARI GEO 


#GEO (Gene Expression Omnibus) adalah database publik milik NCBI
#getGEO(): fungsi untuk mengunduh dataset berdasarkan ID GEO
#GSEMatrix = TRUE -> data diambil dalam format ExpressionSet
#AnnotGPL  = TRUE -> anotasi gen (Gene Symbol) ikut diunduh

gset <- getGEO("GSE68106", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

#PART D. PRE-PROCESSING DATA EKSPRESI 

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel
ex <- exprs(gset)

#Mengapa perlu log2 transformasi?
#Data microarray mentah memiliki rentang nilai sangat besar.
#Log2 digunakan untuk:
#1. Menstabilkan varians
#2. Mendekati asumsi model linear
#3. Memudahkan interpretasi log fold change

#quantile(): menghitung nilai kuantil (persentil)
#as.numeric(): mengubah hasil quantile (yang berupa named vector)
#menjadi vektor numerik biasa agar mudah dibandingkan
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

#LogTransform adalah variabel logika (TRUE / FALSE)
#Operator logika:
#>  : lebih besar dari
#| | : OR (atau)
#&& : AND (dan)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# PART E. DEFINISI KELOMPOK SAMPEL 
# 1. Mengambil metadata
group_info <- pData(gset)[["source_name_ch1"]]

# 2. Mengubah menjadi format valid R (menghilangkan spasi jadi titik)
groups <- make.names(group_info)

# 3. Mengubah menjadi Factor
gset$group <- factor(groups)

# 4. PROSES PENYINGKATAN NAMA (Renaming Levels)
# Kita definisikan nama baru sesuai urutan abjad dari levels yang lama
# Level [1]: bacteria...ciprofloxacin... (Perlakuan)
# Level [2]: bacteria...for.30.hours (Kontrol)

levels(gset$group) <- c("Cipro_Treatment", "Control_30h")

# 5. Verifikasi hasil
cat("Daftar grup yang sudah disingkat:\n")
print(levels(gset$group))

cat("\nJumlah sampel per grup:\n")
print(table(gset$group))

# 6. (Opsional) Menyimpan informasi grup ke variabel sederhana untuk analisis limma
group_list <- gset$group


# PART F. DESIGN MATRIX (KERANGKA STATISTIK)
# 1. Membuat matriks desain
# Menggunakan gset$group yang sudah kita singkat levelnya (Cipro_Treatment & Control_30h)
design <- model.matrix(~0 + gset$group)

# 2. Memberi nama kolom agar rapi
# colnames(design) sekarang akan berisi "Cipro_Treatment" dan "Control_30h"
colnames(design) <- levels(gset$group)

# 3. Menentukan perbandingan biologis (Kontras)
# Pastikan urutannya: [Grup Perlakuan] - [Grup Kontrol]
grup_treatment <- colnames(design)[1] # Cipro_Treatment
grup_control   <- colnames(design)[2] # Control_30h

# 4. Membuat formula kontras
contrast_formula <- paste(grup_treatment, "-", grup_control, sep="")

cat("Design Matrix (5 baris pertama):\n")
print(head(design, 5))

cat("\nKontras yang akan dianalisis:\n")
print(contrast_formula)


#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
# 1. Membuat Design Matrix
# Menggunakan variabel gset$group yang sudah kita singkat namanya tadi
design <- model.matrix(~ 0 + gset$group)

# Memberikan nama kolom agar lebih mudah dibaca (menghapus prefix "gset$group")
colnames(design) <- levels(gset$group)

# 2. Membangun Model Linear
# ex adalah matriks ekspresi gen yang sudah dinormalisasi
fit <- lmFit(ex, design)

# 3. Mendefinisikan Kontras (Perbandingan)
# Di sini kita membandingkan Cipro_Treatment vs Control_30h
# Format: Treatment - Control
contrast_formula <- "Cipro_Treatment - Control_30h"
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

# 4. Menerapkan Kontras ke Model
fit2 <- contrasts.fit(fit, contrast_matrix)

# 5. Empirical Bayes (Menghitung Statistik Signifikansi)
fit2 <- eBayes(fit2)

# 6. Mengambil Hasil Akhir DEG
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",    # Berdasarkan B-statistic (log-odds gen terekspresi diferensial)
  number = Inf,     # Mengambil semua gen
  p.value = 0.01    # Hanya yang FDR < 0.01
)

# Menampilkan 6 gen teratas
head(topTableResults)

#PART H. ANOTASI NAMA GEN 

#Penting:
#Pada data microarray Affymetrix, unit analisis awal adalah PROBE,
#bukan gen. Oleh karena itu, anotasi ulang diperlukan menggunakan
#database resmi Bioconductor.

#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  ecoli2.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#PART I.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI 

#Boxplot digunakan untuk:
#- Mengecek distribusi nilai ekspresi antar sampel
#- Melihat apakah ada batch effect
#- Mengevaluasi apakah normalisasi/log-transform sudah wajar

#Set warna berdasarkan grup
group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)


#PART I.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT) 

#Density plot menunjukkan sebaran global nilai ekspresi gen
#Digunakan untuk:
#- Mengecek efek log-transform
#- Membandingkan distribusi antar grup

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )



# PART I.3 UMAP

# 1. Transpose matriks ekspresi
umap_input <- t(ex)

# 2. Hitung jumlah sampel secara otomatis
n_samples <- nrow(umap_input)

# 3. Tentukan parameter n_neighbors secara dinamis
# n_neighbors harus lebih kecil dari jumlah sampel (minimal 2)
custom_config <- umap.defaults
custom_config$n_neighbors <- min(n_samples - 1, 15)

# 4. Jalankan UMAP dengan konfigurasi baru
umap_result <- umap(umap_input, config = custom_config)

# 5. Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

# 6. Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) + # Ukuran titik diperbesar agar jelas
  theme_minimal() +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) + # Konsisten dengan Density Plot
  labs(
    title = "UMAP Plot Sampel",
    subtitle = paste("n_neighbors disesuaikan ke:", custom_config$n_neighbors),
    x = "UMAP 1",
    y = "UMAP 2"
  )

#PART J.1 VISUALISASI VOLCANO PLOT 

#Volcano plot menggabungkan:
#- Log fold change (efek biologis)
#- Signifikansi statistik

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG")


#PART J.2 VISUALISASI HEATMAP 

#Heatmap digunakan untuk melihat pola ekspresi gen
#antar sampel berdasarkan gen-gen paling signifikan

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap 
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)


#PART K. MENYIMPAN HASIL 

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE68106_DEG.csv")

message("Analisis selesai. File hasil telah disimpan.")

# PART L. ENRICHMENT ANALYSIS (GO & KEGG)

library(clusterProfiler)
library(org.EcK12.eg.db)
library(enrichplot)
library(dplyr)

# 2. Persiapan Daftar Gen (Gene List)
# Kita ambil gen yang signifikan (FDR < 0.05) dan memiliki Symbol
deg_list <- topTableResults %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  filter(adj.P.Val < 0.05)

# 3. Mapping ID ke Entrez ID
# KEGG membutuhkan Entrez ID (numerik), bukan Probe ID atau Symbol
gene_convert <- bitr(deg_list$SYMBOL, 
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.EcK12.eg.db)

# Menyiapkan vektor LogFC untuk visualisasi (opsional untuk cnetplot)
gene_with_fc <- merge(gene_convert, deg_list[, c("SYMBOL", "logFC")], by = "SYMBOL")
entrez_vector <- gene_with_fc$logFC
names(entrez_vector) <- gene_with_fc$ENTREZID


# --- L.1 GENE ONTOLOGY (GO) ENRICHMENT ---

message("Menjalankan Analisis GO...")
ego <- enrichGO(gene          = gene_convert$ENTREZID,
                OrgDb         = org.EcK12.eg.db,
                ont           = "BP",         # BP: Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE)         # Mengubah output Entrez kembali ke Symbol

# Visualisasi GO
# Dotplot: Melihat rasio gen dan signifikansi
dotplot(ego, showCategory = 15) + ggtitle("GO Enrichment: Biological Process")


# Barplot: Melihat jumlah gen per kategori
barplot(ego, showCategory = 15)


# --- L.2 KEGG PATHWAY ENRICHMENT ---

message("Menjalankan Analisis KEGG...")
# 'eco' adalah kode organisme KEGG untuk E. coli K-12 MG1655
ekegg <- enrichKEGG(gene         = gene_convert$ENTREZID,
                    organism     = 'eco',
                    pvalueCutoff = 0.05)

# Visualisasi KEGG
if(!is.null(ekegg) && nrow(ekegg) > 0) {
  dotplot(ekegg, showCategory = 15) + ggtitle("KEGG Pathway Enrichment")
} else {
  message("Peringatan: Tidak ditemukan pathway KEGG yang signifikan.")
}


# --- L.3 VISUALISASI JARINGAN (CNETPLOT) ---

# Menghubungkan gen spesifik dengan jalur biologi yang relevan
ego_sim <- setReadable(ego, 'org.EcK12.eg.db', 'ENTREZID') 
cnetplot(ego_sim, 
         foldChange = entrez_vector, 
         circular = TRUE, 
         colorEdge = TRUE) + 
  ggtitle("Gene-Concept Network (GO)")



# --- PART M. MENYIMPAN HASIL ---

write.csv(as.data.frame(ego), "Hasil_GO_Enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "Hasil_KEGG_Enrichment.csv", row.names = FALSE)

message("Analisis selesai! Cek file CSV dan jendela Plot Anda.")




