PENDAHULUAN<br>
Studi ini bertujuan untuk menyelidiki respons transkriptomik dari Escherichia coli CFT073 terhadap perlakuan antibiotik ciprofloxacin dalam model blood stream infection.<br>
Penelitian dilakukan dengan menumbuhkan bakteri pada pooled human blood dan memberikan tantangan antibiotik pada fase pertumbuhan mid-logaritmik untuk mensimulasikan skenario klinis. Profil ekspresi gen global kemudian dianalisis menggunakan microarray DNA untuk membandingkan respons antara kelompok yang diberi ciprofloxacin dengan kelompok kontrol.<br>
METODE<br>
Dataset yang digunakan adalah GSE68106. Sampel dikelompokkan menjadi Cipro_Treatment dan Control_30h. Data ekspresi memiliki rentang nilai yang sudah kecil dan terdistribusi secara normal, sehingga tidak memerlukan transformasi log2. Fitur gen dianotasi menggunakan ecoli2.db.<br>
Differential expression dianalisis menggunakan limma dengan metode eBayes dan FDR < 0.01. Data divisualisasi menggunakan box dan density plot untuk verifikasi distribusi data, UMAP untuk melihat clustering sampel, dan volcano plot untuk menampilkan hubungan antara log fold change dan signifikansi statistik. Gene onthology dan KEGG pathway digunakan untuk analisis pengayaan fungsional.<br>
HASIL DAN INTERPRETASI<br>
Boxplot menunjukkan data ternormalisasi dengan baik antar sampel. Ini diperkuat dengan overlapping distribusi nilai ekspresi gen antar perlakuan.<br>
Beberapa sampel (terutama pada grup kontrol) terlihat sedikit menyebar, namun secara umum kedua grup membentuk cluster yang berbeda pada UMAP.<br>
Analisis differentially expressed gene menunjukkan 48 gen terupregulasi dan 2 gen terdownregulasi pada kelompok Cipro_Treatment. Gen-gen terupregulasi tersebut berfungsi dalam cell cycle arrest (sulA, recA, dan lexA), DNA repair (umuC, recN, dan cho), serta stress response (rzpD dan yebF), sedangkan downregulasi cydA terkait penurunan respirasi aerobik.<br>
GO enrichment menunjukkan fungsi biologis yang meningkat pada treatment ciprofloxacin yaitu response to external stimuli, DNA metabolic process dan damage response, serta cellular stress response. Biological function ini berkaitan dengan gen-gen yang terupregulasi secara signifikan di volcano plot dan heatmap. Gene concept network menunjukkan sulA dan recN memiliki fold change tertinggi dan berperan penting dalam respon stress akibat treatment ciprofloxacin.<br>
Gen-gen yang terpengaruh oleh ciprofloxacin dalam eksperimen ini tidak terkumpul dalam sebuah KEGG pathway, melainkan lebih tersebar di berbagai fungsi biologis yang lebih luas. Sehingga enrichment KEGG belum dapat dilakukan.<br>
KESIMPULAN<br>
Treatment ciprofloxacin menyebabkan respon stress pada E. coli CFT073. Mekanisme DNA damage dan SOS response utamanya meningkat oleh regulasi gen sulA dan recN.
