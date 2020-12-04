R Notebook
================

  - [Partie dada2 :](#partie-dada2)
  - [Inspecter les profils de qualité de
    lecture](#inspecter-les-profils-de-qualité-de-lecture)
  - [Filtrer et couper :](#filtrer-et-couper)
  - [Apprendre les taux d’erreurs :](#apprendre-les-taux-derreurs)
      - [Création d’un modèle d’erreur paramétrique
        :](#création-dun-modèle-derreur-paramétrique)
  - [Inférence d’échantillon](#inférence-déchantillon)
  - [Fusionner les paires de reads :](#fusionner-les-paires-de-reads)
  - [Construire une table de séquence
    :](#construire-une-table-de-séquence)
  - [Supprimer les chimères](#supprimer-les-chimères)
  - [Suivre les lectures dans le pipeline
    :](#suivre-les-lectures-dans-le-pipeline)
  - [Attribuer une taxonomie :](#attribuer-une-taxonomie)
      - [On va prendre la base Silva pour lui attribuer une taxonomie
        :](#on-va-prendre-la-base-silva-pour-lui-attribuer-une-taxonomie)
  - [Evaluer la précision :](#evaluer-la-précision)
      - [Evalution de la précision de DADA2 sur la communauté estimée
        :](#evalution-de-la-précision-de-dada2-sur-la-communauté-estimée)
  - [Bonus : Transfert à phyloseq :](#bonus-transfert-à-phyloseq)
      - [Visualisez l’alpha-diversité.](#visualisez-lalpha-diversité.)
      - [Création d’un graphique d’ordination
        :](#création-dun-graphique-dordination)
  - [Histogrammes :](#histogrammes)
  - [Tutoriel phyloseq](#tutoriel-phyloseq)
  - [Import des données phyloseq du tutoriel phyloseq
    :](#import-des-données-phyloseq-du-tutoriel-phyloseq)
  - [Filtration taxonomique :](#filtration-taxonomique)
  - [Filtration de prévalence :](#filtration-de-prévalence)
  - [Taxons agglomérés :](#taxons-agglomérés)
  - [Transformation de la valeur d’abondance
    :](#transformation-de-la-valeur-dabondance)
  - [Sous-ensemble par taxonomie :](#sous-ensemble-par-taxonomie)
  - [Prétraitement :](#prétraitement)
  - [Différents projections d’ordination
    :](#différents-projections-dordination)
  - [PCA sur les rangs :](#pca-sur-les-rangs)
  - [Analyse de correspondance canonique
    :](#analyse-de-correspondance-canonique)
  - [Enseignement supervisé :](#enseignement-supervisé)
  - [Analyses basées sur des graphiques
    :](#analyses-basées-sur-des-graphiques)
      - [Créer et tracer des graphiques
        :](#créer-et-tracer-des-graphiques)
      - [Tests à deux échantillons basée sur des graphiques
        :](#tests-à-deux-échantillons-basée-sur-des-graphiques)
      - [Minimum Spanning Tree (MST) :](#minimum-spanning-tree-mst)
      - [Voisins les plus proches :](#voisins-les-plus-proches)
      - [Modélisation linéaire :](#modélisation-linéaire)
      - [Tests multiples hierarchiques
        :](#tests-multiples-hierarchiques)
  - [Multitable techniques :](#multitable-techniques)

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk
or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.

# Partie dada2 :

Chargement du package dada2 :

``` r
library("dada2")
```

    ## Loading required package: Rcpp

Création d’un chemin vers le répertoire MiSeq SOP contenant les données
fastq :

``` r
path <- "~/MiSeq_SOP"
list.files(path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

Lire des fichiers fastq et affecter les séquences Forward aà la variable
fnFs et les séquences Reverse à la variable fnRs : Comment cela
fonctionne ? La variable fnFs est une nouvelle variable. “path” est
relié dans le dossier MiSeq SOP. Le “pattern” est une expression
régulière qui fait qu’ici tout ce qui est avec "\_R1\_001.fastq" va
être chercher et ranger dans le fichier. Grâce à la fonction “sort”
elle range par ordre alphabétique. Ceci fonctionne de la même façon avec
fnRs qui affecte tous les "\_R2\_OO1.fastq“. De plus, grâce à la
fonction”strsplit" on peut isoler le nom des séquences et ne garder que
l’identifiant. “sapply” permet de réaliser ceci sur toutes les
séquences. Puis on va renvoyer tout ceci à la variable “sample.names”.

``` r
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
```

# Inspecter les profils de qualité de lecture

Visualiser les profils de qualité des lectures des reads que l’on a
rangé dans la variable fnFs. La fonction “plot” permet de faire un
graphique et le chiffre 4 renvoie aux quatre premiers reads de la
variable fnFs.

``` r
plotQualityProfile(fnFs[1:4])
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Analyse des figures : l’axe des ordonnées permet de renvoyer les scores
de qualités. l’axe des abscisses donne le nombre de nucléotides : ici
avec Illumina il y a 250 pb. La ligne rouge indique la proportion de
reads qui sont à cette position. La couleur grise renvoie aux heat map
ce qui permet de déterminer la fréquence. Lorsque les couleurs sont
foncées cela renvoie à une fréquence plus élevée. La ligne verte
représente la moyenne. La ligne orange est la médiane et lorsqu’elle
est en pointillée elle indique les 25ème et 75ème quartile. On observe
globalement un bon score de qualité au dessus d’un score de 30. En
revanche on peut tout de même noter une diminution du score de qualité
pour les 10 derniers nucléotides. Il peut donc être nécessaire de
tronquer les 10 dernières positions.

# Filtrer et couper :

Ici on doit filtrer les données Forward et Reverse. Pour cela on doit
constuire un chemin d’accès à un fichier ici : filtFs er filtRs. Ceci
est permis par la fonction “file.path”. “filtered” va permettre de
filtrer.

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

On définit une variable out avec nos données Reverse et Forward, et
Reverse et Forward filtrés. “truncLen” va permettre de tronquer les
reads entre les positions 240 et 160. On fait ceci car dans les scores
de qualités des reads des Forward on avait un mauvais score de qualité,
et pour les reads des Reverse il s’agissait des 90 derniers. Puis on va
utiliser des paramètres de filtration standard avec maxEE qui renvoie au
nombre maximum d’erreur attendue. Puis grâce à la fonction “head” on
peut montrer la tête du fichier out.

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

# Apprendre les taux d’erreurs :

## Création d’un modèle d’erreur paramétrique :

Dada2 va créer un modèle grâce à la fonction “learnErrors” qui va
alterner entre estimation des taux d’erruers et inférence de la
compositon jusqu’à trouver une solution cohérente. Ici errF représente
les reads Forward.

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

Ici c’est la même chose pour les reads Reverse.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

Ici avec la fonction “plotErrors” on visualise les taux d’erreurs
estimés.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Analyse des figures : Sur l’axe des ordonnées on peut voir la fréquence
des erreurs en log10. Sur l’axe des abscisses on a le score de qualité.
Les taux d’erreurs pour chaque transition possible sont indiqués. La
ligne noire indique les taux d’erreur que l’on a estimé. La ligne rouge
indique les taux d’erreur attendus. En diagonale cela indique la
probabilité que A réalise une transition en A, C en C, G en G… Plus le
score de qualité est haute plus la fréquence d’erreur est basse. C’est
bien ce que l’on observe ici un diminution des taux observés avec une
augmentation de la qualité. De plus, les taux d’erreur estimés sont en
corrélation avec les taux observés.

# Inférence d’échantillon

Ici la fonction dada utilise un test statistique : une séquence qui a
été vue trop de fois ne sera pas jugée comme si elle a été causé par
des erreurs d’amplicon. Ainsi, les séquences abondantes sont gardées et
les séquences qui ne le sont pas sont enlevées à partir du seuil que
l’on a déterminé précédemment. Ici il s’agissait des reads Forward.

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

Ici, c’est la même fonction pour les reads Reverse.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

Ici dadaFs renvoie l’objet dada class qui va déduire 128 vraies
variantes de séquence à partir des 1979 uniques, et ce dans le premier
échantillon.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 128 sequence variants were inferred from 1979 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Fusionner les paires de reads :

Ici on doit fusionner les reads Forward et Reverse. La fonction
“mergePairs” va permettre de rassembler par paire pour former des
contigs. Tout ceci va se faire dans l’objet mergers. Puis avec “head” on
montre le début.

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 6551 paired-reads (in 106 unique pairings) successfully merged out of 6907 (in 199 pairings) input.

    ## 5025 paired-reads (in 100 unique pairings) successfully merged out of 5188 (in 156 pairings) input.

    ## 4973 paired-reads (in 80 unique pairings) successfully merged out of 5268 (in 166 pairings) input.

    ## 2595 paired-reads (in 52 unique pairings) successfully merged out of 2756 (in 109 pairings) input.

    ## 2553 paired-reads (in 60 unique pairings) successfully merged out of 2785 (in 119 pairings) input.

    ## 3622 paired-reads (in 53 unique pairings) successfully merged out of 4103 (in 157 pairings) input.

    ## 6079 paired-reads (in 81 unique pairings) successfully merged out of 6515 (in 198 pairings) input.

    ## 3961 paired-reads (in 90 unique pairings) successfully merged out of 4384 (in 188 pairings) input.

    ## 14231 paired-reads (in 143 unique pairings) successfully merged out of 15358 (in 351 pairings) input.

    ## 10526 paired-reads (in 120 unique pairings) successfully merged out of 11166 (in 279 pairings) input.

    ## 11156 paired-reads (in 137 unique pairings) successfully merged out of 11799 (in 298 pairings) input.

    ## 4329 paired-reads (in 84 unique pairings) successfully merged out of 4788 (in 180 pairings) input.

    ## 17431 paired-reads (in 153 unique pairings) successfully merged out of 17812 (in 272 pairings) input.

    ## 5850 paired-reads (in 81 unique pairings) successfully merged out of 6095 (in 159 pairings) input.

    ## 3716 paired-reads (in 86 unique pairings) successfully merged out of 3894 (in 147 pairings) input.

    ## 6865 paired-reads (in 99 unique pairings) successfully merged out of 7193 (in 187 pairings) input.

    ## 4430 paired-reads (in 67 unique pairings) successfully merged out of 4605 (in 127 pairings) input.

    ## 4574 paired-reads (in 100 unique pairings) successfully merged out of 4736 (in 172 pairings) input.

    ## 6094 paired-reads (in 109 unique pairings) successfully merged out of 6314 (in 172 pairings) input.

    ## 4269 paired-reads (in 20 unique pairings) successfully merged out of 4281 (in 28 pairings) input.

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                       sequence
    ## 1 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACAGG
    ## 2 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACAGG
    ## 3 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATGTCGGGGCTCAACCCCGGCCTGCCGTTGAAACTGGCGGCCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ## 4 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTTTTAAGTCAGCGGTAAAAATTCGGGGCTCAACCCCGTCCGGCCGTTGAAACTGGGGGCCTTGAGTGGGCGAGAAGAAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCCTTCCGGCGCCCTACTGACGCTGAGGCACGAAAGTGCGGGGATCGAACAGG
    ## 5 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGACTCTCAAGTCAGCGGTCAAATCGCGGGGCTCAACCCCGTTCCGCCGTTGAAACTGGGAGCCTTGAGTGCGCGAGAAGTAGGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCCTACCGGCGCGCAACTGACGCTCATGCACGAAAGCGTGGGTATCGAACAGG
    ## 6 TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAAGCGGTGCTCAACGCCGTCGAGCCGTTGAAACTGGCGTTCTTGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAGGCACGAAAGCGTGGGTATCGAACAGG
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1       579       1       1    148         0      0      1   TRUE
    ## 2       470       2       2    148         0      0      2   TRUE
    ## 3       449       3       4    148         0      0      1   TRUE
    ## 4       430       4       3    148         0      0      2   TRUE
    ## 5       345       5       6    148         0      0      1   TRUE
    ## 6       282       6       5    148         0      0      2   TRUE

# Construire une table de séquence :

La fonction “makeSequenceTable” va construire une table de séquence à
partir de l’échantillon que l’on a fusionné avant, contenant les Forward
et les Reverse. Cette table de séquence appelé amplicon sequence variant
table permet une plus meilleure résolution de la table OTU. dim permet,
quant à lui, la dimense de notre table de séquence.

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]  20 293

Ici on inspecte la distribution de la longueur de la séquence, c’est à
dire on regarde le nombre de caractère et la distribution. On a une
séquence qui a 251 pb, 88 qui en ont 252 pb.

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  88 196   6   2

# Supprimer les chimères

La fonction “removeBimeraDenovo” va permettre d’éliminer les chimères.
Cela permet de renvoyer des variants de séquences sans chimères dans
seqtab.nochim.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 61 bimeras out of 293 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]  20 232

Ici on fait un rapport de la somme des séquences après avoir enlevé les
séquences chimériques sur les séquences avant cette supression. Cela
permet de voir que les séquences chimériques représentent moins de 4%.

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.964263

# Suivre les lectures dans le pipeline :

Ici le code examine le nombre de lectures effectuées à chaque étape. La
fonction “cbind” permet de “coller” des tableaux, c’est à dire qu’elle
va regrouper les colonnes. De plus elle va extraire les séquences
uniques de l’objet dada. Avec “head” on affiche le début du tableau. On
peut voir que la majorité a été conservé ce qui montre a nouveau une
cohérence des données.

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##        input filtered denoisedF denoisedR merged nonchim
    ## F3D0    7793     7113      6996      6978   6551    6539
    ## F3D1    5869     5299      5227      5239   5025    5014
    ## F3D141  5958     5463      5339      5351   4973    4850
    ## F3D142  3183     2914      2799      2833   2595    2521
    ## F3D143  3178     2941      2822      2868   2553    2519
    ## F3D144  4827     4312      4146      4224   3622    3483

Dans le tableau input indique le nombre de séquences. Filtered revoie a
ce qui a été filtré avec un score de qualité. denoisedF : le nombre de
séquences après l’analyse de dada2 pour les séquences Forward.
denoisedR : idem avec les séquences Reverse. mergerd indique le nombres
de séquences qui forment des contigs. Enfin nonchim indique le nombre de
séquences sans les chimères.

# Attribuer une taxonomie :

## On va prendre la base Silva pour lui attribuer une taxonomie :

“wget” permet de télécharger un fichier sur Internet. Attention pour
cette étape il faut être en bash.

``` bash
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
```

    ## --2020-12-04 14:23:00--  https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137973851 (132M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138_train_set.fa.gz.8’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 4.02M 33s
    ##     50K .......... .......... .......... .......... ..........  0% 9.10M 24s
    ##    100K .......... .......... .......... .......... ..........  0% 13.6M 19s
    ##    150K .......... .......... .......... .......... ..........  0% 46.7M 15s
    ##    200K .......... .......... .......... .......... ..........  0% 16.9M 13s
    ##    250K .......... .......... .......... .......... ..........  0% 49.8M 12s
    ##    300K .......... .......... .......... .......... ..........  0% 14.8M 11s
    ##    350K .......... .......... .......... .......... ..........  0% 82.8M 10s
    ##    400K .......... .......... .......... .......... ..........  0% 18.2M 10s
    ##    450K .......... .......... .......... .......... ..........  0% 72.7M 9s
    ##    500K .......... .......... .......... .......... ..........  0% 12.9M 9s
    ##    550K .......... .......... .......... .......... ..........  0% 28.6M 9s
    ##    600K .......... .......... .......... .......... ..........  0% 24.0M 8s
    ##    650K .......... .......... .......... .......... ..........  0% 32.5M 8s
    ##    700K .......... .......... .......... .......... ..........  0% 22.6M 8s
    ##    750K .......... .......... .......... .......... ..........  0% 47.1M 8s
    ##    800K .......... .......... .......... .......... ..........  0% 19.9M 8s
    ##    850K .......... .......... .......... .......... ..........  0% 36.4M 7s
    ##    900K .......... .......... .......... .......... ..........  0% 23.3M 7s
    ##    950K .......... .......... .......... .......... ..........  0% 37.2M 7s
    ##   1000K .......... .......... .......... .......... ..........  0% 21.8M 7s
    ##   1050K .......... .......... .......... .......... ..........  0% 38.3M 7s
    ##   1100K .......... .......... .......... .......... ..........  0% 19.8M 7s
    ##   1150K .......... .......... .......... .......... ..........  0% 15.5M 7s
    ##   1200K .......... .......... .......... .......... ..........  0% 65.2M 7s
    ##   1250K .......... .......... .......... .......... ..........  0% 67.2M 6s
    ##   1300K .......... .......... .......... .......... ..........  1% 16.1M 7s
    ##   1350K .......... .......... .......... .......... ..........  1% 12.4M 7s
    ##   1400K .......... .......... .......... .......... ..........  1% 66.9M 7s
    ##   1450K .......... .......... .......... .......... ..........  1% 14.1M 7s
    ##   1500K .......... .......... .......... .......... ..........  1% 78.3M 6s
    ##   1550K .......... .......... .......... .......... ..........  1% 12.8M 7s
    ##   1600K .......... .......... .......... .......... ..........  1% 78.7M 6s
    ##   1650K .......... .......... .......... .......... ..........  1% 88.8M 6s
    ##   1700K .......... .......... .......... .......... ..........  1% 14.1M 6s
    ##   1750K .......... .......... .......... .......... ..........  1% 98.0M 6s
    ##   1800K .......... .......... .......... .......... ..........  1% 19.5M 6s
    ##   1850K .......... .......... .......... .......... ..........  1% 72.1M 6s
    ##   1900K .......... .......... .......... .......... ..........  1% 20.6M 6s
    ##   1950K .......... .......... .......... .......... ..........  1% 62.8M 6s
    ##   2000K .......... .......... .......... .......... ..........  1% 74.4M 6s
    ##   2050K .......... .......... .......... .......... ..........  1% 18.5M 6s
    ##   2100K .......... .......... .......... .......... ..........  1% 84.8M 6s
    ##   2150K .......... .......... .......... .......... ..........  1% 22.4M 6s
    ##   2200K .......... .......... .......... .......... ..........  1% 27.8M 6s
    ##   2250K .......... .......... .......... .......... ..........  1% 85.5M 6s
    ##   2300K .......... .......... .......... .......... ..........  1% 16.9M 6s
    ##   2350K .......... .......... .......... .......... ..........  1% 70.8M 6s
    ##   2400K .......... .......... .......... .......... ..........  1% 77.7M 6s
    ##   2450K .......... .......... .......... .......... ..........  1% 14.9M 6s
    ##   2500K .......... .......... .......... .......... ..........  1% 80.8M 6s
    ##   2550K .......... .......... .......... .......... ..........  1% 15.8M 6s
    ##   2600K .......... .......... .......... .......... ..........  1% 60.0M 6s
    ##   2650K .......... .......... .......... .......... ..........  2%  101M 5s
    ##   2700K .......... .......... .......... .......... ..........  2% 20.7M 5s
    ##   2750K .......... .......... .......... .......... ..........  2% 71.2M 5s
    ##   2800K .......... .......... .......... .......... ..........  2% 23.7M 5s
    ##   2850K .......... .......... .......... .......... ..........  2% 41.4M 5s
    ##   2900K .......... .......... .......... .......... ..........  2% 70.7M 5s
    ##   2950K .......... .......... .......... .......... ..........  2% 25.3M 5s
    ##   3000K .......... .......... .......... .......... ..........  2% 28.0M 5s
    ##   3050K .......... .......... .......... .......... ..........  2% 70.6M 5s
    ##   3100K .......... .......... .......... .......... ..........  2% 24.2M 5s
    ##   3150K .......... .......... .......... .......... ..........  2% 35.2M 5s
    ##   3200K .......... .......... .......... .......... ..........  2% 58.9M 5s
    ##   3250K .......... .......... .......... .......... ..........  2% 18.8M 5s
    ##   3300K .......... .......... .......... .......... ..........  2% 83.7M 5s
    ##   3350K .......... .......... .......... .......... ..........  2% 85.2M 5s
    ##   3400K .......... .......... .......... .......... ..........  2% 14.5M 5s
    ##   3450K .......... .......... .......... .......... ..........  2% 80.3M 5s
    ##   3500K .......... .......... .......... .......... ..........  2% 15.5M 5s
    ##   3550K .......... .......... .......... .......... ..........  2% 81.2M 5s
    ##   3600K .......... .......... .......... .......... ..........  2% 79.0M 5s
    ##   3650K .......... .......... .......... .......... ..........  2% 13.5M 5s
    ##   3700K .......... .......... .......... .......... ..........  2% 82.0M 5s
    ##   3750K .......... .......... .......... .......... ..........  2% 48.4M 5s
    ##   3800K .......... .......... .......... .......... ..........  2% 30.3M 5s
    ##   3850K .......... .......... .......... .......... ..........  2% 74.8M 5s
    ##   3900K .......... .......... .......... .......... ..........  2% 29.6M 5s
    ##   3950K .......... .......... .......... .......... ..........  2% 52.0M 5s
    ##   4000K .......... .......... .......... .......... ..........  3% 17.4M 5s
    ##   4050K .......... .......... .......... .......... ..........  3% 92.3M 5s
    ##   4100K .......... .......... .......... .......... ..........  3% 90.9M 5s
    ##   4150K .......... .......... .......... .......... ..........  3% 32.8M 5s
    ##   4200K .......... .......... .......... .......... ..........  3% 32.9M 5s
    ##   4250K .......... .......... .......... .......... ..........  3% 80.9M 5s
    ##   4300K .......... .......... .......... .......... ..........  3% 23.8M 5s
    ##   4350K .......... .......... .......... .......... ..........  3% 30.4M 5s
    ##   4400K .......... .......... .......... .......... ..........  3% 30.3M 5s
    ##   4450K .......... .......... .......... .......... ..........  3%  120M 5s
    ##   4500K .......... .......... .......... .......... ..........  3% 23.7M 5s
    ##   4550K .......... .......... .......... .......... ..........  3% 27.2M 5s
    ##   4600K .......... .......... .......... .......... ..........  3% 21.6M 5s
    ##   4650K .......... .......... .......... .......... ..........  3% 33.3M 5s
    ##   4700K .......... .......... .......... .......... ..........  3% 81.5M 5s
    ##   4750K .......... .......... .......... .......... ..........  3% 26.3M 5s
    ##   4800K .......... .......... .......... .......... ..........  3% 30.4M 5s
    ##   4850K .......... .......... .......... .......... ..........  3% 9.14M 5s
    ##   4900K .......... .......... .......... .......... ..........  3% 72.2M 5s
    ##   4950K .......... .......... .......... .......... ..........  3% 85.8M 5s
    ##   5000K .......... .......... .......... .......... ..........  3% 94.8M 5s
    ##   5050K .......... .......... .......... .......... ..........  3% 14.8M 5s
    ##   5100K .......... .......... .......... .......... ..........  3% 93.6M 5s
    ##   5150K .......... .......... .......... .......... ..........  3% 31.4M 5s
    ##   5200K .......... .......... .......... .......... ..........  3% 78.4M 5s
    ##   5250K .......... .......... .......... .......... ..........  3% 57.3M 5s
    ##   5300K .......... .......... .......... .......... ..........  3% 60.9M 5s
    ##   5350K .......... .......... .......... .......... ..........  4% 29.9M 5s
    ##   5400K .......... .......... .......... .......... ..........  4% 37.5M 5s
    ##   5450K .......... .......... .......... .......... ..........  4% 67.2M 5s
    ##   5500K .......... .......... .......... .......... ..........  4% 72.6M 4s
    ##   5550K .......... .......... .......... .......... ..........  4% 35.4M 4s
    ##   5600K .......... .......... .......... .......... ..........  4% 60.2M 4s
    ##   5650K .......... .......... .......... .......... ..........  4% 82.2M 4s
    ##   5700K .......... .......... .......... .......... ..........  4% 72.9M 4s
    ##   5750K .......... .......... .......... .......... ..........  4% 28.3M 4s
    ##   5800K .......... .......... .......... .......... ..........  4% 34.0M 4s
    ##   5850K .......... .......... .......... .......... ..........  4%  124M 4s
    ##   5900K .......... .......... .......... .......... ..........  4%  102M 4s
    ##   5950K .......... .......... .......... .......... ..........  4% 23.2M 4s
    ##   6000K .......... .......... .......... .......... ..........  4% 68.3M 4s
    ##   6050K .......... .......... .......... .......... ..........  4% 37.3M 4s
    ##   6100K .......... .......... .......... .......... ..........  4% 32.7M 4s
    ##   6150K .......... .......... .......... .......... ..........  4% 60.3M 4s
    ##   6200K .......... .......... .......... .......... ..........  4%  102M 4s
    ##   6250K .......... .......... .......... .......... ..........  4% 39.2M 4s
    ##   6300K .......... .......... .......... .......... ..........  4% 44.7M 4s
    ##   6350K .......... .......... .......... .......... ..........  4% 62.4M 4s
    ##   6400K .......... .......... .......... .......... ..........  4% 68.7M 4s
    ##   6450K .......... .......... .......... .......... ..........  4% 46.2M 4s
    ##   6500K .......... .......... .......... .......... ..........  4% 58.5M 4s
    ##   6550K .......... .......... .......... .......... ..........  4% 24.8M 4s
    ##   6600K .......... .......... .......... .......... ..........  4%  104M 4s
    ##   6650K .......... .......... .......... .......... ..........  4% 64.5M 4s
    ##   6700K .......... .......... .......... .......... ..........  5% 47.9M 4s
    ##   6750K .......... .......... .......... .......... ..........  5% 5.85M 4s
    ##   6800K .......... .......... .......... .......... ..........  5% 98.0M 4s
    ##   6850K .......... .......... .......... .......... ..........  5%  102M 4s
    ##   6900K .......... .......... .......... .......... ..........  5%  110M 4s
    ##   6950K .......... .......... .......... .......... ..........  5% 20.5M 4s
    ##   7000K .......... .......... .......... .......... ..........  5% 51.5M 4s
    ##   7050K .......... .......... .......... .......... ..........  5% 79.1M 4s
    ##   7100K .......... .......... .......... .......... ..........  5%  116M 4s
    ##   7150K .......... .......... .......... .......... ..........  5% 26.0M 4s
    ##   7200K .......... .......... .......... .......... ..........  5% 22.8M 4s
    ##   7250K .......... .......... .......... .......... ..........  5% 69.2M 4s
    ##   7300K .......... .......... .......... .......... ..........  5%  102M 4s
    ##   7350K .......... .......... .......... .......... ..........  5% 67.7M 4s
    ##   7400K .......... .......... .......... .......... ..........  5% 18.4M 4s
    ##   7450K .......... .......... .......... .......... ..........  5% 87.8M 4s
    ##   7500K .......... .......... .......... .......... ..........  5% 64.7M 4s
    ##   7550K .......... .......... .......... .......... ..........  5% 89.9M 4s
    ##   7600K .......... .......... .......... .......... ..........  5% 36.8M 4s
    ##   7650K .......... .......... .......... .......... ..........  5% 38.6M 4s
    ##   7700K .......... .......... .......... .......... ..........  5% 97.2M 4s
    ##   7750K .......... .......... .......... .......... ..........  5% 28.4M 4s
    ##   7800K .......... .......... .......... .......... ..........  5% 90.8M 4s
    ##   7850K .......... .......... .......... .......... ..........  5% 33.4M 4s
    ##   7900K .......... .......... .......... .......... ..........  5% 59.5M 4s
    ##   7950K .......... .......... .......... .......... ..........  5% 87.6M 4s
    ##   8000K .......... .......... .......... .......... ..........  5% 20.7M 4s
    ##   8050K .......... .......... .......... .......... ..........  6%  114M 4s
    ##   8100K .......... .......... .......... .......... ..........  6% 78.8M 4s
    ##   8150K .......... .......... .......... .......... ..........  6% 36.0M 4s
    ##   8200K .......... .......... .......... .......... ..........  6% 24.9M 4s
    ##   8250K .......... .......... .......... .......... ..........  6% 33.7M 4s
    ##   8300K .......... .......... .......... .......... ..........  6%  106M 4s
    ##   8350K .......... .......... .......... .......... ..........  6%  103M 4s
    ##   8400K .......... .......... .......... .......... ..........  6% 16.4M 4s
    ##   8450K .......... .......... .......... .......... ..........  6%  114M 4s
    ##   8500K .......... .......... .......... .......... ..........  6% 96.0M 4s
    ##   8550K .......... .......... .......... .......... ..........  6%  117M 4s
    ##   8600K .......... .......... .......... .......... ..........  6% 19.6M 4s
    ##   8650K .......... .......... .......... .......... ..........  6%  108M 4s
    ##   8700K .......... .......... .......... .......... ..........  6% 84.5M 4s
    ##   8750K .......... .......... .......... .......... ..........  6%  125M 4s
    ##   8800K .......... .......... .......... .......... ..........  6% 22.2M 4s
    ##   8850K .......... .......... .......... .......... ..........  6% 41.1M 4s
    ##   8900K .......... .......... .......... .......... ..........  6%  107M 4s
    ##   8950K .......... .......... .......... .......... ..........  6% 20.7M 4s
    ##   9000K .......... .......... .......... .......... ..........  6% 78.2M 4s
    ##   9050K .......... .......... .......... .......... ..........  6% 61.7M 4s
    ##   9100K .......... .......... .......... .......... ..........  6%  117M 4s
    ##   9150K .......... .......... .......... .......... ..........  6% 28.1M 4s
    ##   9200K .......... .......... .......... .......... ..........  6% 81.6M 4s
    ##   9250K .......... .......... .......... .......... ..........  6% 46.7M 4s
    ##   9300K .......... .......... .......... .......... ..........  6%  102M 4s
    ##   9350K .......... .......... .......... .......... ..........  6% 49.0M 4s
    ##   9400K .......... .......... .......... .......... ..........  7% 29.7M 4s
    ##   9450K .......... .......... .......... .......... ..........  7% 34.8M 4s
    ##   9500K .......... .......... .......... .......... ..........  7%  108M 4s
    ##   9550K .......... .......... .......... .......... ..........  7% 37.9M 4s
    ##   9600K .......... .......... .......... .......... ..........  7% 24.5M 4s
    ##   9650K .......... .......... .......... .......... ..........  7% 75.8M 4s
    ##   9700K .......... .......... .......... .......... ..........  7%  122M 4s
    ##   9750K .......... .......... .......... .......... ..........  7% 43.7M 4s
    ##   9800K .......... .......... .......... .......... ..........  7% 22.7M 4s
    ##   9850K .......... .......... .......... .......... ..........  7% 94.4M 4s
    ##   9900K .......... .......... .......... .......... ..........  7%  117M 4s
    ##   9950K .......... .......... .......... .......... ..........  7% 62.6M 4s
    ##  10000K .......... .......... .......... .......... ..........  7% 29.7M 4s
    ##  10050K .......... .......... .......... .......... ..........  7% 55.3M 4s
    ##  10100K .......... .......... .......... .......... ..........  7% 73.0M 4s
    ##  10150K .......... .......... .......... .......... ..........  7% 99.9M 4s
    ##  10200K .......... .......... .......... .......... ..........  7% 25.3M 4s
    ##  10250K .......... .......... .......... .......... ..........  7% 93.1M 4s
    ##  10300K .......... .......... .......... .......... ..........  7% 60.6M 4s
    ##  10350K .......... .......... .......... .......... ..........  7% 24.8M 4s
    ##  10400K .......... .......... .......... .......... ..........  7% 64.8M 4s
    ##  10450K .......... .......... .......... .......... ..........  7%  120M 4s
    ##  10500K .......... .......... .......... .......... ..........  7% 46.1M 4s
    ##  10550K .......... .......... .......... .......... ..........  7% 29.8M 4s
    ##  10600K .......... .......... .......... .......... ..........  7% 66.5M 4s
    ##  10650K .......... .......... .......... .......... ..........  7%  108M 4s
    ##  10700K .......... .......... .......... .......... ..........  7% 31.0M 4s
    ##  10750K .......... .......... .......... .......... ..........  8% 69.9M 4s
    ##  10800K .......... .......... .......... .......... ..........  8% 73.3M 4s
    ##  10850K .......... .......... .......... .......... ..........  8%  106M 4s
    ##  10900K .......... .......... .......... .......... ..........  8% 28.9M 4s
    ##  10950K .......... .......... .......... .......... ..........  8% 48.4M 4s
    ##  11000K .......... .......... .......... .......... ..........  8% 78.4M 4s
    ##  11050K .......... .......... .......... .......... ..........  8% 48.6M 4s
    ##  11100K .......... .......... .......... .......... ..........  8% 82.9M 3s
    ##  11150K .......... .......... .......... .......... ..........  8% 29.9M 3s
    ##  11200K .......... .......... .......... .......... ..........  8% 94.8M 3s
    ##  11250K .......... .......... .......... .......... ..........  8% 32.2M 3s
    ##  11300K .......... .......... .......... .......... ..........  8% 48.0M 3s
    ##  11350K .......... .......... .......... .......... ..........  8% 50.5M 3s
    ##  11400K .......... .......... .......... .......... ..........  8% 69.6M 3s
    ##  11450K .......... .......... .......... .......... ..........  8% 39.8M 3s
    ##  11500K .......... .......... .......... .......... ..........  8% 38.8M 3s
    ##  11550K .......... .......... .......... .......... ..........  8% 68.1M 3s
    ##  11600K .......... .......... .......... .......... ..........  8% 49.4M 3s
    ##  11650K .......... .......... .......... .......... ..........  8% 58.6M 3s
    ##  11700K .......... .......... .......... .......... ..........  8% 33.4M 3s
    ##  11750K .......... .......... .......... .......... ..........  8% 93.9M 3s
    ##  11800K .......... .......... .......... .......... ..........  8% 44.3M 3s
    ##  11850K .......... .......... .......... .......... ..........  8% 44.8M 3s
    ##  11900K .......... .......... .......... .......... ..........  8% 51.0M 3s
    ##  11950K .......... .......... .......... .......... ..........  8% 55.9M 3s
    ##  12000K .......... .......... .......... .......... ..........  8% 50.9M 3s
    ##  12050K .......... .......... .......... .......... ..........  8%  110M 3s
    ##  12100K .......... .......... .......... .......... ..........  9% 28.7M 3s
    ##  12150K .......... .......... .......... .......... ..........  9% 73.3M 3s
    ##  12200K .......... .......... .......... .......... ..........  9% 41.3M 3s
    ##  12250K .......... .......... .......... .......... ..........  9% 28.3M 3s
    ##  12300K .......... .......... .......... .......... ..........  9% 98.9M 3s
    ##  12350K .......... .......... .......... .......... ..........  9% 47.0M 3s
    ##  12400K .......... .......... .......... .......... ..........  9% 54.9M 3s
    ##  12450K .......... .......... .......... .......... ..........  9% 44.1M 3s
    ##  12500K .......... .......... .......... .......... ..........  9%  108M 3s
    ##  12550K .......... .......... .......... .......... ..........  9% 26.3M 3s
    ##  12600K .......... .......... .......... .......... ..........  9%  121M 3s
    ##  12650K .......... .......... .......... .......... ..........  9% 30.0M 3s
    ##  12700K .......... .......... .......... .......... ..........  9% 45.2M 3s
    ##  12750K .......... .......... .......... .......... ..........  9% 80.4M 3s
    ##  12800K .......... .......... .......... .......... ..........  9% 52.3M 3s
    ##  12850K .......... .......... .......... .......... ..........  9% 44.7M 3s
    ##  12900K .......... .......... .......... .......... ..........  9% 42.3M 3s
    ##  12950K .......... .......... .......... .......... ..........  9% 83.4M 3s
    ##  13000K .......... .......... .......... .......... ..........  9% 26.4M 3s
    ##  13050K .......... .......... .......... .......... ..........  9% 37.5M 3s
    ##  13100K .......... .......... .......... .......... ..........  9% 38.1M 3s
    ##  13150K .......... .......... .......... .......... ..........  9% 26.1M 3s
    ##  13200K .......... .......... .......... .......... ..........  9% 91.1M 3s
    ##  13250K .......... .......... .......... .......... ..........  9% 25.8M 3s
    ##  13300K .......... .......... .......... .......... ..........  9%  233K 5s
    ##  13350K .......... .......... .......... .......... ..........  9% 11.0M 5s
    ##  13400K .......... .......... .......... .......... ..........  9%  103M 5s
    ##  13450K .......... .......... .......... .......... .......... 10%  100M 5s
    ##  13500K .......... .......... .......... .......... .......... 10%  103M 5s
    ##  13550K .......... .......... .......... .......... .......... 10% 87.8M 5s
    ##  13600K .......... .......... .......... .......... .......... 10% 87.7M 5s
    ##  13650K .......... .......... .......... .......... .......... 10%  113M 5s
    ##  13700K .......... .......... .......... .......... .......... 10%  104M 5s
    ##  13750K .......... .......... .......... .......... .......... 10% 73.2M 5s
    ##  13800K .......... .......... .......... .......... .......... 10%  104M 5s
    ##  13850K .......... .......... .......... .......... .......... 10%  110M 5s
    ##  13900K .......... .......... .......... .......... .......... 10%  108M 5s
    ##  13950K .......... .......... .......... .......... .......... 10%  105M 5s
    ##  14000K .......... .......... .......... .......... .......... 10% 98.5M 5s
    ##  14050K .......... .......... .......... .......... .......... 10%  100M 5s
    ##  14100K .......... .......... .......... .......... .......... 10% 12.2M 5s
    ##  14150K .......... .......... .......... .......... .......... 10% 99.8M 5s
    ##  14200K .......... .......... .......... .......... .......... 10% 85.3M 5s
    ##  14250K .......... .......... .......... .......... .......... 10% 18.6M 5s
    ##  14300K .......... .......... .......... .......... .......... 10% 76.9M 5s
    ##  14350K .......... .......... .......... .......... .......... 10% 17.7M 5s
    ##  14400K .......... .......... .......... .......... .......... 10% 72.8M 5s
    ##  14450K .......... .......... .......... .......... .......... 10% 73.2M 5s
    ##  14500K .......... .......... .......... .......... .......... 10% 22.5M 5s
    ##  14550K .......... .......... .......... .......... .......... 10% 39.4M 5s
    ##  14600K .......... .......... .......... .......... .......... 10% 78.2M 5s
    ##  14650K .......... .......... .......... .......... .......... 10% 29.2M 5s
    ##  14700K .......... .......... .......... .......... .......... 10% 31.1M 5s
    ##  14750K .......... .......... .......... .......... .......... 10% 29.5M 5s
    ##  14800K .......... .......... .......... .......... .......... 11% 41.6M 5s
    ##  14850K .......... .......... .......... .......... .......... 11% 46.0M 5s
    ##  14900K .......... .......... .......... .......... .......... 11% 29.0M 5s
    ##  14950K .......... .......... .......... .......... .......... 11% 14.5M 5s
    ##  15000K .......... .......... .......... .......... .......... 11% 45.2M 5s
    ##  15050K .......... .......... .......... .......... .......... 11% 18.0M 5s
    ##  15100K .......... .......... .......... .......... .......... 11% 38.3M 5s
    ##  15150K .......... .......... .......... .......... .......... 11% 56.0M 5s
    ##  15200K .......... .......... .......... .......... .......... 11% 11.0M 5s
    ##  15250K .......... .......... .......... .......... .......... 11% 41.7M 5s
    ##  15300K .......... .......... .......... .......... .......... 11% 48.6M 5s
    ##  15350K .......... .......... .......... .......... .......... 11% 53.6M 5s
    ##  15400K .......... .......... .......... .......... .......... 11% 27.6M 5s
    ##  15450K .......... .......... .......... .......... .......... 11% 24.6M 5s
    ##  15500K .......... .......... .......... .......... .......... 11% 15.4M 5s
    ##  15550K .......... .......... .......... .......... .......... 11% 41.9M 5s
    ##  15600K .......... .......... .......... .......... .......... 11% 48.6M 5s
    ##  15650K .......... .......... .......... .......... .......... 11% 51.9M 5s
    ##  15700K .......... .......... .......... .......... .......... 11% 15.4M 5s
    ##  15750K .......... .......... .......... .......... .......... 11% 24.7M 5s
    ##  15800K .......... .......... .......... .......... .......... 11% 49.9M 5s
    ##  15850K .......... .......... .......... .......... .......... 11% 46.9M 5s
    ##  15900K .......... .......... .......... .......... .......... 11% 25.6M 5s
    ##  15950K .......... .......... .......... .......... .......... 11% 36.4M 5s
    ##  16000K .......... .......... .......... .......... .......... 11% 47.8M 5s
    ##  16050K .......... .......... .......... .......... .......... 11% 57.4M 5s
    ##  16100K .......... .......... .......... .......... .......... 11% 33.8M 5s
    ##  16150K .......... .......... .......... .......... .......... 12% 45.1M 5s
    ##  16200K .......... .......... .......... .......... .......... 12% 44.0M 5s
    ##  16250K .......... .......... .......... .......... .......... 12% 43.8M 5s
    ##  16300K .......... .......... .......... .......... .......... 12% 51.8M 5s
    ##  16350K .......... .......... .......... .......... .......... 12% 53.8M 5s
    ##  16400K .......... .......... .......... .......... .......... 12% 53.5M 5s
    ##  16450K .......... .......... .......... .......... .......... 12% 6.37M 5s
    ##  16500K .......... .......... .......... .......... .......... 12% 47.8M 5s
    ##  16550K .......... .......... .......... .......... .......... 12% 49.4M 5s
    ##  16600K .......... .......... .......... .......... .......... 12% 53.4M 5s
    ##  16650K .......... .......... .......... .......... .......... 12% 45.1M 5s
    ##  16700K .......... .......... .......... .......... .......... 12% 17.4M 5s
    ##  16750K .......... .......... .......... .......... .......... 12% 43.0M 5s
    ##  16800K .......... .......... .......... .......... .......... 12% 44.2M 5s
    ##  16850K .......... .......... .......... .......... .......... 12% 42.2M 5s
    ##  16900K .......... .......... .......... .......... .......... 12% 49.1M 5s
    ##  16950K .......... .......... .......... .......... .......... 12% 52.5M 5s
    ##  17000K .......... .......... .......... .......... .......... 12% 50.8M 5s
    ##  17050K .......... .......... .......... .......... .......... 12% 57.4M 5s
    ##  17100K .......... .......... .......... .......... .......... 12% 23.0M 5s
    ##  17150K .......... .......... .......... .......... .......... 12% 48.8M 5s
    ##  17200K .......... .......... .......... .......... .......... 12% 56.2M 5s
    ##  17250K .......... .......... .......... .......... .......... 12% 66.4M 5s
    ##  17300K .......... .......... .......... .......... .......... 12% 16.0M 5s
    ##  17350K .......... .......... .......... .......... .......... 12% 49.7M 5s
    ##  17400K .......... .......... .......... .......... .......... 12% 59.5M 5s
    ##  17450K .......... .......... .......... .......... .......... 12% 60.4M 5s
    ##  17500K .......... .......... .......... .......... .......... 13% 24.9M 5s
    ##  17550K .......... .......... .......... .......... .......... 13% 28.4M 5s
    ##  17600K .......... .......... .......... .......... .......... 13% 52.5M 5s
    ##  17650K .......... .......... .......... .......... .......... 13% 56.4M 5s
    ##  17700K .......... .......... .......... .......... .......... 13% 65.9M 5s
    ##  17750K .......... .......... .......... .......... .......... 13% 22.8M 5s
    ##  17800K .......... .......... .......... .......... .......... 13% 41.8M 5s
    ##  17850K .......... .......... .......... .......... .......... 13% 48.4M 5s
    ##  17900K .......... .......... .......... .......... .......... 13% 40.2M 5s
    ##  17950K .......... .......... .......... .......... .......... 13% 51.3M 5s
    ##  18000K .......... .......... .......... .......... .......... 13% 26.9M 5s
    ##  18050K .......... .......... .......... .......... .......... 13% 48.6M 5s
    ##  18100K .......... .......... .......... .......... .......... 13% 42.9M 5s
    ##  18150K .......... .......... .......... .......... .......... 13% 48.5M 5s
    ##  18200K .......... .......... .......... .......... .......... 13% 42.5M 5s
    ##  18250K .......... .......... .......... .......... .......... 13% 51.2M 5s
    ##  18300K .......... .......... .......... .......... .......... 13% 42.0M 5s
    ##  18350K .......... .......... .......... .......... .......... 13% 46.4M 5s
    ##  18400K .......... .......... .......... .......... .......... 13% 50.3M 5s
    ##  18450K .......... .......... .......... .......... .......... 13% 35.3M 4s
    ##  18500K .......... .......... .......... .......... .......... 13% 43.0M 4s
    ##  18550K .......... .......... .......... .......... .......... 13% 53.6M 4s
    ##  18600K .......... .......... .......... .......... .......... 13% 51.2M 4s
    ##  18650K .......... .......... .......... .......... .......... 13% 22.3M 4s
    ##  18700K .......... .......... .......... .......... .......... 13% 45.9M 4s
    ##  18750K .......... .......... .......... .......... .......... 13% 47.3M 4s
    ##  18800K .......... .......... .......... .......... .......... 13% 26.5M 4s
    ##  18850K .......... .......... .......... .......... .......... 14% 41.7M 4s
    ##  18900K .......... .......... .......... .......... .......... 14% 45.0M 4s
    ##  18950K .......... .......... .......... .......... .......... 14% 57.5M 4s
    ##  19000K .......... .......... .......... .......... .......... 14% 52.6M 4s
    ##  19050K .......... .......... .......... .......... .......... 14% 29.9M 4s
    ##  19100K .......... .......... .......... .......... .......... 14% 40.1M 4s
    ##  19150K .......... .......... .......... .......... .......... 14% 47.6M 4s
    ##  19200K .......... .......... .......... .......... .......... 14% 55.6M 4s
    ##  19250K .......... .......... .......... .......... .......... 14% 63.8M 4s
    ##  19300K .......... .......... .......... .......... .......... 14% 45.1M 4s
    ##  19350K .......... .......... .......... .......... .......... 14% 45.9M 4s
    ##  19400K .......... .......... .......... .......... .......... 14% 47.6M 4s
    ##  19450K .......... .......... .......... .......... .......... 14% 58.3M 4s
    ##  19500K .......... .......... .......... .......... .......... 14% 53.5M 4s
    ##  19550K .......... .......... .......... .......... .......... 14% 50.5M 4s
    ##  19600K .......... .......... .......... .......... .......... 14% 52.0M 4s
    ##  19650K .......... .......... .......... .......... .......... 14% 57.8M 4s
    ##  19700K .......... .......... .......... .......... .......... 14% 54.3M 4s
    ##  19750K .......... .......... .......... .......... .......... 14% 52.2M 4s
    ##  19800K .......... .......... .......... .......... .......... 14% 60.7M 4s
    ##  19850K .......... .......... .......... .......... .......... 14% 63.0M 4s
    ##  19900K .......... .......... .......... .......... .......... 14% 60.4M 4s
    ##  19950K .......... .......... .......... .......... .......... 14% 52.6M 4s
    ##  20000K .......... .......... .......... .......... .......... 14% 45.7M 4s
    ##  20050K .......... .......... .......... .......... .......... 14% 59.5M 4s
    ##  20100K .......... .......... .......... .......... .......... 14% 53.1M 4s
    ##  20150K .......... .......... .......... .......... .......... 14% 62.4M 4s
    ##  20200K .......... .......... .......... .......... .......... 15% 56.7M 4s
    ##  20250K .......... .......... .......... .......... .......... 15% 61.6M 4s
    ##  20300K .......... .......... .......... .......... .......... 15% 62.3M 4s
    ##  20350K .......... .......... .......... .......... .......... 15% 74.2M 4s
    ##  20400K .......... .......... .......... .......... .......... 15% 62.4M 4s
    ##  20450K .......... .......... .......... .......... .......... 15% 52.5M 4s
    ##  20500K .......... .......... .......... .......... .......... 15% 59.4M 4s
    ##  20550K .......... .......... .......... .......... .......... 15% 66.5M 4s
    ##  20600K .......... .......... .......... .......... .......... 15% 51.9M 4s
    ##  20650K .......... .......... .......... .......... .......... 15% 64.3M 4s
    ##  20700K .......... .......... .......... .......... .......... 15% 53.6M 4s
    ##  20750K .......... .......... .......... .......... .......... 15% 60.9M 4s
    ##  20800K .......... .......... .......... .......... .......... 15% 50.0M 4s
    ##  20850K .......... .......... .......... .......... .......... 15% 78.6M 4s
    ##  20900K .......... .......... .......... .......... .......... 15% 63.0M 4s
    ##  20950K .......... .......... .......... .......... .......... 15% 79.3M 4s
    ##  21000K .......... .......... .......... .......... .......... 15% 67.5M 4s
    ##  21050K .......... .......... .......... .......... .......... 15% 59.1M 4s
    ##  21100K .......... .......... .......... .......... .......... 15% 58.9M 4s
    ##  21150K .......... .......... .......... .......... .......... 15% 76.3M 4s
    ##  21200K .......... .......... .......... .......... .......... 15% 61.7M 4s
    ##  21250K .......... .......... .......... .......... .......... 15% 65.4M 4s
    ##  21300K .......... .......... .......... .......... .......... 15% 17.1M 4s
    ##  21350K .......... .......... .......... .......... .......... 15% 59.2M 4s
    ##  21400K .......... .......... .......... .......... .......... 15% 51.7M 4s
    ##  21450K .......... .......... .......... .......... .......... 15% 62.3M 4s
    ##  21500K .......... .......... .......... .......... .......... 15% 69.3M 4s
    ##  21550K .......... .......... .......... .......... .......... 16% 24.4M 4s
    ##  21600K .......... .......... .......... .......... .......... 16% 53.1M 4s
    ##  21650K .......... .......... .......... .......... .......... 16% 59.0M 4s
    ##  21700K .......... .......... .......... .......... .......... 16% 56.4M 4s
    ##  21750K .......... .......... .......... .......... .......... 16% 76.5M 4s
    ##  21800K .......... .......... .......... .......... .......... 16% 36.2M 4s
    ##  21850K .......... .......... .......... .......... .......... 16% 47.7M 4s
    ##  21900K .......... .......... .......... .......... .......... 16% 44.0M 4s
    ##  21950K .......... .......... .......... .......... .......... 16% 44.9M 4s
    ##  22000K .......... .......... .......... .......... .......... 16% 48.3M 4s
    ##  22050K .......... .......... .......... .......... .......... 16% 48.6M 4s
    ##  22100K .......... .......... .......... .......... .......... 16% 45.0M 4s
    ##  22150K .......... .......... .......... .......... .......... 16% 57.6M 4s
    ##  22200K .......... .......... .......... .......... .......... 16% 49.9M 4s
    ##  22250K .......... .......... .......... .......... .......... 16% 7.35M 4s
    ##  22300K .......... .......... .......... .......... .......... 16% 48.2M 4s
    ##  22350K .......... .......... .......... .......... .......... 16% 61.4M 4s
    ##  22400K .......... .......... .......... .......... .......... 16% 53.4M 4s
    ##  22450K .......... .......... .......... .......... .......... 16% 57.8M 4s
    ##  22500K .......... .......... .......... .......... .......... 16% 46.3M 4s
    ##  22550K .......... .......... .......... .......... .......... 16% 43.9M 4s
    ##  22600K .......... .......... .......... .......... .......... 16% 50.1M 4s
    ##  22650K .......... .......... .......... .......... .......... 16% 54.2M 4s
    ##  22700K .......... .......... .......... .......... .......... 16% 54.3M 4s
    ##  22750K .......... .......... .......... .......... .......... 16% 58.6M 4s
    ##  22800K .......... .......... .......... .......... .......... 16% 48.6M 4s
    ##  22850K .......... .......... .......... .......... .......... 16% 15.8M 4s
    ##  22900K .......... .......... .......... .......... .......... 17% 56.9M 4s
    ##  22950K .......... .......... .......... .......... .......... 17% 51.5M 4s
    ##  23000K .......... .......... .......... .......... .......... 17% 56.6M 4s
    ##  23050K .......... .......... .......... .......... .......... 17%  221K 5s
    ##  23100K .......... .......... .......... .......... .......... 17% 45.7M 5s
    ##  23150K .......... .......... .......... .......... .......... 17% 51.8M 5s
    ##  23200K .......... .......... .......... .......... .......... 17% 54.2M 5s
    ##  23250K .......... .......... .......... .......... .......... 17% 57.5M 5s
    ##  23300K .......... .......... .......... .......... .......... 17% 53.0M 5s
    ##  23350K .......... .......... .......... .......... .......... 17% 59.1M 5s
    ##  23400K .......... .......... .......... .......... .......... 17% 55.8M 5s
    ##  23450K .......... .......... .......... .......... .......... 17% 34.2M 5s
    ##  23500K .......... .......... .......... .......... .......... 17% 52.8M 5s
    ##  23550K .......... .......... .......... .......... .......... 17% 60.9M 5s
    ##  23600K .......... .......... .......... .......... .......... 17% 57.7M 5s
    ##  23650K .......... .......... .......... .......... .......... 17% 60.7M 5s
    ##  23700K .......... .......... .......... .......... .......... 17% 57.5M 5s
    ##  23750K .......... .......... .......... .......... .......... 17% 54.3M 5s
    ##  23800K .......... .......... .......... .......... .......... 17% 51.9M 5s
    ##  23850K .......... .......... .......... .......... .......... 17% 61.6M 5s
    ##  23900K .......... .......... .......... .......... .......... 17% 58.5M 5s
    ##  23950K .......... .......... .......... .......... .......... 17% 65.1M 5s
    ##  24000K .......... .......... .......... .......... .......... 17% 54.4M 5s
    ##  24050K .......... .......... .......... .......... .......... 17% 60.0M 5s
    ##  24100K .......... .......... .......... .......... .......... 17% 61.3M 5s
    ##  24150K .......... .......... .......... .......... .......... 17% 65.8M 5s
    ##  24200K .......... .......... .......... .......... .......... 17% 62.5M 5s
    ##  24250K .......... .......... .......... .......... .......... 18% 59.5M 5s
    ##  24300K .......... .......... .......... .......... .......... 18% 59.1M 5s
    ##  24350K .......... .......... .......... .......... .......... 18% 59.1M 5s
    ##  24400K .......... .......... .......... .......... .......... 18% 61.4M 5s
    ##  24450K .......... .......... .......... .......... .......... 18% 51.1M 5s
    ##  24500K .......... .......... .......... .......... .......... 18% 56.6M 5s
    ##  24550K .......... .......... .......... .......... .......... 18% 70.3M 5s
    ##  24600K .......... .......... .......... .......... .......... 18% 64.3M 5s
    ##  24650K .......... .......... .......... .......... .......... 18% 51.2M 5s
    ##  24700K .......... .......... .......... .......... .......... 18% 55.7M 5s
    ##  24750K .......... .......... .......... .......... .......... 18% 58.2M 5s
    ##  24800K .......... .......... .......... .......... .......... 18% 64.6M 5s
    ##  24850K .......... .......... .......... .......... .......... 18% 73.6M 5s
    ##  24900K .......... .......... .......... .......... .......... 18% 64.7M 5s
    ##  24950K .......... .......... .......... .......... .......... 18% 60.1M 5s
    ##  25000K .......... .......... .......... .......... .......... 18% 54.9M 5s
    ##  25050K .......... .......... .......... .......... .......... 18% 71.1M 5s
    ##  25100K .......... .......... .......... .......... .......... 18% 66.7M 5s
    ##  25150K .......... .......... .......... .......... .......... 18% 70.8M 5s
    ##  25200K .......... .......... .......... .......... .......... 18% 66.8M 5s
    ##  25250K .......... .......... .......... .......... .......... 18% 63.5M 5s
    ##  25300K .......... .......... .......... .......... .......... 18% 64.7M 5s
    ##  25350K .......... .......... .......... .......... .......... 18% 70.2M 5s
    ##  25400K .......... .......... .......... .......... .......... 18% 75.4M 5s
    ##  25450K .......... .......... .......... .......... .......... 18% 1.21M 5s
    ##  25500K .......... .......... .......... .......... .......... 18% 2.85M 5s
    ##  25550K .......... .......... .......... .......... .......... 18% 20.6M 5s
    ##  25600K .......... .......... .......... .......... .......... 19% 37.4M 5s
    ##  25650K .......... .......... .......... .......... .......... 19% 50.8M 5s
    ##  25700K .......... .......... .......... .......... .......... 19% 11.5M 5s
    ##  25750K .......... .......... .......... .......... .......... 19% 41.9M 5s
    ##  25800K .......... .......... .......... .......... .......... 19% 13.9M 5s
    ##  25850K .......... .......... .......... .......... .......... 19% 9.21M 5s
    ##  25900K .......... .......... .......... .......... .......... 19% 50.0M 5s
    ##  25950K .......... .......... .......... .......... .......... 19% 8.26M 5s
    ##  26000K .......... .......... .......... .......... .......... 19% 14.4M 5s
    ##  26050K .......... .......... .......... .......... .......... 19% 10.5M 5s
    ##  26100K .......... .......... .......... .......... .......... 19% 6.24M 5s
    ##  26150K .......... .......... .......... .......... .......... 19% 5.68M 5s
    ##  26200K .......... .......... .......... .......... .......... 19% 1.62M 5s
    ##  26250K .......... .......... .......... .......... .......... 19% 46.4M 5s
    ##  26300K .......... .......... .......... .......... .......... 19% 67.9M 5s
    ##  26350K .......... .......... .......... .......... .......... 19% 61.3M 5s
    ##  26400K .......... .......... .......... .......... .......... 19% 61.3M 5s
    ##  26450K .......... .......... .......... .......... .......... 19% 70.2M 5s
    ##  26500K .......... .......... .......... .......... .......... 19% 50.3M 5s
    ##  26550K .......... .......... .......... .......... .......... 19% 6.67M 5s
    ##  26600K .......... .......... .......... .......... .......... 19% 2.81M 5s
    ##  26650K .......... .......... .......... .......... .......... 19% 1.28M 5s
    ##  26700K .......... .......... .......... .......... .......... 19% 72.5M 5s
    ##  26750K .......... .......... .......... .......... .......... 19% 74.3M 5s
    ##  26800K .......... .......... .......... .......... .......... 19% 80.0M 5s
    ##  26850K .......... .......... .......... .......... .......... 19% 5.39M 5s
    ##  26900K .......... .......... .......... .......... .......... 20% 2.67M 5s
    ##  26950K .......... .......... .......... .......... .......... 20% 10.3M 5s
    ##  27000K .......... .......... .......... .......... .......... 20% 5.94M 5s
    ##  27050K .......... .......... .......... .......... .......... 20% 4.17M 5s
    ##  27100K .......... .......... .......... .......... .......... 20% 8.17M 5s
    ##  27150K .......... .......... .......... .......... .......... 20% 10.7M 5s
    ##  27200K .......... .......... .......... .......... .......... 20% 3.64M 5s
    ##  27250K .......... .......... .......... .......... .......... 20% 7.98M 5s
    ##  27300K .......... .......... .......... .......... .......... 20% 5.29M 5s
    ##  27350K .......... .......... .......... .......... .......... 20% 7.36M 5s
    ##  27400K .......... .......... .......... .......... .......... 20% 4.48M 5s
    ##  27450K .......... .......... .......... .......... .......... 20% 8.81M 5s
    ##  27500K .......... .......... .......... .......... .......... 20% 9.14M 5s
    ##  27550K .......... .......... .......... .......... .......... 20% 2.72M 6s
    ##  27600K .......... .......... .......... .......... .......... 20% 3.85M 6s
    ##  27650K .......... .......... .......... .......... .......... 20% 12.0M 6s
    ##  27700K .......... .......... .......... .......... .......... 20% 10.8M 6s
    ##  27750K .......... .......... .......... .......... .......... 20% 12.1M 6s
    ##  27800K .......... .......... .......... .......... .......... 20% 5.20M 6s
    ##  27850K .......... .......... .......... .......... .......... 20%  881K 6s
    ##  27900K .......... .......... .......... .......... .......... 20% 7.72M 6s
    ##  27950K .......... .......... .......... .......... .......... 20% 80.5M 6s
    ##  28000K .......... .......... .......... .......... .......... 20% 74.5M 6s
    ##  28050K .......... .......... .......... .......... .......... 20% 21.6M 6s
    ##  28100K .......... .......... .......... .......... .......... 20% 1.64M 6s
    ##  28150K .......... .......... .......... .......... .......... 20% 7.04M 6s
    ##  28200K .......... .......... .......... .......... .......... 20% 11.5M 6s
    ##  28250K .......... .......... .......... .......... .......... 21% 6.21M 6s
    ##  28300K .......... .......... .......... .......... .......... 21%  100K 8s
    ##  28350K .......... .......... .......... .......... .......... 21% 88.0M 8s
    ##  28400K .......... .......... .......... .......... .......... 21% 84.4M 8s
    ##  28450K .......... .......... .......... .......... .......... 21% 83.0M 8s
    ##  28500K .......... .......... .......... .......... .......... 21% 73.7M 8s
    ##  28550K .......... .......... .......... .......... .......... 21% 83.7M 8s
    ##  28600K .......... .......... .......... .......... .......... 21% 79.3M 8s
    ##  28650K .......... .......... .......... .......... .......... 21% 99.6M 8s
    ##  28700K .......... .......... .......... .......... .......... 21% 78.1M 8s
    ##  28750K .......... .......... .......... .......... .......... 21% 96.3M 8s
    ##  28800K .......... .......... .......... .......... .......... 21% 73.5M 8s
    ##  28850K .......... .......... .......... .......... .......... 21% 65.3M 8s
    ##  28900K .......... .......... .......... .......... .......... 21% 75.0M 8s
    ##  28950K .......... .......... .......... .......... .......... 21% 98.9M 8s
    ##  29000K .......... .......... .......... .......... .......... 21% 70.7M 8s
    ##  29050K .......... .......... .......... .......... .......... 21% 92.2M 8s
    ##  29100K .......... .......... .......... .......... .......... 21% 80.5M 8s
    ##  29150K .......... .......... .......... .......... .......... 21% 83.1M 8s
    ##  29200K .......... .......... .......... .......... .......... 21% 76.6M 8s
    ##  29250K .......... .......... .......... .......... .......... 21%  111M 8s
    ##  29300K .......... .......... .......... .......... .......... 21% 3.57M 8s
    ##  29350K .......... .......... .......... .......... .......... 21% 1.65M 8s
    ##  29400K .......... .......... .......... .......... .......... 21% 2.33M 8s
    ##  29450K .......... .......... .......... .......... .......... 21% 3.54M 8s
    ##  29500K .......... .......... .......... .......... .......... 21% 1.66M 8s
    ##  29550K .......... .......... .......... .......... .......... 21% 4.36M 8s
    ##  29600K .......... .......... .......... .......... .......... 22% 6.30M 8s
    ##  29650K .......... .......... .......... .......... .......... 22% 4.78M 8s
    ##  29700K .......... .......... .......... .......... .......... 22% 2.09M 8s
    ##  29750K .......... .......... .......... .......... .......... 22% 4.46M 8s
    ##  29800K .......... .......... .......... .......... .......... 22% 3.46M 8s
    ##  29850K .......... .......... .......... .......... .......... 22% 11.6M 8s
    ##  29900K .......... .......... .......... .......... .......... 22% 6.80M 8s
    ##  29950K .......... .......... .......... .......... .......... 22% 7.37M 8s
    ##  30000K .......... .......... .......... .......... .......... 22% 5.25M 8s
    ##  30050K .......... .......... .......... .......... .......... 22% 3.45M 8s
    ##  30100K .......... .......... .......... .......... .......... 22% 10.8M 8s
    ##  30150K .......... .......... .......... .......... .......... 22% 3.43M 8s
    ##  30200K .......... .......... .......... .......... .......... 22% 1.07M 8s
    ##  30250K .......... .......... .......... .......... .......... 22% 95.6M 8s
    ##  30300K .......... .......... .......... .......... .......... 22% 77.2M 8s
    ##  30350K .......... .......... .......... .......... .......... 22% 94.8M 8s
    ##  30400K .......... .......... .......... .......... .......... 22% 19.3M 8s
    ##  30450K .......... .......... .......... .......... .......... 22% 3.37M 8s
    ##  30500K .......... .......... .......... .......... .......... 22% 8.96M 8s
    ##  30550K .......... .......... .......... .......... .......... 22% 10.2M 8s
    ##  30600K .......... .......... .......... .......... .......... 22% 5.26M 8s
    ##  30650K .......... .......... .......... .......... .......... 22% 2.65M 8s
    ##  30700K .......... .......... .......... .......... .......... 22% 7.94M 8s
    ##  30750K .......... .......... .......... .......... .......... 22%  222K 9s
    ##  30800K .......... .......... .......... .......... .......... 22% 10.7M 9s
    ##  30850K .......... .......... .......... .......... .......... 22% 10.2M 9s
    ##  30900K .......... .......... .......... .......... .......... 22% 11.2M 9s
    ##  30950K .......... .......... .......... .......... .......... 23% 13.2M 9s
    ##  31000K .......... .......... .......... .......... .......... 23% 7.66M 9s
    ##  31050K .......... .......... .......... .......... .......... 23% 5.83M 9s
    ##  31100K .......... .......... .......... .......... .......... 23% 10.2M 9s
    ##  31150K .......... .......... .......... .......... .......... 23% 7.25M 9s
    ##  31200K .......... .......... .......... .......... .......... 23% 5.50M 9s
    ##  31250K .......... .......... .......... .......... .......... 23% 6.82M 9s
    ##  31300K .......... .......... .......... .......... .......... 23% 3.56M 9s
    ##  31350K .......... .......... .......... .......... .......... 23% 9.62M 9s
    ##  31400K .......... .......... .......... .......... .......... 23% 12.9M 9s
    ##  31450K .......... .......... .......... .......... .......... 23% 12.8M 9s
    ##  31500K .......... .......... .......... .......... .......... 23% 5.91M 9s
    ##  31550K .......... .......... .......... .......... .......... 23% 6.24M 9s
    ##  31600K .......... .......... .......... .......... .......... 23% 5.84M 9s
    ##  31650K .......... .......... .......... .......... .......... 23% 6.84M 9s
    ##  31700K .......... .......... .......... .......... .......... 23% 68.9K 11s
    ##  31750K .......... .......... .......... .......... .......... 23% 2.46M 11s
    ##  31800K .......... .......... .......... .......... .......... 23% 1.11M 12s
    ##  31850K .......... .......... .......... .......... .......... 23% 3.86M 12s
    ##  31900K .......... .......... .......... .......... .......... 23% 2.44M 12s
    ##  31950K .......... .......... .......... .......... .......... 23% 4.14M 12s
    ##  32000K .......... .......... .......... .......... .......... 23% 1.82M 12s
    ##  32050K .......... .......... .......... .......... .......... 23% 3.35M 12s
    ##  32100K .......... .......... .......... .......... .......... 23% 5.08M 12s
    ##  32150K .......... .......... .......... .......... .......... 23% 4.98M 12s
    ##  32200K .......... .......... .......... .......... .......... 23% 5.57M 12s
    ##  32250K .......... .......... .......... .......... .......... 23% 2.45M 12s
    ##  32300K .......... .......... .......... .......... .......... 24% 3.62M 12s
    ##  32350K .......... .......... .......... .......... .......... 24% 7.06M 12s
    ##  32400K .......... .......... .......... .......... .......... 24% 6.14M 12s
    ##  32450K .......... .......... .......... .......... .......... 24% 3.47M 12s
    ##  32500K .......... .......... .......... .......... .......... 24% 1.90M 12s
    ##  32550K .......... .......... .......... .......... .......... 24%  816K 12s
    ##  32600K .......... .......... .......... .......... .......... 24% 89.4M 12s
    ##  32650K .......... .......... .......... .......... .......... 24%  109M 12s
    ##  32700K .......... .......... .......... .......... .......... 24%  180K 13s
    ##  32750K .......... .......... .......... .......... .......... 24% 71.1M 13s
    ##  32800K .......... .......... .......... .......... .......... 24% 88.9M 13s
    ##  32850K .......... .......... .......... .......... .......... 24% 20.3M 13s
    ##  32900K .......... .......... .......... .......... .......... 24% 7.24M 13s
    ##  32950K .......... .......... .......... .......... .......... 24% 4.31M 13s
    ##  33000K .......... .......... .......... .......... .......... 24% 3.60M 13s
    ##  33050K .......... .......... .......... .......... .......... 24% 6.11M 13s
    ##  33100K .......... .......... .......... .......... .......... 24% 4.70M 13s
    ##  33150K .......... .......... .......... .......... .......... 24% 2.10M 13s
    ##  33200K .......... .......... .......... .......... .......... 24% 2.19M 13s
    ##  33250K .......... .......... .......... .......... .......... 24% 10.3M 13s
    ##  33300K .......... .......... .......... .......... .......... 24% 6.82M 13s
    ##  33350K .......... .......... .......... .......... .......... 24% 11.7M 13s
    ##  33400K .......... .......... .......... .......... .......... 24% 7.89M 13s
    ##  33450K .......... .......... .......... .......... .......... 24% 8.23M 13s
    ##  33500K .......... .......... .......... .......... .......... 24% 3.89M 13s
    ##  33550K .......... .......... .......... .......... .......... 24% 5.25M 13s
    ##  33600K .......... .......... .......... .......... .......... 24% 3.55M 13s
    ##  33650K .......... .......... .......... .......... .......... 25% 6.39M 13s
    ##  33700K .......... .......... .......... .......... .......... 25% 5.66M 13s
    ##  33750K .......... .......... .......... .......... .......... 25% 12.1M 13s
    ##  33800K .......... .......... .......... .......... .......... 25% 12.7M 13s
    ##  33850K .......... .......... .......... .......... .......... 25% 13.7M 13s
    ##  33900K .......... .......... .......... .......... .......... 25% 3.46M 13s
    ##  33950K .......... .......... .......... .......... .......... 25% 5.79M 13s
    ##  34000K .......... .......... .......... .......... .......... 25% 5.78M 13s
    ##  34050K .......... .......... .......... .......... .......... 25% 5.33M 13s
    ##  34100K .......... .......... .......... .......... .......... 25% 9.30M 13s
    ##  34150K .......... .......... .......... .......... .......... 25% 14.1M 13s
    ##  34200K .......... .......... .......... .......... .......... 25% 13.0M 13s
    ##  34250K .......... .......... .......... .......... .......... 25% 13.3M 13s
    ##  34300K .......... .......... .......... .......... .......... 25% 9.82M 13s
    ##  34350K .......... .......... .......... .......... .......... 25% 8.18M 13s
    ##  34400K .......... .......... .......... .......... .......... 25% 13.6M 13s
    ##  34450K .......... .......... .......... .......... .......... 25% 16.2M 13s
    ##  34500K .......... .......... .......... .......... .......... 25% 63.5M 13s
    ##  34550K .......... .......... .......... .......... .......... 25% 13.1M 13s
    ##  34600K .......... .......... .......... .......... .......... 25% 13.1M 13s
    ##  34650K .......... .......... .......... .......... .......... 25% 14.1M 13s
    ##  34700K .......... .......... .......... .......... .......... 25% 61.8M 13s
    ##  34750K .......... .......... .......... .......... .......... 25% 14.8M 13s
    ##  34800K .......... .......... .......... .......... .......... 25% 19.5M 13s
    ##  34850K .......... .......... .......... .......... .......... 25% 16.8M 13s
    ##  34900K .......... .......... .......... .......... .......... 25% 24.8M 13s
    ##  34950K .......... .......... .......... .......... .......... 25% 16.1M 13s
    ##  35000K .......... .......... .......... .......... .......... 26% 14.0M 13s
    ##  35050K .......... .......... .......... .......... .......... 26% 31.0M 13s
    ##  35100K .......... .......... .......... .......... .......... 26% 16.0M 12s
    ##  35150K .......... .......... .......... .......... .......... 26% 14.5M 12s
    ##  35200K .......... .......... .......... .......... .......... 26% 15.5M 12s
    ##  35250K .......... .......... .......... .......... .......... 26% 90.6M 12s
    ##  35300K .......... .......... .......... .......... .......... 26% 13.4M 12s
    ##  35350K .......... .......... .......... .......... .......... 26% 13.1M 12s
    ##  35400K .......... .......... .......... .......... .......... 26% 13.7M 12s
    ##  35450K .......... .......... .......... .......... .......... 26% 12.9M 12s
    ##  35500K .......... .......... .......... .......... .......... 26% 84.9M 12s
    ##  35550K .......... .......... .......... .......... .......... 26% 12.8M 12s
    ##  35600K .......... .......... .......... .......... .......... 26% 13.5M 12s
    ##  35650K .......... .......... .......... .......... .......... 26% 14.1M 12s
    ##  35700K .......... .......... .......... .......... .......... 26% 30.8M 12s
    ##  35750K .......... .......... .......... .......... .......... 26% 19.8M 12s
    ##  35800K .......... .......... .......... .......... .......... 26% 26.5M 12s
    ##  35850K .......... .......... .......... .......... .......... 26% 30.6M 12s
    ##  35900K .......... .......... .......... .......... .......... 26% 11.8M 12s
    ##  35950K .......... .......... .......... .......... .......... 26% 39.5M 12s
    ##  36000K .......... .......... .......... .......... .......... 26% 16.2M 12s
    ##  36050K .......... .......... .......... .......... .......... 26% 35.4M 12s
    ##  36100K .......... .......... .......... .......... .......... 26% 23.5M 12s
    ##  36150K .......... .......... .......... .......... .......... 26% 33.9M 12s
    ##  36200K .......... .......... .......... .......... .......... 26% 11.8M 12s
    ##  36250K .......... .......... .......... .......... .......... 26% 98.7M 12s
    ##  36300K .......... .......... .......... .......... .......... 26% 11.7M 12s
    ##  36350K .......... .......... .......... .......... .......... 27% 85.5M 12s
    ##  36400K .......... .......... .......... .......... .......... 27% 21.0M 12s
    ##  36450K .......... .......... .......... .......... .......... 27% 82.7M 12s
    ##  36500K .......... .......... .......... .......... .......... 27% 16.3M 12s
    ##  36550K .......... .......... .......... .......... .......... 27% 83.1M 12s
    ##  36600K .......... .......... .......... .......... .......... 27% 14.4M 12s
    ##  36650K .......... .......... .......... .......... .......... 27% 94.5M 12s
    ##  36700K .......... .......... .......... .......... .......... 27% 13.3M 12s
    ##  36750K .......... .......... .......... .......... .......... 27% 83.2M 12s
    ##  36800K .......... .......... .......... .......... .......... 27% 11.5M 12s
    ##  36850K .......... .......... .......... .......... .......... 27% 29.7M 12s
    ##  36900K .......... .......... .......... .......... .......... 27% 19.6M 12s
    ##  36950K .......... .......... .......... .......... .......... 27% 62.0M 12s
    ##  37000K .......... .......... .......... .......... .......... 27% 83.9M 12s
    ##  37050K .......... .......... .......... .......... .......... 27% 13.6M 12s
    ##  37100K .......... .......... .......... .......... .......... 27% 77.1M 12s
    ##  37150K .......... .......... .......... .......... .......... 27% 14.8M 12s
    ##  37200K .......... .......... .......... .......... .......... 27% 70.6M 12s
    ##  37250K .......... .......... .......... .......... .......... 27% 78.9M 12s
    ##  37300K .......... .......... .......... .......... .......... 27% 18.7M 12s
    ##  37350K .......... .......... .......... .......... .......... 27% 89.0M 12s
    ##  37400K .......... .......... .......... .......... .......... 27% 14.2M 12s
    ##  37450K .......... .......... .......... .......... .......... 27% 81.6M 12s
    ##  37500K .......... .......... .......... .......... .......... 27% 13.6M 12s
    ##  37550K .......... .......... .......... .......... .......... 27% 75.4M 12s
    ##  37600K .......... .......... .......... .......... .......... 27% 28.1M 12s
    ##  37650K .......... .......... .......... .......... .......... 27% 22.5M 12s
    ##  37700K .......... .......... .......... .......... .......... 28% 79.0M 12s
    ##  37750K .......... .......... .......... .......... .......... 28% 57.0M 12s
    ##  37800K .......... .......... .......... .......... .......... 28% 12.1M 12s
    ##  37850K .......... .......... .......... .......... .......... 28% 99.7M 12s
    ##  37900K .......... .......... .......... .......... .......... 28% 16.4M 12s
    ##  37950K .......... .......... .......... .......... .......... 28% 74.1M 12s
    ##  38000K .......... .......... .......... .......... .......... 28% 14.7M 12s
    ##  38050K .......... .......... .......... .......... .......... 28% 77.0M 12s
    ##  38100K .......... .......... .......... .......... .......... 28% 87.9M 11s
    ##  38150K .......... .......... .......... .......... .......... 28% 14.0M 11s
    ##  38200K .......... .......... .......... .......... .......... 28% 78.0M 11s
    ##  38250K .......... .......... .......... .......... .......... 28% 18.1M 11s
    ##  38300K .......... .......... .......... .......... .......... 28% 52.4M 11s
    ##  38350K .......... .......... .......... .......... .......... 28%  105M 11s
    ##  38400K .......... .......... .......... .......... .......... 28% 15.2M 11s
    ##  38450K .......... .......... .......... .......... .......... 28% 82.2M 11s
    ##  38500K .......... .......... .......... .......... .......... 28% 17.8M 11s
    ##  38550K .......... .......... .......... .......... .......... 28% 35.2M 11s
    ##  38600K .......... .......... .......... .......... .......... 28% 86.4M 11s
    ##  38650K .......... .......... .......... .......... .......... 28% 9.42M 11s
    ##  38700K .......... .......... .......... .......... .......... 28% 60.2M 11s
    ##  38750K .......... .......... .......... .......... .......... 28% 23.2M 11s
    ##  38800K .......... .......... .......... .......... .......... 28% 39.7M 11s
    ##  38850K .......... .......... .......... .......... .......... 28% 73.3M 11s
    ##  38900K .......... .......... .......... .......... .......... 28% 21.7M 11s
    ##  38950K .......... .......... .......... .......... .......... 28% 46.8M 11s
    ##  39000K .......... .......... .......... .......... .......... 28% 25.5M 11s
    ##  39050K .......... .......... .......... .......... .......... 29% 83.5M 11s
    ##  39100K .......... .......... .......... .......... .......... 29% 34.1M 11s
    ##  39150K .......... .......... .......... .......... .......... 29% 14.4M 11s
    ##  39200K .......... .......... .......... .......... .......... 29% 76.7M 11s
    ##  39250K .......... .......... .......... .......... .......... 29% 16.6M 11s
    ##  39300K .......... .......... .......... .......... .......... 29% 60.4M 11s
    ##  39350K .......... .......... .......... .......... .......... 29%  120M 11s
    ##  39400K .......... .......... .......... .......... .......... 29% 8.42M 11s
    ##  39450K .......... .......... .......... .......... .......... 29% 88.3M 11s
    ##  39500K .......... .......... .......... .......... .......... 29%  100M 11s
    ##  39550K .......... .......... .......... .......... .......... 29% 11.1M 11s
    ##  39600K .......... .......... .......... .......... .......... 29% 92.8M 11s
    ##  39650K .......... .......... .......... .......... .......... 29%  127M 11s
    ##  39700K .......... .......... .......... .......... .......... 29% 16.1M 11s
    ##  39750K .......... .......... .......... .......... .......... 29% 68.1M 11s
    ##  39800K .......... .......... .......... .......... .......... 29% 94.8M 11s
    ##  39850K .......... .......... .......... .......... .......... 29% 16.5M 11s
    ##  39900K .......... .......... .......... .......... .......... 29% 87.0M 11s
    ##  39950K .......... .......... .......... .......... .......... 29%  130M 11s
    ##  40000K .......... .......... .......... .......... .......... 29% 15.7M 11s
    ##  40050K .......... .......... .......... .......... .......... 29% 81.4M 11s
    ##  40100K .......... .......... .......... .......... .......... 29% 20.5M 11s
    ##  40150K .......... .......... .......... .......... .......... 29% 92.8M 11s
    ##  40200K .......... .......... .......... .......... .......... 29% 92.2M 11s
    ##  40250K .......... .......... .......... .......... .......... 29% 22.9M 11s
    ##  40300K .......... .......... .......... .......... .......... 29% 30.7M 11s
    ##  40350K .......... .......... .......... .......... .......... 29% 73.6M 11s
    ##  40400K .......... .......... .......... .......... .......... 30% 22.1M 11s
    ##  40450K .......... .......... .......... .......... .......... 30%  105M 11s
    ##  40500K .......... .......... .......... .......... .......... 30% 28.5M 11s
    ##  40550K .......... .......... .......... .......... .......... 30% 94.6M 11s
    ##  40600K .......... .......... .......... .......... .......... 30% 56.0M 11s
    ##  40650K .......... .......... .......... .......... .......... 30% 43.0M 11s
    ##  40700K .......... .......... .......... .......... .......... 30% 18.4M 11s
    ##  40750K .......... .......... .......... .......... .......... 30% 92.3M 11s
    ##  40800K .......... .......... .......... .......... .......... 30% 83.7M 11s
    ##  40850K .......... .......... .......... .......... .......... 30%  110M 11s
    ##  40900K .......... .......... .......... .......... .......... 30% 12.9M 11s
    ##  40950K .......... .......... .......... .......... .......... 30%  107M 11s
    ##  41000K .......... .......... .......... .......... .......... 30% 63.6M 11s
    ##  41050K .......... .......... .......... .......... .......... 30% 10.9M 11s
    ##  41100K .......... .......... .......... .......... .......... 30% 81.8M 11s
    ##  41150K .......... .......... .......... .......... .......... 30% 31.9M 11s
    ##  41200K .......... .......... .......... .......... .......... 30% 6.62M 11s
    ##  41250K .......... .......... .......... .......... .......... 30%  112M 10s
    ##  41300K .......... .......... .......... .......... .......... 30% 98.1M 10s
    ##  41350K .......... .......... .......... .......... .......... 30%  128M 10s
    ##  41400K .......... .......... .......... .......... .......... 30%  109M 10s
    ##  41450K .......... .......... .......... .......... .......... 30%  125M 10s
    ##  41500K .......... .......... .......... .......... .......... 30% 98.5M 10s
    ##  41550K .......... .......... .......... .......... .......... 30%  104M 10s
    ##  41600K .......... .......... .......... .......... .......... 30% 71.3M 10s
    ##  41650K .......... .......... .......... .......... .......... 30%  120M 10s
    ##  41700K .......... .......... .......... .......... .......... 30% 20.9M 10s
    ##  41750K .......... .......... .......... .......... .......... 31% 25.2M 10s
    ##  41800K .......... .......... .......... .......... .......... 31% 58.7M 10s
    ##  41850K .......... .......... .......... .......... .......... 31% 58.6M 10s
    ##  41900K .......... .......... .......... .......... .......... 31% 18.2M 10s
    ##  41950K .......... .......... .......... .......... .......... 31% 12.4M 10s
    ##  42000K .......... .......... .......... .......... .......... 31% 14.1M 10s
    ##  42050K .......... .......... .......... .......... .......... 31% 35.6M 10s
    ##  42100K .......... .......... .......... .......... .......... 31% 14.5M 10s
    ##  42150K .......... .......... .......... .......... .......... 31% 64.2M 10s
    ##  42200K .......... .......... .......... .......... .......... 31% 79.2M 10s
    ##  42250K .......... .......... .......... .......... .......... 31% 23.8M 10s
    ##  42300K .......... .......... .......... .......... .......... 31% 28.9M 10s
    ##  42350K .......... .......... .......... .......... .......... 31%  101M 10s
    ##  42400K .......... .......... .......... .......... .......... 31% 24.0M 10s
    ##  42450K .......... .......... .......... .......... .......... 31% 33.9M 10s
    ##  42500K .......... .......... .......... .......... .......... 31% 18.9M 10s
    ##  42550K .......... .......... .......... .......... .......... 31% 50.7M 10s
    ##  42600K .......... .......... .......... .......... .......... 31% 98.2M 10s
    ##  42650K .......... .......... .......... .......... .......... 31% 16.7M 10s
    ##  42700K .......... .......... .......... .......... .......... 31% 64.5M 10s
    ##  42750K .......... .......... .......... .......... .......... 31% 20.5M 10s
    ##  42800K .......... .......... .......... .......... .......... 31% 81.3M 10s
    ##  42850K .......... .......... .......... .......... .......... 31%  107M 10s
    ##  42900K .......... .......... .......... .......... .......... 31% 11.2M 10s
    ##  42950K .......... .......... .......... .......... .......... 31%  109M 10s
    ##  43000K .......... .......... .......... .......... .......... 31% 15.5M 10s
    ##  43050K .......... .......... .......... .......... .......... 31% 75.2M 10s
    ##  43100K .......... .......... .......... .......... .......... 32%  101M 10s
    ##  43150K .......... .......... .......... .......... .......... 32% 17.7M 10s
    ##  43200K .......... .......... .......... .......... .......... 32% 67.5M 10s
    ##  43250K .......... .......... .......... .......... .......... 32% 93.5M 10s
    ##  43300K .......... .......... .......... .......... .......... 32% 17.1M 10s
    ##  43350K .......... .......... .......... .......... .......... 32% 37.9M 10s
    ##  43400K .......... .......... .......... .......... .......... 32% 79.4M 10s
    ##  43450K .......... .......... .......... .......... .......... 32% 14.4M 10s
    ##  43500K .......... .......... .......... .......... .......... 32% 20.2M 10s
    ##  43550K .......... .......... .......... .......... .......... 32% 59.8M 10s
    ##  43600K .......... .......... .......... .......... .......... 32% 28.3M 10s
    ##  43650K .......... .......... .......... .......... .......... 32% 27.8M 10s
    ##  43700K .......... .......... .......... .......... .......... 32% 18.5M 10s
    ##  43750K .......... .......... .......... .......... .......... 32%  111M 10s
    ##  43800K .......... .......... .......... .......... .......... 32% 72.2M 10s
    ##  43850K .......... .......... .......... .......... .......... 32% 15.9M 10s
    ##  43900K .......... .......... .......... .......... .......... 32% 78.7M 10s
    ##  43950K .......... .......... .......... .......... .......... 32%  106M 10s
    ##  44000K .......... .......... .......... .......... .......... 32% 12.6M 10s
    ##  44050K .......... .......... .......... .......... .......... 32% 99.1M 10s
    ##  44100K .......... .......... .......... .......... .......... 32% 84.1M 10s
    ##  44150K .......... .......... .......... .......... .......... 32% 24.2M 10s
    ##  44200K .......... .......... .......... .......... .......... 32% 36.7M 10s
    ##  44250K .......... .......... .......... .......... .......... 32% 94.8M 10s
    ##  44300K .......... .......... .......... .......... .......... 32% 22.7M 10s
    ##  44350K .......... .......... .......... .......... .......... 32% 46.0M 10s
    ##  44400K .......... .......... .......... .......... .......... 32% 16.8M 10s
    ##  44450K .......... .......... .......... .......... .......... 33%  102M 10s
    ##  44500K .......... .......... .......... .......... .......... 33% 33.1M 10s
    ##  44550K .......... .......... .......... .......... .......... 33% 14.1M 10s
    ##  44600K .......... .......... .......... .......... .......... 33% 50.8M 10s
    ##  44650K .......... .......... .......... .......... .......... 33% 60.3M 10s
    ##  44700K .......... .......... .......... .......... .......... 33% 18.7M 10s
    ##  44750K .......... .......... .......... .......... .......... 33% 54.6M 10s
    ##  44800K .......... .......... .......... .......... .......... 33% 36.9M 10s
    ##  44850K .......... .......... .......... .......... .......... 33% 25.1M 10s
    ##  44900K .......... .......... .......... .......... .......... 33% 27.0M 9s
    ##  44950K .......... .......... .......... .......... .......... 33% 62.6M 9s
    ##  45000K .......... .......... .......... .......... .......... 33% 28.9M 9s
    ##  45050K .......... .......... .......... .......... .......... 33% 35.3M 9s
    ##  45100K .......... .......... .......... .......... .......... 33% 28.8M 9s
    ##  45150K .......... .......... .......... .......... .......... 33% 53.6M 9s
    ##  45200K .......... .......... .......... .......... .......... 33% 34.5M 9s
    ##  45250K .......... .......... .......... .......... .......... 33% 24.1M 9s
    ##  45300K .......... .......... .......... .......... .......... 33% 45.0M 9s
    ##  45350K .......... .......... .......... .......... .......... 33% 13.3M 9s
    ##  45400K .......... .......... .......... .......... .......... 33% 34.2M 9s
    ##  45450K .......... .......... .......... .......... .......... 33% 62.1M 9s
    ##  45500K .......... .......... .......... .......... .......... 33% 25.8M 9s
    ##  45550K .......... .......... .......... .......... .......... 33% 38.9M 9s
    ##  45600K .......... .......... .......... .......... .......... 33% 20.8M 9s
    ##  45650K .......... .......... .......... .......... .......... 33% 46.3M 9s
    ##  45700K .......... .......... .......... .......... .......... 33% 38.5M 9s
    ##  45750K .......... .......... .......... .......... .......... 33% 27.0M 9s
    ##  45800K .......... .......... .......... .......... .......... 34% 45.8M 9s
    ##  45850K .......... .......... .......... .......... .......... 34% 19.5M 9s
    ##  45900K .......... .......... .......... .......... .......... 34% 57.2M 9s
    ##  45950K .......... .......... .......... .......... .......... 34% 13.9M 9s
    ##  46000K .......... .......... .......... .......... .......... 34% 20.9M 9s
    ##  46050K .......... .......... .......... .......... .......... 34% 56.7M 9s
    ##  46100K .......... .......... .......... .......... .......... 34% 4.48M 9s
    ##  46150K .......... .......... .......... .......... .......... 34% 51.5M 9s
    ##  46200K .......... .......... .......... .......... .......... 34% 61.7M 9s
    ##  46250K .......... .......... .......... .......... .......... 34% 26.1M 9s
    ##  46300K .......... .......... .......... .......... .......... 34% 10.3M 9s
    ##  46350K .......... .......... .......... .......... .......... 34% 68.2M 9s
    ##  46400K .......... .......... .......... .......... .......... 34% 57.8M 9s
    ##  46450K .......... .......... .......... .......... .......... 34% 28.7M 9s
    ##  46500K .......... .......... .......... .......... .......... 34% 40.2M 9s
    ##  46550K .......... .......... .......... .......... .......... 34% 16.2M 9s
    ##  46600K .......... .......... .......... .......... .......... 34% 54.3M 9s
    ##  46650K .......... .......... .......... .......... .......... 34% 73.5M 9s
    ##  46700K .......... .......... .......... .......... .......... 34% 22.3M 9s
    ##  46750K .......... .......... .......... .......... .......... 34% 37.9M 9s
    ##  46800K .......... .......... .......... .......... .......... 34% 63.1M 9s
    ##  46850K .......... .......... .......... .......... .......... 34% 11.9M 9s
    ##  46900K .......... .......... .......... .......... .......... 34% 19.5M 9s
    ##  46950K .......... .......... .......... .......... .......... 34% 64.8M 9s
    ##  47000K .......... .......... .......... .......... .......... 34% 53.4M 9s
    ##  47050K .......... .......... .......... .......... .......... 34% 28.3M 9s
    ##  47100K .......... .......... .......... .......... .......... 34% 52.7M 9s
    ##  47150K .......... .......... .......... .......... .......... 35% 63.7M 9s
    ##  47200K .......... .......... .......... .......... .......... 35% 21.3M 9s
    ##  47250K .......... .......... .......... .......... .......... 35% 66.9M 9s
    ##  47300K .......... .......... .......... .......... .......... 35% 69.8M 9s
    ##  47350K .......... .......... .......... .......... .......... 35% 22.4M 9s
    ##  47400K .......... .......... .......... .......... .......... 35% 54.4M 9s
    ##  47450K .......... .......... .......... .......... .......... 35% 26.5M 9s
    ##  47500K .......... .......... .......... .......... .......... 35% 30.1M 9s
    ##  47550K .......... .......... .......... .......... .......... 35% 78.9M 9s
    ##  47600K .......... .......... .......... .......... .......... 35% 30.2M 9s
    ##  47650K .......... .......... .......... .......... .......... 35% 28.8M 9s
    ##  47700K .......... .......... .......... .......... .......... 35% 34.1M 9s
    ##  47750K .......... .......... .......... .......... .......... 35% 64.4M 9s
    ##  47800K .......... .......... .......... .......... .......... 35% 22.5M 9s
    ##  47850K .......... .......... .......... .......... .......... 35% 61.0M 9s
    ##  47900K .......... .......... .......... .......... .......... 35% 18.9M 9s
    ##  47950K .......... .......... .......... .......... .......... 35% 54.6M 9s
    ##  48000K .......... .......... .......... .......... .......... 35% 57.1M 9s
    ##  48050K .......... .......... .......... .......... .......... 35% 72.7M 9s
    ##  48100K .......... .......... .......... .......... .......... 35% 19.2M 9s
    ##  48150K .......... .......... .......... .......... .......... 35% 77.6M 9s
    ##  48200K .......... .......... .......... .......... .......... 35% 19.8M 9s
    ##  48250K .......... .......... .......... .......... .......... 35% 79.9M 9s
    ##  48300K .......... .......... .......... .......... .......... 35% 19.2M 9s
    ##  48350K .......... .......... .......... .......... .......... 35% 20.4M 9s
    ##  48400K .......... .......... .......... .......... .......... 35% 60.5M 9s
    ##  48450K .......... .......... .......... .......... .......... 35% 78.2M 9s
    ##  48500K .......... .......... .......... .......... .......... 36% 22.5M 9s
    ##  48550K .......... .......... .......... .......... .......... 36% 71.1M 9s
    ##  48600K .......... .......... .......... .......... .......... 36% 14.0M 9s
    ##  48650K .......... .......... .......... .......... .......... 36% 74.5M 9s
    ##  48700K .......... .......... .......... .......... .......... 36% 8.20M 9s
    ##  48750K .......... .......... .......... .......... .......... 36% 92.3M 9s
    ##  48800K .......... .......... .......... .......... .......... 36% 77.1M 9s
    ##  48850K .......... .......... .......... .......... .......... 36% 17.8M 9s
    ##  48900K .......... .......... .......... .......... .......... 36% 65.7M 9s
    ##  48950K .......... .......... .......... .......... .......... 36% 67.4M 9s
    ##  49000K .......... .......... .......... .......... .......... 36% 6.35M 9s
    ##  49050K .......... .......... .......... .......... .......... 36% 83.6M 9s
    ##  49100K .......... .......... .......... .......... .......... 36% 80.2M 9s
    ##  49150K .......... .......... .......... .......... .......... 36% 21.6M 9s
    ##  49200K .......... .......... .......... .......... .......... 36% 21.6M 9s
    ##  49250K .......... .......... .......... .......... .......... 36% 77.0M 8s
    ##  49300K .......... .......... .......... .......... .......... 36% 29.6M 8s
    ##  49350K .......... .......... .......... .......... .......... 36% 20.3M 8s
    ##  49400K .......... .......... .......... .......... .......... 36% 67.1M 8s
    ##  49450K .......... .......... .......... .......... .......... 36% 97.3M 8s
    ##  49500K .......... .......... .......... .......... .......... 36% 18.1M 8s
    ##  49550K .......... .......... .......... .......... .......... 36% 46.4M 8s
    ##  49600K .......... .......... .......... .......... .......... 36% 74.4M 8s
    ##  49650K .......... .......... .......... .......... .......... 36% 18.0M 8s
    ##  49700K .......... .......... .......... .......... .......... 36% 65.4M 8s
    ##  49750K .......... .......... .......... .......... .......... 36% 83.5M 8s
    ##  49800K .......... .......... .......... .......... .......... 36% 81.5M 8s
    ##  49850K .......... .......... .......... .......... .......... 37% 15.9M 8s
    ##  49900K .......... .......... .......... .......... .......... 37% 74.0M 8s
    ##  49950K .......... .......... .......... .......... .......... 37% 98.9M 8s
    ##  50000K .......... .......... .......... .......... .......... 37% 17.9M 8s
    ##  50050K .......... .......... .......... .......... .......... 37% 85.4M 8s
    ##  50100K .......... .......... .......... .......... .......... 37% 20.5M 8s
    ##  50150K .......... .......... .......... .......... .......... 37% 67.6M 8s
    ##  50200K .......... .......... .......... .......... .......... 37% 80.5M 8s
    ##  50250K .......... .......... .......... .......... .......... 37% 18.5M 8s
    ##  50300K .......... .......... .......... .......... .......... 37% 76.2M 8s
    ##  50350K .......... .......... .......... .......... .......... 37% 76.0M 8s
    ##  50400K .......... .......... .......... .......... .......... 37% 18.8M 8s
    ##  50450K .......... .......... .......... .......... .......... 37% 78.8M 8s
    ##  50500K .......... .......... .......... .......... .......... 37% 28.1M 8s
    ##  50550K .......... .......... .......... .......... .......... 37% 73.8M 8s
    ##  50600K .......... .......... .......... .......... .......... 37% 50.9M 8s
    ##  50650K .......... .......... .......... .......... .......... 37% 22.9M 8s
    ##  50700K .......... .......... .......... .......... .......... 37% 62.9M 8s
    ##  50750K .......... .......... .......... .......... .......... 37% 24.4M 8s
    ##  50800K .......... .......... .......... .......... .......... 37% 45.9M 8s
    ##  50850K .......... .......... .......... .......... .......... 37% 64.8M 8s
    ##  50900K .......... .......... .......... .......... .......... 37% 19.1M 8s
    ##  50950K .......... .......... .......... .......... .......... 37% 88.3M 8s
    ##  51000K .......... .......... .......... .......... .......... 37% 18.0M 8s
    ##  51050K .......... .......... .......... .......... .......... 37% 81.0M 8s
    ##  51100K .......... .......... .......... .......... .......... 37% 63.6M 8s
    ##  51150K .......... .......... .......... .......... .......... 37% 94.4M 8s
    ##  51200K .......... .......... .......... .......... .......... 38% 28.8M 8s
    ##  51250K .......... .......... .......... .......... .......... 38% 33.1M 8s
    ##  51300K .......... .......... .......... .......... .......... 38% 78.0M 8s
    ##  51350K .......... .......... .......... .......... .......... 38% 82.2M 8s
    ##  51400K .......... .......... .......... .......... .......... 38% 59.8M 8s
    ##  51450K .......... .......... .......... .......... .......... 38% 27.9M 8s
    ##  51500K .......... .......... .......... .......... .......... 38% 55.5M 8s
    ##  51550K .......... .......... .......... .......... .......... 38% 94.5M 8s
    ##  51600K .......... .......... .......... .......... .......... 38% 23.3M 8s
    ##  51650K .......... .......... .......... .......... .......... 38% 13.1M 8s
    ##  51700K .......... .......... .......... .......... .......... 38% 45.0M 8s
    ##  51750K .......... .......... .......... .......... .......... 38% 88.2M 8s
    ##  51800K .......... .......... .......... .......... .......... 38% 23.0M 8s
    ##  51850K .......... .......... .......... .......... .......... 38% 6.99M 8s
    ##  51900K .......... .......... .......... .......... .......... 38% 69.4M 8s
    ##  51950K .......... .......... .......... .......... .......... 38% 62.8M 8s
    ##  52000K .......... .......... .......... .......... .......... 38% 71.3M 8s
    ##  52050K .......... .......... .......... .......... .......... 38% 13.1M 8s
    ##  52100K .......... .......... .......... .......... .......... 38% 68.1M 8s
    ##  52150K .......... .......... .......... .......... .......... 38% 82.7M 8s
    ##  52200K .......... .......... .......... .......... .......... 38% 21.1M 8s
    ##  52250K .......... .......... .......... .......... .......... 38% 87.2M 8s
    ##  52300K .......... .......... .......... .......... .......... 38% 44.9M 8s
    ##  52350K .......... .......... .......... .......... .......... 38% 92.3M 8s
    ##  52400K .......... .......... .......... .......... .......... 38% 31.8M 8s
    ##  52450K .......... .......... .......... .......... .......... 38% 95.3M 8s
    ##  52500K .......... .......... .......... .......... .......... 39% 21.0M 8s
    ##  52550K .......... .......... .......... .......... .......... 39% 60.2M 8s
    ##  52600K .......... .......... .......... .......... .......... 39% 27.2M 8s
    ##  52650K .......... .......... .......... .......... .......... 39% 6.83M 8s
    ##  52700K .......... .......... .......... .......... .......... 39% 59.3M 8s
    ##  52750K .......... .......... .......... .......... .......... 39% 84.9M 8s
    ##  52800K .......... .......... .......... .......... .......... 39% 77.4M 8s
    ##  52850K .......... .......... .......... .......... .......... 39% 10.2M 8s
    ##  52900K .......... .......... .......... .......... .......... 39% 67.6M 8s
    ##  52950K .......... .......... .......... .......... .......... 39% 74.0M 8s
    ##  53000K .......... .......... .......... .......... .......... 39% 74.1M 8s
    ##  53050K .......... .......... .......... .......... .......... 39% 20.8M 8s
    ##  53100K .......... .......... .......... .......... .......... 39% 44.8M 8s
    ##  53150K .......... .......... .......... .......... .......... 39% 74.1M 8s
    ##  53200K .......... .......... .......... .......... .......... 39% 80.7M 8s
    ##  53250K .......... .......... .......... .......... .......... 39% 31.4M 8s
    ##  53300K .......... .......... .......... .......... .......... 39% 26.8M 8s
    ##  53350K .......... .......... .......... .......... .......... 39% 62.4M 8s
    ##  53400K .......... .......... .......... .......... .......... 39% 54.3M 8s
    ##  53450K .......... .......... .......... .......... .......... 39% 84.4M 8s
    ##  53500K .......... .......... .......... .......... .......... 39% 21.5M 8s
    ##  53550K .......... .......... .......... .......... .......... 39% 80.5M 8s
    ##  53600K .......... .......... .......... .......... .......... 39% 76.6M 8s
    ##  53650K .......... .......... .......... .......... .......... 39% 15.5M 8s
    ##  53700K .......... .......... .......... .......... .......... 39% 39.9M 8s
    ##  53750K .......... .......... .......... .......... .......... 39% 34.9M 8s
    ##  53800K .......... .......... .......... .......... .......... 39% 38.4M 8s
    ##  53850K .......... .......... .......... .......... .......... 40% 73.7M 8s
    ##  53900K .......... .......... .......... .......... .......... 40% 77.7M 8s
    ##  53950K .......... .......... .......... .......... .......... 40% 25.5M 8s
    ##  54000K .......... .......... .......... .......... .......... 40% 22.8M 8s
    ##  54050K .......... .......... .......... .......... .......... 40% 54.6M 8s
    ##  54100K .......... .......... .......... .......... .......... 40% 52.7M 8s
    ##  54150K .......... .......... .......... .......... .......... 40% 69.8M 7s
    ##  54200K .......... .......... .......... .......... .......... 40% 11.5M 7s
    ##  54250K .......... .......... .......... .......... .......... 40% 68.9M 7s
    ##  54300K .......... .......... .......... .......... .......... 40% 76.8M 7s
    ##  54350K .......... .......... .......... .......... .......... 40% 49.1M 7s
    ##  54400K .......... .......... .......... .......... .......... 40% 28.3M 7s
    ##  54450K .......... .......... .......... .......... .......... 40% 58.4M 7s
    ##  54500K .......... .......... .......... .......... .......... 40% 56.2M 7s
    ##  54550K .......... .......... .......... .......... .......... 40% 20.3M 7s
    ##  54600K .......... .......... .......... .......... .......... 40% 31.4M 7s
    ##  54650K .......... .......... .......... .......... .......... 40% 35.9M 7s
    ##  54700K .......... .......... .......... .......... .......... 40% 15.9M 7s
    ##  54750K .......... .......... .......... .......... .......... 40% 8.80M 7s
    ##  54800K .......... .......... .......... .......... .......... 40% 85.5M 7s
    ##  54850K .......... .......... .......... .......... .......... 40% 86.5M 7s
    ##  54900K .......... .......... .......... .......... .......... 40% 82.8M 7s
    ##  54950K .......... .......... .......... .......... .......... 40% 18.7M 7s
    ##  55000K .......... .......... .......... .......... .......... 40% 63.9M 7s
    ##  55050K .......... .......... .......... .......... .......... 40% 90.7M 7s
    ##  55100K .......... .......... .......... .......... .......... 40% 87.4M 7s
    ##  55150K .......... .......... .......... .......... .......... 40% 12.0M 7s
    ##  55200K .......... .......... .......... .......... .......... 41% 63.5M 7s
    ##  55250K .......... .......... .......... .......... .......... 41% 70.3M 7s
    ##  55300K .......... .......... .......... .......... .......... 41% 90.0M 7s
    ##  55350K .......... .......... .......... .......... .......... 41% 24.2M 7s
    ##  55400K .......... .......... .......... .......... .......... 41% 26.6M 7s
    ##  55450K .......... .......... .......... .......... .......... 41% 80.7M 7s
    ##  55500K .......... .......... .......... .......... .......... 41% 56.6M 7s
    ##  55550K .......... .......... .......... .......... .......... 41% 18.9M 7s
    ##  55600K .......... .......... .......... .......... .......... 41% 70.1M 7s
    ##  55650K .......... .......... .......... .......... .......... 41% 62.8M 7s
    ##  55700K .......... .......... .......... .......... .......... 41% 77.2M 7s
    ##  55750K .......... .......... .......... .......... .......... 41% 28.8M 7s
    ##  55800K .......... .......... .......... .......... .......... 41% 68.7M 7s
    ##  55850K .......... .......... .......... .......... .......... 41% 35.0M 7s
    ##  55900K .......... .......... .......... .......... .......... 41% 80.1M 7s
    ##  55950K .......... .......... .......... .......... .......... 41% 63.3M 7s
    ##  56000K .......... .......... .......... .......... .......... 41% 26.3M 7s
    ##  56050K .......... .......... .......... .......... .......... 41% 71.7M 7s
    ##  56100K .......... .......... .......... .......... .......... 41% 61.7M 7s
    ##  56150K .......... .......... .......... .......... .......... 41% 12.6M 7s
    ##  56200K .......... .......... .......... .......... .......... 41% 10.5M 7s
    ##  56250K .......... .......... .......... .......... .......... 41% 47.1M 7s
    ##  56300K .......... .......... .......... .......... .......... 41% 67.0M 7s
    ##  56350K .......... .......... .......... .......... .......... 41% 85.6M 7s
    ##  56400K .......... .......... .......... .......... .......... 41% 20.7M 7s
    ##  56450K .......... .......... .......... .......... .......... 41% 69.3M 7s
    ##  56500K .......... .......... .......... .......... .......... 41% 74.6M 7s
    ##  56550K .......... .......... .......... .......... .......... 42% 73.2M 7s
    ##  56600K .......... .......... .......... .......... .......... 42% 23.2M 7s
    ##  56650K .......... .......... .......... .......... .......... 42% 70.4M 7s
    ##  56700K .......... .......... .......... .......... .......... 42% 65.9M 7s
    ##  56750K .......... .......... .......... .......... .......... 42% 77.0M 7s
    ##  56800K .......... .......... .......... .......... .......... 42% 18.7M 7s
    ##  56850K .......... .......... .......... .......... .......... 42% 59.5M 7s
    ##  56900K .......... .......... .......... .......... .......... 42% 45.9M 7s
    ##  56950K .......... .......... .......... .......... .......... 42% 93.3M 7s
    ##  57000K .......... .......... .......... .......... .......... 42% 24.8M 7s
    ##  57050K .......... .......... .......... .......... .......... 42% 49.4M 7s
    ##  57100K .......... .......... .......... .......... .......... 42% 73.8M 7s
    ##  57150K .......... .......... .......... .......... .......... 42% 13.0M 7s
    ##  57200K .......... .......... .......... .......... .......... 42% 4.67M 7s
    ##  57250K .......... .......... .......... .......... .......... 42% 77.8M 7s
    ##  57300K .......... .......... .......... .......... .......... 42% 83.1M 7s
    ##  57350K .......... .......... .......... .......... .......... 42% 78.5M 7s
    ##  57400K .......... .......... .......... .......... .......... 42% 8.99M 7s
    ##  57450K .......... .......... .......... .......... .......... 42% 69.0M 7s
    ##  57500K .......... .......... .......... .......... .......... 42% 70.5M 7s
    ##  57550K .......... .......... .......... .......... .......... 42% 94.4M 7s
    ##  57600K .......... .......... .......... .......... .......... 42% 5.87M 7s
    ##  57650K .......... .......... .......... .......... .......... 42% 54.9M 7s
    ##  57700K .......... .......... .......... .......... .......... 42% 74.1M 7s
    ##  57750K .......... .......... .......... .......... .......... 42% 86.6M 7s
    ##  57800K .......... .......... .......... .......... .......... 42% 21.2M 7s
    ##  57850K .......... .......... .......... .......... .......... 42% 55.6M 7s
    ##  57900K .......... .......... .......... .......... .......... 43% 80.2M 7s
    ##  57950K .......... .......... .......... .......... .......... 43% 86.2M 7s
    ##  58000K .......... .......... .......... .......... .......... 43% 13.3M 7s
    ##  58050K .......... .......... .......... .......... .......... 43% 57.2M 7s
    ##  58100K .......... .......... .......... .......... .......... 43% 77.1M 7s
    ##  58150K .......... .......... .......... .......... .......... 43% 80.4M 7s
    ##  58200K .......... .......... .......... .......... .......... 43% 12.9M 7s
    ##  58250K .......... .......... .......... .......... .......... 43% 74.4M 7s
    ##  58300K .......... .......... .......... .......... .......... 43% 59.5M 7s
    ##  58350K .......... .......... .......... .......... .......... 43% 94.1M 7s
    ##  58400K .......... .......... .......... .......... .......... 43% 72.9M 7s
    ##  58450K .......... .......... .......... .......... .......... 43% 20.4M 7s
    ##  58500K .......... .......... .......... .......... .......... 43% 45.8M 7s
    ##  58550K .......... .......... .......... .......... .......... 43% 86.3M 7s
    ##  58600K .......... .......... .......... .......... .......... 43% 78.8M 7s
    ##  58650K .......... .......... .......... .......... .......... 43% 42.9M 7s
    ##  58700K .......... .......... .......... .......... .......... 43% 50.1M 7s
    ##  58750K .......... .......... .......... .......... .......... 43% 38.0M 7s
    ##  58800K .......... .......... .......... .......... .......... 43% 77.2M 7s
    ##  58850K .......... .......... .......... .......... .......... 43% 89.8M 7s
    ##  58900K .......... .......... .......... .......... .......... 43% 24.5M 7s
    ##  58950K .......... .......... .......... .......... .......... 43% 69.8M 7s
    ##  59000K .......... .......... .......... .......... .......... 43% 63.6M 7s
    ##  59050K .......... .......... .......... .......... .......... 43% 26.0M 7s
    ##  59100K .......... .......... .......... .......... .......... 43% 68.2M 7s
    ##  59150K .......... .......... .......... .......... .......... 43% 64.4M 7s
    ##  59200K .......... .......... .......... .......... .......... 43% 55.0M 7s
    ##  59250K .......... .......... .......... .......... .......... 44% 36.0M 7s
    ##  59300K .......... .......... .......... .......... .......... 44% 69.2M 7s
    ##  59350K .......... .......... .......... .......... .......... 44% 62.6M 7s
    ##  59400K .......... .......... .......... .......... .......... 44% 49.6M 7s
    ##  59450K .......... .......... .......... .......... .......... 44% 35.0M 7s
    ##  59500K .......... .......... .......... .......... .......... 44% 5.43M 7s
    ##  59550K .......... .......... .......... .......... .......... 44% 22.5M 7s
    ##  59600K .......... .......... .......... .......... .......... 44% 69.3M 7s
    ##  59650K .......... .......... .......... .......... .......... 44% 58.0M 7s
    ##  59700K .......... .......... .......... .......... .......... 44% 14.4M 7s
    ##  59750K .......... .......... .......... .......... .......... 44% 48.5M 7s
    ##  59800K .......... .......... .......... .......... .......... 44% 60.5M 7s
    ##  59850K .......... .......... .......... .......... .......... 44% 83.9M 7s
    ##  59900K .......... .......... .......... .......... .......... 44% 34.5M 7s
    ##  59950K .......... .......... .......... .......... .......... 44% 57.2M 6s
    ##  60000K .......... .......... .......... .......... .......... 44% 75.4M 6s
    ##  60050K .......... .......... .......... .......... .......... 44% 68.3M 6s
    ##  60100K .......... .......... .......... .......... .......... 44% 29.9M 6s
    ##  60150K .......... .......... .......... .......... .......... 44% 78.0M 6s
    ##  60200K .......... .......... .......... .......... .......... 44% 39.1M 6s
    ##  60250K .......... .......... .......... .......... .......... 44% 7.09M 6s
    ##  60300K .......... .......... .......... .......... .......... 44% 33.6M 6s
    ##  60350K .......... .......... .......... .......... .......... 44% 39.9M 6s
    ##  60400K .......... .......... .......... .......... .......... 44% 50.4M 6s
    ##  60450K .......... .......... .......... .......... .......... 44% 53.1M 6s
    ##  60500K .......... .......... .......... .......... .......... 44% 40.0M 6s
    ##  60550K .......... .......... .......... .......... .......... 44% 58.4M 6s
    ##  60600K .......... .......... .......... .......... .......... 45% 50.9M 6s
    ##  60650K .......... .......... .......... .......... .......... 45% 48.1M 6s
    ##  60700K .......... .......... .......... .......... .......... 45% 46.6M 6s
    ##  60750K .......... .......... .......... .......... .......... 45% 54.3M 6s
    ##  60800K .......... .......... .......... .......... .......... 45% 47.2M 6s
    ##  60850K .......... .......... .......... .......... .......... 45% 56.9M 6s
    ##  60900K .......... .......... .......... .......... .......... 45% 41.6M 6s
    ##  60950K .......... .......... .......... .......... .......... 45% 45.7M 6s
    ##  61000K .......... .......... .......... .......... .......... 45% 46.7M 6s
    ##  61050K .......... .......... .......... .......... .......... 45% 59.3M 6s
    ##  61100K .......... .......... .......... .......... .......... 45% 14.4M 6s
    ##  61150K .......... .......... .......... .......... .......... 45% 53.5M 6s
    ##  61200K .......... .......... .......... .......... .......... 45% 51.2M 6s
    ##  61250K .......... .......... .......... .......... .......... 45% 52.2M 6s
    ##  61300K .......... .......... .......... .......... .......... 45% 30.9M 6s
    ##  61350K .......... .......... .......... .......... .......... 45% 43.8M 6s
    ##  61400K .......... .......... .......... .......... .......... 45% 39.7M 6s
    ##  61450K .......... .......... .......... .......... .......... 45% 53.1M 6s
    ##  61500K .......... .......... .......... .......... .......... 45% 20.5M 6s
    ##  61550K .......... .......... .......... .......... .......... 45% 44.0M 6s
    ##  61600K .......... .......... .......... .......... .......... 45% 52.9M 6s
    ##  61650K .......... .......... .......... .......... .......... 45% 54.3M 6s
    ##  61700K .......... .......... .......... .......... .......... 45% 48.1M 6s
    ##  61750K .......... .......... .......... .......... .......... 45% 46.2M 6s
    ##  61800K .......... .......... .......... .......... .......... 45% 47.7M 6s
    ##  61850K .......... .......... .......... .......... .......... 45% 41.9M 6s
    ##  61900K .......... .......... .......... .......... .......... 45% 52.8M 6s
    ##  61950K .......... .......... .......... .......... .......... 46% 48.8M 6s
    ##  62000K .......... .......... .......... .......... .......... 46% 49.9M 6s
    ##  62050K .......... .......... .......... .......... .......... 46% 27.2M 6s
    ##  62100K .......... .......... .......... .......... .......... 46% 12.4M 6s
    ##  62150K .......... .......... .......... .......... .......... 46% 39.3M 6s
    ##  62200K .......... .......... .......... .......... .......... 46% 57.1M 6s
    ##  62250K .......... .......... .......... .......... .......... 46% 58.4M 6s
    ##  62300K .......... .......... .......... .......... .......... 46% 11.4M 6s
    ##  62350K .......... .......... .......... .......... .......... 46% 47.8M 6s
    ##  62400K .......... .......... .......... .......... .......... 46% 50.8M 6s
    ##  62450K .......... .......... .......... .......... .......... 46% 53.3M 6s
    ##  62500K .......... .......... .......... .......... .......... 46% 38.0M 6s
    ##  62550K .......... .......... .......... .......... .......... 46% 48.9M 6s
    ##  62600K .......... .......... .......... .......... .......... 46% 50.7M 6s
    ##  62650K .......... .......... .......... .......... .......... 46% 55.6M 6s
    ##  62700K .......... .......... .......... .......... .......... 46% 27.9M 6s
    ##  62750K .......... .......... .......... .......... .......... 46% 61.3M 6s
    ##  62800K .......... .......... .......... .......... .......... 46% 55.1M 6s
    ##  62850K .......... .......... .......... .......... .......... 46% 66.3M 6s
    ##  62900K .......... .......... .......... .......... .......... 46% 18.7M 6s
    ##  62950K .......... .......... .......... .......... .......... 46% 25.5M 6s
    ##  63000K .......... .......... .......... .......... .......... 46% 32.9M 6s
    ##  63050K .......... .......... .......... .......... .......... 46% 27.5M 6s
    ##  63100K .......... .......... .......... .......... .......... 46% 52.3M 6s
    ##  63150K .......... .......... .......... .......... .......... 46% 45.1M 6s
    ##  63200K .......... .......... .......... .......... .......... 46% 54.1M 6s
    ##  63250K .......... .......... .......... .......... .......... 46% 52.9M 6s
    ##  63300K .......... .......... .......... .......... .......... 47% 39.5M 6s
    ##  63350K .......... .......... .......... .......... .......... 47% 25.9M 6s
    ##  63400K .......... .......... .......... .......... .......... 47% 54.0M 6s
    ##  63450K .......... .......... .......... .......... .......... 47% 67.1M 6s
    ##  63500K .......... .......... .......... .......... .......... 47% 51.9M 6s
    ##  63550K .......... .......... .......... .......... .......... 47% 51.6M 6s
    ##  63600K .......... .......... .......... .......... .......... 47% 53.3M 6s
    ##  63650K .......... .......... .......... .......... .......... 47% 67.4M 6s
    ##  63700K .......... .......... .......... .......... .......... 47% 31.4M 6s
    ##  63750K .......... .......... .......... .......... .......... 47% 55.4M 6s
    ##  63800K .......... .......... .......... .......... .......... 47% 41.6M 6s
    ##  63850K .......... .......... .......... .......... .......... 47% 60.8M 6s
    ##  63900K .......... .......... .......... .......... .......... 47% 29.3M 6s
    ##  63950K .......... .......... .......... .......... .......... 47% 67.1M 6s
    ##  64000K .......... .......... .......... .......... .......... 47% 54.7M 6s
    ##  64050K .......... .......... .......... .......... .......... 47% 69.2M 6s
    ##  64100K .......... .......... .......... .......... .......... 47% 24.5M 6s
    ##  64150K .......... .......... .......... .......... .......... 47% 53.1M 6s
    ##  64200K .......... .......... .......... .......... .......... 47% 45.4M 6s
    ##  64250K .......... .......... .......... .......... .......... 47% 59.5M 6s
    ##  64300K .......... .......... .......... .......... .......... 47% 57.0M 6s
    ##  64350K .......... .......... .......... .......... .......... 47% 18.2M 6s
    ##  64400K .......... .......... .......... .......... .......... 47% 45.2M 6s
    ##  64450K .......... .......... .......... .......... .......... 47% 67.9M 6s
    ##  64500K .......... .......... .......... .......... .......... 47% 49.9M 6s
    ##  64550K .......... .......... .......... .......... .......... 47% 59.5M 6s
    ##  64600K .......... .......... .......... .......... .......... 47% 55.6M 6s
    ##  64650K .......... .......... .......... .......... .......... 48% 49.6M 6s
    ##  64700K .......... .......... .......... .......... .......... 48% 42.1M 6s
    ##  64750K .......... .......... .......... .......... .......... 48% 52.2M 6s
    ##  64800K .......... .......... .......... .......... .......... 48% 41.4M 6s
    ##  64850K .......... .......... .......... .......... .......... 48% 59.8M 6s
    ##  64900K .......... .......... .......... .......... .......... 48% 43.1M 6s
    ##  64950K .......... .......... .......... .......... .......... 48% 46.6M 6s
    ##  65000K .......... .......... .......... .......... .......... 48% 52.8M 6s
    ##  65050K .......... .......... .......... .......... .......... 48% 69.4M 6s
    ##  65100K .......... .......... .......... .......... .......... 48% 49.2M 6s
    ##  65150K .......... .......... .......... .......... .......... 48% 56.1M 6s
    ##  65200K .......... .......... .......... .......... .......... 48% 64.4M 6s
    ##  65250K .......... .......... .......... .......... .......... 48% 30.1M 6s
    ##  65300K .......... .......... .......... .......... .......... 48% 54.7M 6s
    ##  65350K .......... .......... .......... .......... .......... 48% 49.8M 6s
    ##  65400K .......... .......... .......... .......... .......... 48% 38.0M 6s
    ##  65450K .......... .......... .......... .......... .......... 48% 61.8M 6s
    ##  65500K .......... .......... .......... .......... .......... 48% 28.0M 6s
    ##  65550K .......... .......... .......... .......... .......... 48% 72.0M 6s
    ##  65600K .......... .......... .......... .......... .......... 48% 55.6M 6s
    ##  65650K .......... .......... .......... .......... .......... 48% 22.5M 6s
    ##  65700K .......... .......... .......... .......... .......... 48% 45.4M 6s
    ##  65750K .......... .......... .......... .......... .......... 48% 72.3M 6s
    ##  65800K .......... .......... .......... .......... .......... 48% 62.4M 6s
    ##  65850K .......... .......... .......... .......... .......... 48% 26.7M 6s
    ##  65900K .......... .......... .......... .......... .......... 48% 41.3M 6s
    ##  65950K .......... .......... .......... .......... .......... 48% 63.9M 6s
    ##  66000K .......... .......... .......... .......... .......... 49% 63.4M 6s
    ##  66050K .......... .......... .......... .......... .......... 49% 33.5M 6s
    ##  66100K .......... .......... .......... .......... .......... 49% 58.9M 6s
    ##  66150K .......... .......... .......... .......... .......... 49% 57.2M 6s
    ##  66200K .......... .......... .......... .......... .......... 49% 63.7M 6s
    ##  66250K .......... .......... .......... .......... .......... 49% 35.8M 6s
    ##  66300K .......... .......... .......... .......... .......... 49% 55.3M 6s
    ##  66350K .......... .......... .......... .......... .......... 49% 47.1M 6s
    ##  66400K .......... .......... .......... .......... .......... 49% 62.1M 6s
    ##  66450K .......... .......... .......... .......... .......... 49% 61.3M 6s
    ##  66500K .......... .......... .......... .......... .......... 49% 30.0M 6s
    ##  66550K .......... .......... .......... .......... .......... 49% 56.1M 6s
    ##  66600K .......... .......... .......... .......... .......... 49% 50.4M 5s
    ##  66650K .......... .......... .......... .......... .......... 49% 73.5M 5s
    ##  66700K .......... .......... .......... .......... .......... 49% 30.7M 5s
    ##  66750K .......... .......... .......... .......... .......... 49% 61.9M 5s
    ##  66800K .......... .......... .......... .......... .......... 49% 53.2M 5s
    ##  66850K .......... .......... .......... .......... .......... 49% 78.9M 5s
    ##  66900K .......... .......... .......... .......... .......... 49% 15.5M 5s
    ##  66950K .......... .......... .......... .......... .......... 49% 68.4M 5s
    ##  67000K .......... .......... .......... .......... .......... 49% 74.2M 5s
    ##  67050K .......... .......... .......... .......... .......... 49% 17.9M 5s
    ##  67100K .......... .......... .......... .......... .......... 49% 56.5M 5s
    ##  67150K .......... .......... .......... .......... .......... 49% 80.5M 5s
    ##  67200K .......... .......... .......... .......... .......... 49% 64.1M 5s
    ##  67250K .......... .......... .......... .......... .......... 49% 30.0M 5s
    ##  67300K .......... .......... .......... .......... .......... 49% 64.2M 5s
    ##  67350K .......... .......... .......... .......... .......... 50% 49.1M 5s
    ##  67400K .......... .......... .......... .......... .......... 50% 65.8M 5s
    ##  67450K .......... .......... .......... .......... .......... 50% 71.8M 5s
    ##  67500K .......... .......... .......... .......... .......... 50% 64.1M 5s
    ##  67550K .......... .......... .......... .......... .......... 50% 75.8M 5s
    ##  67600K .......... .......... .......... .......... .......... 50% 39.4M 5s
    ##  67650K .......... .......... .......... .......... .......... 50% 64.0M 5s
    ##  67700K .......... .......... .......... .......... .......... 50% 20.7M 5s
    ##  67750K .......... .......... .......... .......... .......... 50% 61.0M 5s
    ##  67800K .......... .......... .......... .......... .......... 50% 64.7M 5s
    ##  67850K .......... .......... .......... .......... .......... 50% 76.7M 5s
    ##  67900K .......... .......... .......... .......... .......... 50% 24.7M 5s
    ##  67950K .......... .......... .......... .......... .......... 50% 55.2M 5s
    ##  68000K .......... .......... .......... .......... .......... 50% 66.4M 5s
    ##  68050K .......... .......... .......... .......... .......... 50% 75.2M 5s
    ##  68100K .......... .......... .......... .......... .......... 50% 27.0M 5s
    ##  68150K .......... .......... .......... .......... .......... 50% 46.2M 5s
    ##  68200K .......... .......... .......... .......... .......... 50% 59.2M 5s
    ##  68250K .......... .......... .......... .......... .......... 50% 69.3M 5s
    ##  68300K .......... .......... .......... .......... .......... 50% 69.0M 5s
    ##  68350K .......... .......... .......... .......... .......... 50% 33.8M 5s
    ##  68400K .......... .......... .......... .......... .......... 50% 57.6M 5s
    ##  68450K .......... .......... .......... .......... .......... 50% 58.7M 5s
    ##  68500K .......... .......... .......... .......... .......... 50% 64.9M 5s
    ##  68550K .......... .......... .......... .......... .......... 50% 67.9M 5s
    ##  68600K .......... .......... .......... .......... .......... 50% 61.0M 5s
    ##  68650K .......... .......... .......... .......... .......... 50% 64.0M 5s
    ##  68700K .......... .......... .......... .......... .......... 51% 67.1M 5s
    ##  68750K .......... .......... .......... .......... .......... 51% 74.2M 5s
    ##  68800K .......... .......... .......... .......... .......... 51% 60.6M 5s
    ##  68850K .......... .......... .......... .......... .......... 51% 30.6M 5s
    ##  68900K .......... .......... .......... .......... .......... 51% 34.8M 5s
    ##  68950K .......... .......... .......... .......... .......... 51% 67.3M 5s
    ##  69000K .......... .......... .......... .......... .......... 51% 62.4M 5s
    ##  69050K .......... .......... .......... .......... .......... 51% 68.6M 5s
    ##  69100K .......... .......... .......... .......... .......... 51% 19.6M 5s
    ##  69150K .......... .......... .......... .......... .......... 51% 57.0M 5s
    ##  69200K .......... .......... .......... .......... .......... 51% 68.3M 5s
    ##  69250K .......... .......... .......... .......... .......... 51% 66.9M 5s
    ##  69300K .......... .......... .......... .......... .......... 51% 67.1M 5s
    ##  69350K .......... .......... .......... .......... .......... 51% 71.2M 5s
    ##  69400K .......... .......... .......... .......... .......... 51% 50.2M 5s
    ##  69450K .......... .......... .......... .......... .......... 51% 60.2M 5s
    ##  69500K .......... .......... .......... .......... .......... 51% 61.5M 5s
    ##  69550K .......... .......... .......... .......... .......... 51% 26.5M 5s
    ##  69600K .......... .......... .......... .......... .......... 51% 61.0M 5s
    ##  69650K .......... .......... .......... .......... .......... 51% 69.4M 5s
    ##  69700K .......... .......... .......... .......... .......... 51% 66.5M 5s
    ##  69750K .......... .......... .......... .......... .......... 51% 66.8M 5s
    ##  69800K .......... .......... .......... .......... .......... 51% 60.2M 5s
    ##  69850K .......... .......... .......... .......... .......... 51% 23.6M 5s
    ##  69900K .......... .......... .......... .......... .......... 51% 57.8M 5s
    ##  69950K .......... .......... .......... .......... .......... 51% 55.4M 5s
    ##  70000K .......... .......... .......... .......... .......... 51% 59.4M 5s
    ##  70050K .......... .......... .......... .......... .......... 52% 68.1M 5s
    ##  70100K .......... .......... .......... .......... .......... 52% 61.9M 5s
    ##  70150K .......... .......... .......... .......... .......... 52% 60.0M 5s
    ##  70200K .......... .......... .......... .......... .......... 52% 69.8M 5s
    ##  70250K .......... .......... .......... .......... .......... 52% 69.8M 5s
    ##  70300K .......... .......... .......... .......... .......... 52% 70.1M 5s
    ##  70350K .......... .......... .......... .......... .......... 52% 59.1M 5s
    ##  70400K .......... .......... .......... .......... .......... 52% 59.8M 5s
    ##  70450K .......... .......... .......... .......... .......... 52% 21.4M 5s
    ##  70500K .......... .......... .......... .......... .......... 52% 63.3M 5s
    ##  70550K .......... .......... .......... .......... .......... 52% 69.6M 5s
    ##  70600K .......... .......... .......... .......... .......... 52% 71.7M 5s
    ##  70650K .......... .......... .......... .......... .......... 52% 65.0M 5s
    ##  70700K .......... .......... .......... .......... .......... 52% 18.3M 5s
    ##  70750K .......... .......... .......... .......... .......... 52% 76.4M 5s
    ##  70800K .......... .......... .......... .......... .......... 52% 70.1M 5s
    ##  70850K .......... .......... .......... .......... .......... 52% 77.4M 5s
    ##  70900K .......... .......... .......... .......... .......... 52% 71.7M 5s
    ##  70950K .......... .......... .......... .......... .......... 52% 26.0M 5s
    ##  71000K .......... .......... .......... .......... .......... 52% 67.0M 5s
    ##  71050K .......... .......... .......... .......... .......... 52% 61.4M 5s
    ##  71100K .......... .......... .......... .......... .......... 52% 71.7M 5s
    ##  71150K .......... .......... .......... .......... .......... 52% 70.1M 5s
    ##  71200K .......... .......... .......... .......... .......... 52% 42.8M 5s
    ##  71250K .......... .......... .......... .......... .......... 52% 25.8M 5s
    ##  71300K .......... .......... .......... .......... .......... 52% 64.4M 5s
    ##  71350K .......... .......... .......... .......... .......... 52% 53.6M 5s
    ##  71400K .......... .......... .......... .......... .......... 53% 64.5M 5s
    ##  71450K .......... .......... .......... .......... .......... 53% 72.1M 5s
    ##  71500K .......... .......... .......... .......... .......... 53% 30.6M 5s
    ##  71550K .......... .......... .......... .......... .......... 53% 56.7M 5s
    ##  71600K .......... .......... .......... .......... .......... 53% 61.4M 5s
    ##  71650K .......... .......... .......... .......... .......... 53% 73.1M 5s
    ##  71700K .......... .......... .......... .......... .......... 53% 69.8M 5s
    ##  71750K .......... .......... .......... .......... .......... 53% 53.4M 5s
    ##  71800K .......... .......... .......... .......... .......... 53% 46.5M 5s
    ##  71850K .......... .......... .......... .......... .......... 53% 78.5M 5s
    ##  71900K .......... .......... .......... .......... .......... 53% 70.2M 5s
    ##  71950K .......... .......... .......... .......... .......... 53% 71.9M 5s
    ##  72000K .......... .......... .......... .......... .......... 53% 67.1M 5s
    ##  72050K .......... .......... .......... .......... .......... 53% 62.0M 5s
    ##  72100K .......... .......... .......... .......... .......... 53% 43.6M 5s
    ##  72150K .......... .......... .......... .......... .......... 53% 66.6M 5s
    ##  72200K .......... .......... .......... .......... .......... 53% 65.6M 5s
    ##  72250K .......... .......... .......... .......... .......... 53% 69.8M 5s
    ##  72300K .......... .......... .......... .......... .......... 53% 58.3M 5s
    ##  72350K .......... .......... .......... .......... .......... 53% 44.1M 5s
    ##  72400K .......... .......... .......... .......... .......... 53% 53.2M 5s
    ##  72450K .......... .......... .......... .......... .......... 53% 63.7M 5s
    ##  72500K .......... .......... .......... .......... .......... 53% 61.7M 5s
    ##  72550K .......... .......... .......... .......... .......... 53% 69.3M 5s
    ##  72600K .......... .......... .......... .......... .......... 53% 18.2M 5s
    ##  72650K .......... .......... .......... .......... .......... 53% 41.8M 5s
    ##  72700K .......... .......... .......... .......... .......... 53% 43.6M 5s
    ##  72750K .......... .......... .......... .......... .......... 54% 48.1M 5s
    ##  72800K .......... .......... .......... .......... .......... 54% 51.9M 5s
    ##  72850K .......... .......... .......... .......... .......... 54% 51.3M 5s
    ##  72900K .......... .......... .......... .......... .......... 54% 48.4M 5s
    ##  72950K .......... .......... .......... .......... .......... 54% 33.8M 5s
    ##  73000K .......... .......... .......... .......... .......... 54% 47.9M 5s
    ##  73050K .......... .......... .......... .......... .......... 54% 51.1M 5s
    ##  73100K .......... .......... .......... .......... .......... 54% 52.8M 5s
    ##  73150K .......... .......... .......... .......... .......... 54% 35.4M 5s
    ##  73200K .......... .......... .......... .......... .......... 54% 38.7M 5s
    ##  73250K .......... .......... .......... .......... .......... 54% 49.5M 5s
    ##  73300K .......... .......... .......... .......... .......... 54% 60.9M 5s
    ##  73350K .......... .......... .......... .......... .......... 54% 53.2M 5s
    ##  73400K .......... .......... .......... .......... .......... 54% 54.3M 5s
    ##  73450K .......... .......... .......... .......... .......... 54% 62.3M 5s
    ##  73500K .......... .......... .......... .......... .......... 54% 55.7M 5s
    ##  73550K .......... .......... .......... .......... .......... 54% 53.7M 5s
    ##  73600K .......... .......... .......... .......... .......... 54% 54.8M 5s
    ##  73650K .......... .......... .......... .......... .......... 54% 56.2M 5s
    ##  73700K .......... .......... .......... .......... .......... 54% 55.5M 5s
    ##  73750K .......... .......... .......... .......... .......... 54% 59.9M 5s
    ##  73800K .......... .......... .......... .......... .......... 54% 54.7M 5s
    ##  73850K .......... .......... .......... .......... .......... 54% 53.2M 5s
    ##  73900K .......... .......... .......... .......... .......... 54% 57.3M 5s
    ##  73950K .......... .......... .......... .......... .......... 54% 59.2M 5s
    ##  74000K .......... .......... .......... .......... .......... 54% 60.0M 5s
    ##  74050K .......... .......... .......... .......... .......... 54% 52.4M 5s
    ##  74100K .......... .......... .......... .......... .......... 55% 55.5M 5s
    ##  74150K .......... .......... .......... .......... .......... 55% 62.6M 5s
    ##  74200K .......... .......... .......... .......... .......... 55% 56.6M 5s
    ##  74250K .......... .......... .......... .......... .......... 55% 70.7M 4s
    ##  74300K .......... .......... .......... .......... .......... 55% 58.4M 4s
    ##  74350K .......... .......... .......... .......... .......... 55% 70.4M 4s
    ##  74400K .......... .......... .......... .......... .......... 55% 56.2M 4s
    ##  74450K .......... .......... .......... .......... .......... 55% 53.5M 4s
    ##  74500K .......... .......... .......... .......... .......... 55% 49.6M 4s
    ##  74550K .......... .......... .......... .......... .......... 55% 71.5M 4s
    ##  74600K .......... .......... .......... .......... .......... 55% 55.4M 4s
    ##  74650K .......... .......... .......... .......... .......... 55% 75.8M 4s
    ##  74700K .......... .......... .......... .......... .......... 55% 64.4M 4s
    ##  74750K .......... .......... .......... .......... .......... 55% 57.9M 4s
    ##  74800K .......... .......... .......... .......... .......... 55% 55.3M 4s
    ##  74850K .......... .......... .......... .......... .......... 55% 67.3M 4s
    ##  74900K .......... .......... .......... .......... .......... 55% 63.6M 4s
    ##  74950K .......... .......... .......... .......... .......... 55% 66.0M 4s
    ##  75000K .......... .......... .......... .......... .......... 55% 61.3M 4s
    ##  75050K .......... .......... .......... .......... .......... 55% 68.2M 4s
    ##  75100K .......... .......... .......... .......... .......... 55% 60.6M 4s
    ##  75150K .......... .......... .......... .......... .......... 55% 73.1M 4s
    ##  75200K .......... .......... .......... .......... .......... 55% 67.9M 4s
    ##  75250K .......... .......... .......... .......... .......... 55% 80.8M 4s
    ##  75300K .......... .......... .......... .......... .......... 55% 44.3M 4s
    ##  75350K .......... .......... .......... .......... .......... 55% 25.6M 4s
    ##  75400K .......... .......... .......... .......... .......... 55% 56.4M 4s
    ##  75450K .......... .......... .......... .......... .......... 56% 76.1M 4s
    ##  75500K .......... .......... .......... .......... .......... 56% 73.3M 4s
    ##  75550K .......... .......... .......... .......... .......... 56% 76.2M 4s
    ##  75600K .......... .......... .......... .......... .......... 56% 57.8M 4s
    ##  75650K .......... .......... .......... .......... .......... 56% 67.0M 4s
    ##  75700K .......... .......... .......... .......... .......... 56% 56.2M 4s
    ##  75750K .......... .......... .......... .......... .......... 56% 77.3M 4s
    ##  75800K .......... .......... .......... .......... .......... 56% 35.2M 4s
    ##  75850K .......... .......... .......... .......... .......... 56% 64.9M 4s
    ##  75900K .......... .......... .......... .......... .......... 56% 23.5M 4s
    ##  75950K .......... .......... .......... .......... .......... 56% 61.7M 4s
    ##  76000K .......... .......... .......... .......... .......... 56% 60.2M 4s
    ##  76050K .......... .......... .......... .......... .......... 56% 70.3M 4s
    ##  76100K .......... .......... .......... .......... .......... 56% 70.2M 4s
    ##  76150K .......... .......... .......... .......... .......... 56% 39.3M 4s
    ##  76200K .......... .......... .......... .......... .......... 56% 48.9M 4s
    ##  76250K .......... .......... .......... .......... .......... 56% 49.9M 4s
    ##  76300K .......... .......... .......... .......... .......... 56% 62.7M 4s
    ##  76350K .......... .......... .......... .......... .......... 56% 81.3M 4s
    ##  76400K .......... .......... .......... .......... .......... 56% 42.0M 4s
    ##  76450K .......... .......... .......... .......... .......... 56% 18.5M 4s
    ##  76500K .......... .......... .......... .......... .......... 56% 72.4M 4s
    ##  76550K .......... .......... .......... .......... .......... 56% 68.0M 4s
    ##  76600K .......... .......... .......... .......... .......... 56% 73.4M 4s
    ##  76650K .......... .......... .......... .......... .......... 56% 84.6M 4s
    ##  76700K .......... .......... .......... .......... .......... 56% 31.5M 4s
    ##  76750K .......... .......... .......... .......... .......... 56% 63.4M 4s
    ##  76800K .......... .......... .......... .......... .......... 57% 1.09M 4s
    ##  76850K .......... .......... .......... .......... .......... 57% 41.9M 4s
    ##  76900K .......... .......... .......... .......... .......... 57% 44.2M 4s
    ##  76950K .......... .......... .......... .......... .......... 57% 42.9M 4s
    ##  77000K .......... .......... .......... .......... .......... 57% 41.6M 4s
    ##  77050K .......... .......... .......... .......... .......... 57% 47.3M 4s
    ##  77100K .......... .......... .......... .......... .......... 57% 41.8M 4s
    ##  77150K .......... .......... .......... .......... .......... 57% 40.4M 4s
    ##  77200K .......... .......... .......... .......... .......... 57% 46.8M 4s
    ##  77250K .......... .......... .......... .......... .......... 57% 41.8M 4s
    ##  77300K .......... .......... .......... .......... .......... 57% 45.1M 4s
    ##  77350K .......... .......... .......... .......... .......... 57% 52.5M 4s
    ##  77400K .......... .......... .......... .......... .......... 57% 45.2M 4s
    ##  77450K .......... .......... .......... .......... .......... 57% 47.5M 4s
    ##  77500K .......... .......... .......... .......... .......... 57% 43.1M 4s
    ##  77550K .......... .......... .......... .......... .......... 57% 49.1M 4s
    ##  77600K .......... .......... .......... .......... .......... 57% 43.8M 4s
    ##  77650K .......... .......... .......... .......... .......... 57% 42.5M 4s
    ##  77700K .......... .......... .......... .......... .......... 57% 36.7M 4s
    ##  77750K .......... .......... .......... .......... .......... 57% 45.7M 4s
    ##  77800K .......... .......... .......... .......... .......... 57% 47.7M 4s
    ##  77850K .......... .......... .......... .......... .......... 57% 54.1M 4s
    ##  77900K .......... .......... .......... .......... .......... 57% 33.6M 4s
    ##  77950K .......... .......... .......... .......... .......... 57% 47.3M 4s
    ##  78000K .......... .......... .......... .......... .......... 57% 43.3M 4s
    ##  78050K .......... .......... .......... .......... .......... 57% 39.2M 4s
    ##  78100K .......... .......... .......... .......... .......... 58% 40.2M 4s
    ##  78150K .......... .......... .......... .......... .......... 58% 42.5M 4s
    ##  78200K .......... .......... .......... .......... .......... 58% 40.3M 4s
    ##  78250K .......... .......... .......... .......... .......... 58% 46.5M 4s
    ##  78300K .......... .......... .......... .......... .......... 58% 45.2M 4s
    ##  78350K .......... .......... .......... .......... .......... 58% 40.6M 4s
    ##  78400K .......... .......... .......... .......... .......... 58% 45.6M 4s
    ##  78450K .......... .......... .......... .......... .......... 58% 46.2M 4s
    ##  78500K .......... .......... .......... .......... .......... 58% 49.3M 4s
    ##  78550K .......... .......... .......... .......... .......... 58% 63.8M 4s
    ##  78600K .......... .......... .......... .......... .......... 58% 52.1M 4s
    ##  78650K .......... .......... .......... .......... .......... 58% 52.7M 4s
    ##  78700K .......... .......... .......... .......... .......... 58% 61.3M 4s
    ##  78750K .......... .......... .......... .......... .......... 58% 61.2M 4s
    ##  78800K .......... .......... .......... .......... .......... 58% 47.9M 4s
    ##  78850K .......... .......... .......... .......... .......... 58% 57.4M 4s
    ##  78900K .......... .......... .......... .......... .......... 58% 59.0M 4s
    ##  78950K .......... .......... .......... .......... .......... 58% 55.6M 4s
    ##  79000K .......... .......... .......... .......... .......... 58% 61.3M 4s
    ##  79050K .......... .......... .......... .......... .......... 58% 59.6M 4s
    ##  79100K .......... .......... .......... .......... .......... 58% 49.0M 4s
    ##  79150K .......... .......... .......... .......... .......... 58% 56.4M 4s
    ##  79200K .......... .......... .......... .......... .......... 58% 48.8M 4s
    ##  79250K .......... .......... .......... .......... .......... 58% 47.6M 4s
    ##  79300K .......... .......... .......... .......... .......... 58% 54.7M 4s
    ##  79350K .......... .......... .......... .......... .......... 58% 67.2M 4s
    ##  79400K .......... .......... .......... .......... .......... 58% 59.9M 4s
    ##  79450K .......... .......... .......... .......... .......... 59% 70.5M 4s
    ##  79500K .......... .......... .......... .......... .......... 59% 61.0M 4s
    ##  79550K .......... .......... .......... .......... .......... 59% 16.9M 4s
    ##  79600K .......... .......... .......... .......... .......... 59% 51.7M 4s
    ##  79650K .......... .......... .......... .......... .......... 59% 54.1M 4s
    ##  79700K .......... .......... .......... .......... .......... 59% 57.3M 4s
    ##  79750K .......... .......... .......... .......... .......... 59% 72.6M 4s
    ##  79800K .......... .......... .......... .......... .......... 59% 24.4M 4s
    ##  79850K .......... .......... .......... .......... .......... 59% 48.5M 4s
    ##  79900K .......... .......... .......... .......... .......... 59% 54.0M 4s
    ##  79950K .......... .......... .......... .......... .......... 59% 69.0M 4s
    ##  80000K .......... .......... .......... .......... .......... 59% 64.6M 4s
    ##  80050K .......... .......... .......... .......... .......... 59% 48.6M 4s
    ##  80100K .......... .......... .......... .......... .......... 59% 54.5M 4s
    ##  80150K .......... .......... .......... .......... .......... 59% 58.5M 4s
    ##  80200K .......... .......... .......... .......... .......... 59% 52.1M 4s
    ##  80250K .......... .......... .......... .......... .......... 59% 72.2M 4s
    ##  80300K .......... .......... .......... .......... .......... 59% 62.2M 4s
    ##  80350K .......... .......... .......... .......... .......... 59% 72.5M 4s
    ##  80400K .......... .......... .......... .......... .......... 59% 60.8M 4s
    ##  80450K .......... .......... .......... .......... .......... 59% 62.5M 4s
    ##  80500K .......... .......... .......... .......... .......... 59% 69.8M 4s
    ##  80550K .......... .......... .......... .......... .......... 59% 57.2M 4s
    ##  80600K .......... .......... .......... .......... .......... 59% 46.1M 4s
    ##  80650K .......... .......... .......... .......... .......... 59% 54.2M 4s
    ##  80700K .......... .......... .......... .......... .......... 59% 49.1M 4s
    ##  80750K .......... .......... .......... .......... .......... 59% 62.2M 4s
    ##  80800K .......... .......... .......... .......... .......... 60% 54.9M 4s
    ##  80850K .......... .......... .......... .......... .......... 60% 23.3M 4s
    ##  80900K .......... .......... .......... .......... .......... 60% 33.3M 4s
    ##  80950K .......... .......... .......... .......... .......... 60% 43.4M 4s
    ##  81000K .......... .......... .......... .......... .......... 60% 43.0M 4s
    ##  81050K .......... .......... .......... .......... .......... 60% 40.5M 4s
    ##  81100K .......... .......... .......... .......... .......... 60% 36.5M 4s
    ##  81150K .......... .......... .......... .......... .......... 60% 51.7M 4s
    ##  81200K .......... .......... .......... .......... .......... 60% 45.8M 4s
    ##  81250K .......... .......... .......... .......... .......... 60% 45.4M 4s
    ##  81300K .......... .......... .......... .......... .......... 60% 41.5M 4s
    ##  81350K .......... .......... .......... .......... .......... 60% 46.2M 4s
    ##  81400K .......... .......... .......... .......... .......... 60% 38.6M 4s
    ##  81450K .......... .......... .......... .......... .......... 60% 55.1M 4s
    ##  81500K .......... .......... .......... .......... .......... 60% 42.9M 4s
    ##  81550K .......... .......... .......... .......... .......... 60% 42.5M 4s
    ##  81600K .......... .......... .......... .......... .......... 60% 45.5M 4s
    ##  81650K .......... .......... .......... .......... .......... 60% 58.4M 4s
    ##  81700K .......... .......... .......... .......... .......... 60% 52.7M 4s
    ##  81750K .......... .......... .......... .......... .......... 60% 64.6M 4s
    ##  81800K .......... .......... .......... .......... .......... 60% 54.6M 4s
    ##  81850K .......... .......... .......... .......... .......... 60% 45.3M 4s
    ##  81900K .......... .......... .......... .......... .......... 60% 42.9M 4s
    ##  81950K .......... .......... .......... .......... .......... 60% 60.3M 4s
    ##  82000K .......... .......... .......... .......... .......... 60% 48.6M 4s
    ##  82050K .......... .......... .......... .......... .......... 60% 55.9M 4s
    ##  82100K .......... .......... .......... .......... .......... 60% 54.4M 4s
    ##  82150K .......... .......... .......... .......... .......... 61% 56.3M 4s
    ##  82200K .......... .......... .......... .......... .......... 61% 56.8M 4s
    ##  82250K .......... .......... .......... .......... .......... 61% 51.9M 4s
    ##  82300K .......... .......... .......... .......... .......... 61% 51.9M 4s
    ##  82350K .......... .......... .......... .......... .......... 61% 60.4M 4s
    ##  82400K .......... .......... .......... .......... .......... 61% 57.6M 4s
    ##  82450K .......... .......... .......... .......... .......... 61% 63.5M 4s
    ##  82500K .......... .......... .......... .......... .......... 61% 50.5M 4s
    ##  82550K .......... .......... .......... .......... .......... 61% 60.5M 4s
    ##  82600K .......... .......... .......... .......... .......... 61% 51.5M 4s
    ##  82650K .......... .......... .......... .......... .......... 61% 65.4M 4s
    ##  82700K .......... .......... .......... .......... .......... 61% 51.7M 4s
    ##  82750K .......... .......... .......... .......... .......... 61% 9.40M 4s
    ##  82800K .......... .......... .......... .......... .......... 61% 62.9M 4s
    ##  82850K .......... .......... .......... .......... .......... 61% 63.8M 4s
    ##  82900K .......... .......... .......... .......... .......... 61% 62.4M 4s
    ##  82950K .......... .......... .......... .......... .......... 61% 68.1M 4s
    ##  83000K .......... .......... .......... .......... .......... 61% 60.4M 4s
    ##  83050K .......... .......... .......... .......... .......... 61% 62.1M 4s
    ##  83100K .......... .......... .......... .......... .......... 61% 61.9M 4s
    ##  83150K .......... .......... .......... .......... .......... 61% 71.3M 4s
    ##  83200K .......... .......... .......... .......... .......... 61% 18.2M 4s
    ##  83250K .......... .......... .......... .......... .......... 61% 73.3M 4s
    ##  83300K .......... .......... .......... .......... .......... 61% 58.6M 4s
    ##  83350K .......... .......... .......... .......... .......... 61% 49.3M 4s
    ##  83400K .......... .......... .......... .......... .......... 61% 51.3M 4s
    ##  83450K .......... .......... .......... .......... .......... 61% 69.1M 4s
    ##  83500K .......... .......... .......... .......... .......... 62% 55.6M 4s
    ##  83550K .......... .......... .......... .......... .......... 62% 72.0M 4s
    ##  83600K .......... .......... .......... .......... .......... 62% 56.6M 4s
    ##  83650K .......... .......... .......... .......... .......... 62% 60.4M 4s
    ##  83700K .......... .......... .......... .......... .......... 62% 61.0M 4s
    ##  83750K .......... .......... .......... .......... .......... 62% 72.6M 3s
    ##  83800K .......... .......... .......... .......... .......... 62% 64.3M 3s
    ##  83850K .......... .......... .......... .......... .......... 62% 59.6M 3s
    ##  83900K .......... .......... .......... .......... .......... 62% 55.4M 3s
    ##  83950K .......... .......... .......... .......... .......... 62% 11.7M 3s
    ##  84000K .......... .......... .......... .......... .......... 62% 64.1M 3s
    ##  84050K .......... .......... .......... .......... .......... 62% 22.1M 3s
    ##  84100K .......... .......... .......... .......... .......... 62% 53.0M 3s
    ##  84150K .......... .......... .......... .......... .......... 62% 73.2M 3s
    ##  84200K .......... .......... .......... .......... .......... 62% 13.5M 3s
    ##  84250K .......... .......... .......... .......... .......... 62% 69.3M 3s
    ##  84300K .......... .......... .......... .......... .......... 62% 62.2M 3s
    ##  84350K .......... .......... .......... .......... .......... 62% 21.4M 3s
    ##  84400K .......... .......... .......... .......... .......... 62% 63.3M 3s
    ##  84450K .......... .......... .......... .......... .......... 62% 69.7M 3s
    ##  84500K .......... .......... .......... .......... .......... 62% 17.2M 3s
    ##  84550K .......... .......... .......... .......... .......... 62% 43.5M 3s
    ##  84600K .......... .......... .......... .......... .......... 62% 72.5M 3s
    ##  84650K .......... .......... .......... .......... .......... 62% 13.6M 3s
    ##  84700K .......... .......... .......... .......... .......... 62% 69.9M 3s
    ##  84750K .......... .......... .......... .......... .......... 62% 86.9M 3s
    ##  84800K .......... .......... .......... .......... .......... 62% 19.6M 3s
    ##  84850K .......... .......... .......... .......... .......... 63% 65.7M 3s
    ##  84900K .......... .......... .......... .......... .......... 63% 20.0M 3s
    ##  84950K .......... .......... .......... .......... .......... 63% 49.9M 3s
    ##  85000K .......... .......... .......... .......... .......... 63% 27.2M 3s
    ##  85050K .......... .......... .......... .......... .......... 63% 61.1M 3s
    ##  85100K .......... .......... .......... .......... .......... 63% 15.2M 3s
    ##  85150K .......... .......... .......... .......... .......... 63% 18.4M 3s
    ##  85200K .......... .......... .......... .......... .......... 63% 29.1M 3s
    ##  85250K .......... .......... .......... .......... .......... 63% 72.2M 3s
    ##  85300K .......... .......... .......... .......... .......... 63% 28.7M 3s
    ##  85350K .......... .......... .......... .......... .......... 63% 27.8M 3s
    ##  85400K .......... .......... .......... .......... .......... 63% 21.9M 3s
    ##  85450K .......... .......... .......... .......... .......... 63% 71.6M 3s
    ##  85500K .......... .......... .......... .......... .......... 63% 68.1M 3s
    ##  85550K .......... .......... .......... .......... .......... 63% 13.4M 3s
    ##  85600K .......... .......... .......... .......... .......... 63% 80.6M 3s
    ##  85650K .......... .......... .......... .......... .......... 63% 18.0M 3s
    ##  85700K .......... .......... .......... .......... .......... 63% 57.7M 3s
    ##  85750K .......... .......... .......... .......... .......... 63% 91.4M 3s
    ##  85800K .......... .......... .......... .......... .......... 63% 21.0M 3s
    ##  85850K .......... .......... .......... .......... .......... 63% 51.4M 3s
    ##  85900K .......... .......... .......... .......... .......... 63% 15.6M 3s
    ##  85950K .......... .......... .......... .......... .......... 63% 90.6M 3s
    ##  86000K .......... .......... .......... .......... .......... 63% 69.0M 3s
    ##  86050K .......... .......... .......... .......... .......... 63% 9.19M 3s
    ##  86100K .......... .......... .......... .......... .......... 63% 16.2M 3s
    ##  86150K .......... .......... .......... .......... .......... 63% 49.6M 3s
    ##  86200K .......... .......... .......... .......... .......... 64% 20.0M 3s
    ##  86250K .......... .......... .......... .......... .......... 64% 85.6M 3s
    ##  86300K .......... .......... .......... .......... .......... 64% 51.7M 3s
    ##  86350K .......... .......... .......... .......... .......... 64% 28.9M 3s
    ##  86400K .......... .......... .......... .......... .......... 64% 26.3M 3s
    ##  86450K .......... .......... .......... .......... .......... 64% 79.0M 3s
    ##  86500K .......... .......... .......... .......... .......... 64% 21.6M 3s
    ##  86550K .......... .......... .......... .......... .......... 64% 70.9M 3s
    ##  86600K .......... .......... .......... .......... .......... 64% 16.1M 3s
    ##  86650K .......... .......... .......... .......... .......... 64% 48.1M 3s
    ##  86700K .......... .......... .......... .......... .......... 64% 82.1M 3s
    ##  86750K .......... .......... .......... .......... .......... 64% 20.6M 3s
    ##  86800K .......... .......... .......... .......... .......... 64% 70.9M 3s
    ##  86850K .......... .......... .......... .......... .......... 64% 59.5M 3s
    ##  86900K .......... .......... .......... .......... .......... 64% 73.2M 3s
    ##  86950K .......... .......... .......... .......... .......... 64% 19.0M 3s
    ##  87000K .......... .......... .......... .......... .......... 64% 62.6M 3s
    ##  87050K .......... .......... .......... .......... .......... 64% 30.4M 3s
    ##  87100K .......... .......... .......... .......... .......... 64% 26.1M 3s
    ##  87150K .......... .......... .......... .......... .......... 64% 70.9M 3s
    ##  87200K .......... .......... .......... .......... .......... 64% 35.3M 3s
    ##  87250K .......... .......... .......... .......... .......... 64% 15.3M 3s
    ##  87300K .......... .......... .......... .......... .......... 64% 70.4M 3s
    ##  87350K .......... .......... .......... .......... .......... 64% 88.4M 3s
    ##  87400K .......... .......... .......... .......... .......... 64% 14.1M 3s
    ##  87450K .......... .......... .......... .......... .......... 64% 61.2M 3s
    ##  87500K .......... .......... .......... .......... .......... 64% 18.8M 3s
    ##  87550K .......... .......... .......... .......... .......... 65% 58.8M 3s
    ##  87600K .......... .......... .......... .......... .......... 65% 73.6M 3s
    ##  87650K .......... .......... .......... .......... .......... 65% 15.8M 3s
    ##  87700K .......... .......... .......... .......... .......... 65% 59.4M 3s
    ##  87750K .......... .......... .......... .......... .......... 65% 76.4M 3s
    ##  87800K .......... .......... .......... .......... .......... 65% 24.7M 3s
    ##  87850K .......... .......... .......... .......... .......... 65% 25.9M 3s
    ##  87900K .......... .......... .......... .......... .......... 65% 74.3M 3s
    ##  87950K .......... .......... .......... .......... .......... 65% 40.7M 3s
    ##  88000K .......... .......... .......... .......... .......... 65% 31.4M 3s
    ##  88050K .......... .......... .......... .......... .......... 65% 88.2M 3s
    ##  88100K .......... .......... .......... .......... .......... 65% 23.1M 3s
    ##  88150K .......... .......... .......... .......... .......... 65% 34.5M 3s
    ##  88200K .......... .......... .......... .......... .......... 65% 20.7M 3s
    ##  88250K .......... .......... .......... .......... .......... 65% 70.2M 3s
    ##  88300K .......... .......... .......... .......... .......... 65% 19.1M 3s
    ##  88350K .......... .......... .......... .......... .......... 65% 85.1M 3s
    ##  88400K .......... .......... .......... .......... .......... 65% 34.9M 3s
    ##  88450K .......... .......... .......... .......... .......... 65% 23.0M 3s
    ##  88500K .......... .......... .......... .......... .......... 65% 27.7M 3s
    ##  88550K .......... .......... .......... .......... .......... 65% 81.4M 3s
    ##  88600K .......... .......... .......... .......... .......... 65% 39.8M 3s
    ##  88650K .......... .......... .......... .......... .......... 65% 21.6M 3s
    ##  88700K .......... .......... .......... .......... .......... 65% 54.5M 3s
    ##  88750K .......... .......... .......... .......... .......... 65% 19.7M 3s
    ##  88800K .......... .......... .......... .......... .......... 65% 74.9M 3s
    ##  88850K .......... .......... .......... .......... .......... 65% 54.4M 3s
    ##  88900K .......... .......... .......... .......... .......... 66% 16.4M 3s
    ##  88950K .......... .......... .......... .......... .......... 66% 43.2M 3s
    ##  89000K .......... .......... .......... .......... .......... 66% 52.2M 3s
    ##  89050K .......... .......... .......... .......... .......... 66% 46.3M 3s
    ##  89100K .......... .......... .......... .......... .......... 66% 29.4M 3s
    ##  89150K .......... .......... .......... .......... .......... 66% 46.2M 3s
    ##  89200K .......... .......... .......... .......... .......... 66% 17.6M 3s
    ##  89250K .......... .......... .......... .......... .......... 66% 44.5M 3s
    ##  89300K .......... .......... .......... .......... .......... 66% 53.7M 3s
    ##  89350K .......... .......... .......... .......... .......... 66% 17.9M 3s
    ##  89400K .......... .......... .......... .......... .......... 66% 46.4M 3s
    ##  89450K .......... .......... .......... .......... .......... 66% 16.1M 3s
    ##  89500K .......... .......... .......... .......... .......... 66% 42.4M 3s
    ##  89550K .......... .......... .......... .......... .......... 66% 47.3M 3s
    ##  89600K .......... .......... .......... .......... .......... 66% 12.5M 3s
    ##  89650K .......... .......... .......... .......... .......... 66% 45.9M 3s
    ##  89700K .......... .......... .......... .......... .......... 66% 50.6M 3s
    ##  89750K .......... .......... .......... .......... .......... 66% 27.5M 3s
    ##  89800K .......... .......... .......... .......... .......... 66% 38.3M 3s
    ##  89850K .......... .......... .......... .......... .......... 66% 51.8M 3s
    ##  89900K .......... .......... .......... .......... .......... 66% 24.6M 3s
    ##  89950K .......... .......... .......... .......... .......... 66% 37.5M 3s
    ##  90000K .......... .......... .......... .......... .......... 66% 53.4M 3s
    ##  90050K .......... .......... .......... .......... .......... 66% 13.4M 3s
    ##  90100K .......... .......... .......... .......... .......... 66% 55.2M 3s
    ##  90150K .......... .......... .......... .......... .......... 66% 25.2M 3s
    ##  90200K .......... .......... .......... .......... .......... 66% 44.7M 3s
    ##  90250K .......... .......... .......... .......... .......... 67% 51.9M 3s
    ##  90300K .......... .......... .......... .......... .......... 67% 17.0M 3s
    ##  90350K .......... .......... .......... .......... .......... 67% 58.2M 3s
    ##  90400K .......... .......... .......... .......... .......... 67% 59.0M 3s
    ##  90450K .......... .......... .......... .......... .......... 67% 58.6M 3s
    ##  90500K .......... .......... .......... .......... .......... 67% 24.6M 3s
    ##  90550K .......... .......... .......... .......... .......... 67% 28.0M 3s
    ##  90600K .......... .......... .......... .......... .......... 67% 33.7M 3s
    ##  90650K .......... .......... .......... .......... .......... 67% 60.3M 3s
    ##  90700K .......... .......... .......... .......... .......... 67% 25.3M 3s
    ##  90750K .......... .......... .......... .......... .......... 67% 41.6M 3s
    ##  90800K .......... .......... .......... .......... .......... 67% 23.9M 3s
    ##  90850K .......... .......... .......... .......... .......... 67% 33.7M 3s
    ##  90900K .......... .......... .......... .......... .......... 67% 57.5M 3s
    ##  90950K .......... .......... .......... .......... .......... 67% 24.5M 3s
    ##  91000K .......... .......... .......... .......... .......... 67% 43.5M 3s
    ##  91050K .......... .......... .......... .......... .......... 67% 22.2M 3s
    ##  91100K .......... .......... .......... .......... .......... 67% 28.0M 3s
    ##  91150K .......... .......... .......... .......... .......... 67% 61.1M 3s
    ##  91200K .......... .......... .......... .......... .......... 67% 11.3M 3s
    ##  91250K .......... .......... .......... .......... .......... 67% 26.2M 3s
    ##  91300K .......... .......... .......... .......... .......... 67% 26.6M 3s
    ##  91350K .......... .......... .......... .......... .......... 67% 32.0M 3s
    ##  91400K .......... .......... .......... .......... .......... 67% 36.8M 3s
    ##  91450K .......... .......... .......... .......... .......... 67% 47.1M 3s
    ##  91500K .......... .......... .......... .......... .......... 67% 23.1M 3s
    ##  91550K .......... .......... .......... .......... .......... 67% 23.8M 3s
    ##  91600K .......... .......... .......... .......... .......... 68% 18.2M 3s
    ##  91650K .......... .......... .......... .......... .......... 68% 43.4M 3s
    ##  91700K .......... .......... .......... .......... .......... 68% 41.9M 3s
    ##  91750K .......... .......... .......... .......... .......... 68% 24.0M 3s
    ##  91800K .......... .......... .......... .......... .......... 68% 44.2M 3s
    ##  91850K .......... .......... .......... .......... .......... 68% 63.6M 3s
    ##  91900K .......... .......... .......... .......... .......... 68% 26.2M 3s
    ##  91950K .......... .......... .......... .......... .......... 68% 31.5M 3s
    ##  92000K .......... .......... .......... .......... .......... 68% 18.3M 3s
    ##  92050K .......... .......... .......... .......... .......... 68% 51.2M 3s
    ##  92100K .......... .......... .......... .......... .......... 68% 59.3M 3s
    ##  92150K .......... .......... .......... .......... .......... 68% 53.5M 3s
    ##  92200K .......... .......... .......... .......... .......... 68% 19.2M 3s
    ##  92250K .......... .......... .......... .......... .......... 68% 26.1M 3s
    ##  92300K .......... .......... .......... .......... .......... 68% 20.2M 3s
    ##  92350K .......... .......... .......... .......... .......... 68% 60.3M 3s
    ##  92400K .......... .......... .......... .......... .......... 68% 50.1M 3s
    ##  92450K .......... .......... .......... .......... .......... 68% 19.4M 3s
    ##  92500K .......... .......... .......... .......... .......... 68% 47.1M 3s
    ##  92550K .......... .......... .......... .......... .......... 68% 43.6M 3s
    ##  92600K .......... .......... .......... .......... .......... 68% 49.2M 3s
    ##  92650K .......... .......... .......... .......... .......... 68% 21.5M 3s
    ##  92700K .......... .......... .......... .......... .......... 68% 51.2M 3s
    ##  92750K .......... .......... .......... .......... .......... 68% 29.9M 3s
    ##  92800K .......... .......... .......... .......... .......... 68% 26.6M 3s
    ##  92850K .......... .......... .......... .......... .......... 68% 63.6M 3s
    ##  92900K .......... .......... .......... .......... .......... 68% 13.8M 3s
    ##  92950K .......... .......... .......... .......... .......... 69% 60.6M 3s
    ##  93000K .......... .......... .......... .......... .......... 69% 9.40M 3s
    ##  93050K .......... .......... .......... .......... .......... 69% 58.8M 3s
    ##  93100K .......... .......... .......... .......... .......... 69% 67.8M 3s
    ##  93150K .......... .......... .......... .......... .......... 69% 17.0M 3s
    ##  93200K .......... .......... .......... .......... .......... 69% 66.6M 3s
    ##  93250K .......... .......... .......... .......... .......... 69% 45.6M 3s
    ##  93300K .......... .......... .......... .......... .......... 69% 52.4M 3s
    ##  93350K .......... .......... .......... .......... .......... 69% 26.0M 3s
    ##  93400K .......... .......... .......... .......... .......... 69% 27.3M 3s
    ##  93450K .......... .......... .......... .......... .......... 69% 14.8M 3s
    ##  93500K .......... .......... .......... .......... .......... 69% 61.1M 3s
    ##  93550K .......... .......... .......... .......... .......... 69% 73.6M 3s
    ##  93600K .......... .......... .......... .......... .......... 69% 16.1M 3s
    ##  93650K .......... .......... .......... .......... .......... 69% 76.6M 3s
    ##  93700K .......... .......... .......... .......... .......... 69% 14.8M 3s
    ##  93750K .......... .......... .......... .......... .......... 69% 52.6M 3s
    ##  93800K .......... .......... .......... .......... .......... 69% 76.0M 3s
    ##  93850K .......... .......... .......... .......... .......... 69% 11.2M 3s
    ##  93900K .......... .......... .......... .......... .......... 69% 70.4M 3s
    ##  93950K .......... .......... .......... .......... .......... 69% 75.2M 3s
    ##  94000K .......... .......... .......... .......... .......... 69% 73.5M 3s
    ##  94050K .......... .......... .......... .......... .......... 69% 13.8M 3s
    ##  94100K .......... .......... .......... .......... .......... 69% 60.9M 3s
    ##  94150K .......... .......... .......... .......... .......... 69% 13.8M 3s
    ##  94200K .......... .......... .......... .......... .......... 69% 6.40M 3s
    ##  94250K .......... .......... .......... .......... .......... 69% 65.1M 3s
    ##  94300K .......... .......... .......... .......... .......... 70% 72.2M 3s
    ##  94350K .......... .......... .......... .......... .......... 70% 79.1M 3s
    ##  94400K .......... .......... .......... .......... .......... 70% 77.7M 3s
    ##  94450K .......... .......... .......... .......... .......... 70% 68.8M 3s
    ##  94500K .......... .......... .......... .......... .......... 70% 69.8M 3s
    ##  94550K .......... .......... .......... .......... .......... 70% 59.7M 3s
    ##  94600K .......... .......... .......... .......... .......... 70% 74.8M 3s
    ##  94650K .......... .......... .......... .......... .......... 70% 83.9M 3s
    ##  94700K .......... .......... .......... .......... .......... 70% 76.2M 3s
    ##  94750K .......... .......... .......... .......... .......... 70% 66.1M 3s
    ##  94800K .......... .......... .......... .......... .......... 70% 57.3M 3s
    ##  94850K .......... .......... .......... .......... .......... 70% 11.8M 3s
    ##  94900K .......... .......... .......... .......... .......... 70% 57.5M 3s
    ##  94950K .......... .......... .......... .......... .......... 70% 83.6M 3s
    ##  95000K .......... .......... .......... .......... .......... 70% 16.8M 3s
    ##  95050K .......... .......... .......... .......... .......... 70% 69.0M 3s
    ##  95100K .......... .......... .......... .......... .......... 70% 33.6M 3s
    ##  95150K .......... .......... .......... .......... .......... 70% 20.4M 3s
    ##  95200K .......... .......... .......... .......... .......... 70% 65.9M 3s
    ##  95250K .......... .......... .......... .......... .......... 70% 22.9M 3s
    ##  95300K .......... .......... .......... .......... .......... 70% 23.5M 3s
    ##  95350K .......... .......... .......... .......... .......... 70% 90.8M 3s
    ##  95400K .......... .......... .......... .......... .......... 70% 16.4M 3s
    ##  95450K .......... .......... .......... .......... .......... 70% 56.1M 3s
    ##  95500K .......... .......... .......... .......... .......... 70% 20.6M 3s
    ##  95550K .......... .......... .......... .......... .......... 70% 33.1M 3s
    ##  95600K .......... .......... .......... .......... .......... 70% 27.6M 3s
    ##  95650K .......... .......... .......... .......... .......... 71% 26.0M 3s
    ##  95700K .......... .......... .......... .......... .......... 71% 22.9M 2s
    ##  95750K .......... .......... .......... .......... .......... 71% 20.8M 2s
    ##  95800K .......... .......... .......... .......... .......... 71% 17.6M 2s
    ##  95850K .......... .......... .......... .......... .......... 71% 85.1M 2s
    ##  95900K .......... .......... .......... .......... .......... 71% 19.1M 2s
    ##  95950K .......... .......... .......... .......... .......... 71% 14.5M 2s
    ##  96000K .......... .......... .......... .......... .......... 71% 54.9M 2s
    ##  96050K .......... .......... .......... .......... .......... 71% 36.7M 2s
    ##  96100K .......... .......... .......... .......... .......... 71% 19.0M 2s
    ##  96150K .......... .......... .......... .......... .......... 71% 81.1M 2s
    ##  96200K .......... .......... .......... .......... .......... 71% 16.8M 2s
    ##  96250K .......... .......... .......... .......... .......... 71% 59.0M 2s
    ##  96300K .......... .......... .......... .......... .......... 71% 40.2M 2s
    ##  96350K .......... .......... .......... .......... .......... 71% 19.5M 2s
    ##  96400K .......... .......... .......... .......... .......... 71% 68.8M 2s
    ##  96450K .......... .......... .......... .......... .......... 71% 14.7M 2s
    ##  96500K .......... .......... .......... .......... .......... 71% 44.2M 2s
    ##  96550K .......... .......... .......... .......... .......... 71% 24.2M 2s
    ##  96600K .......... .......... .......... .......... .......... 71% 33.4M 2s
    ##  96650K .......... .......... .......... .......... .......... 71% 51.8M 2s
    ##  96700K .......... .......... .......... .......... .......... 71% 34.4M 2s
    ##  96750K .......... .......... .......... .......... .......... 71% 22.0M 2s
    ##  96800K .......... .......... .......... .......... .......... 71% 12.5M 2s
    ##  96850K .......... .......... .......... .......... .......... 71% 37.4M 2s
    ##  96900K .......... .......... .......... .......... .......... 71% 43.2M 2s
    ##  96950K .......... .......... .......... .......... .......... 71% 7.29M 2s
    ##  97000K .......... .......... .......... .......... .......... 72% 45.4M 2s
    ##  97050K .......... .......... .......... .......... .......... 72% 56.0M 2s
    ##  97100K .......... .......... .......... .......... .......... 72% 3.78M 2s
    ##  97150K .......... .......... .......... .......... .......... 72% 52.6M 2s
    ##  97200K .......... .......... .......... .......... .......... 72% 46.6M 2s
    ##  97250K .......... .......... .......... .......... .......... 72% 51.7M 2s
    ##  97300K .......... .......... .......... .......... .......... 72% 54.1M 2s
    ##  97350K .......... .......... .......... .......... .......... 72% 50.4M 2s
    ##  97400K .......... .......... .......... .......... .......... 72% 48.7M 2s
    ##  97450K .......... .......... .......... .......... .......... 72% 46.2M 2s
    ##  97500K .......... .......... .......... .......... .......... 72% 51.6M 2s
    ##  97550K .......... .......... .......... .......... .......... 72% 47.1M 2s
    ##  97600K .......... .......... .......... .......... .......... 72% 56.2M 2s
    ##  97650K .......... .......... .......... .......... .......... 72% 19.7M 2s
    ##  97700K .......... .......... .......... .......... .......... 72% 41.0M 2s
    ##  97750K .......... .......... .......... .......... .......... 72% 19.8M 2s
    ##  97800K .......... .......... .......... .......... .......... 72% 40.8M 2s
    ##  97850K .......... .......... .......... .......... .......... 72% 17.3M 2s
    ##  97900K .......... .......... .......... .......... .......... 72% 38.1M 2s
    ##  97950K .......... .......... .......... .......... .......... 72% 17.5M 2s
    ##  98000K .......... .......... .......... .......... .......... 72% 45.1M 2s
    ##  98050K .......... .......... .......... .......... .......... 72% 14.6M 2s
    ##  98100K .......... .......... .......... .......... .......... 72% 47.8M 2s
    ##  98150K .......... .......... .......... .......... .......... 72% 16.6M 2s
    ##  98200K .......... .......... .......... .......... .......... 72% 51.6M 2s
    ##  98250K .......... .......... .......... .......... .......... 72% 18.2M 2s
    ##  98300K .......... .......... .......... .......... .......... 72% 49.8M 2s
    ##  98350K .......... .......... .......... .......... .......... 73% 15.4M 2s
    ##  98400K .......... .......... .......... .......... .......... 73% 29.1M 2s
    ##  98450K .......... .......... .......... .......... .......... 73% 22.8M 2s
    ##  98500K .......... .......... .......... .......... .......... 73% 32.6M 2s
    ##  98550K .......... .......... .......... .......... .......... 73% 19.3M 2s
    ##  98600K .......... .......... .......... .......... .......... 73% 47.9M 2s
    ##  98650K .......... .......... .......... .......... .......... 73% 24.3M 2s
    ##  98700K .......... .......... .......... .......... .......... 73% 44.8M 2s
    ##  98750K .......... .......... .......... .......... .......... 73% 20.4M 2s
    ##  98800K .......... .......... .......... .......... .......... 73% 35.1M 2s
    ##  98850K .......... .......... .......... .......... .......... 73% 25.4M 2s
    ##  98900K .......... .......... .......... .......... .......... 73% 18.9M 2s
    ##  98950K .......... .......... .......... .......... .......... 73% 37.8M 2s
    ##  99000K .......... .......... .......... .......... .......... 73% 24.1M 2s
    ##  99050K .......... .......... .......... .......... .......... 73% 36.6M 2s
    ##  99100K .......... .......... .......... .......... .......... 73% 21.1M 2s
    ##  99150K .......... .......... .......... .......... .......... 73% 22.6M 2s
    ##  99200K .......... .......... .......... .......... .......... 73% 41.9M 2s
    ##  99250K .......... .......... .......... .......... .......... 73% 36.4M 2s
    ##  99300K .......... .......... .......... .......... .......... 73% 15.9M 2s
    ##  99350K .......... .......... .......... .......... .......... 73% 24.4M 2s
    ##  99400K .......... .......... .......... .......... .......... 73% 40.9M 2s
    ##  99450K .......... .......... .......... .......... .......... 73% 16.4M 2s
    ##  99500K .......... .......... .......... .......... .......... 73% 47.6M 2s
    ##  99550K .......... .......... .......... .......... .......... 73% 21.2M 2s
    ##  99600K .......... .......... .......... .......... .......... 73% 15.6M 2s
    ##  99650K .......... .......... .......... .......... .......... 73% 49.8M 2s
    ##  99700K .......... .......... .......... .......... .......... 74% 22.3M 2s
    ##  99750K .......... .......... .......... .......... .......... 74% 35.7M 2s
    ##  99800K .......... .......... .......... .......... .......... 74% 19.8M 2s
    ##  99850K .......... .......... .......... .......... .......... 74% 42.5M 2s
    ##  99900K .......... .......... .......... .......... .......... 74% 48.3M 2s
    ##  99950K .......... .......... .......... .......... .......... 74% 13.2M 2s
    ## 100000K .......... .......... .......... .......... .......... 74% 53.0M 2s
    ## 100050K .......... .......... .......... .......... .......... 74% 20.4M 2s
    ## 100100K .......... .......... .......... .......... .......... 74% 18.0M 2s
    ## 100150K .......... .......... .......... .......... .......... 74% 33.2M 2s
    ## 100200K .......... .......... .......... .......... .......... 74% 52.5M 2s
    ## 100250K .......... .......... .......... .......... .......... 74% 14.9M 2s
    ## 100300K .......... .......... .......... .......... .......... 74% 49.3M 2s
    ## 100350K .......... .......... .......... .......... .......... 74% 18.7M 2s
    ## 100400K .......... .......... .......... .......... .......... 74% 42.0M 2s
    ## 100450K .......... .......... .......... .......... .......... 74% 15.7M 2s
    ## 100500K .......... .......... .......... .......... .......... 74% 46.7M 2s
    ## 100550K .......... .......... .......... .......... .......... 74% 14.5M 2s
    ## 100600K .......... .......... .......... .......... .......... 74% 46.8M 2s
    ## 100650K .......... .......... .......... .......... .......... 74% 63.3M 2s
    ## 100700K .......... .......... .......... .......... .......... 74% 8.86M 2s
    ## 100750K .......... .......... .......... .......... .......... 74% 21.2M 2s
    ## 100800K .......... .......... .......... .......... .......... 74% 13.7M 2s
    ## 100850K .......... .......... .......... .......... .......... 74% 51.1M 2s
    ## 100900K .......... .......... .......... .......... .......... 74% 12.6M 2s
    ## 100950K .......... .......... .......... .......... .......... 74% 47.6M 2s
    ## 101000K .......... .......... .......... .......... .......... 74% 19.3M 2s
    ## 101050K .......... .......... .......... .......... .......... 75% 15.1M 2s
    ## 101100K .......... .......... .......... .......... .......... 75% 37.0M 2s
    ## 101150K .......... .......... .......... .......... .......... 75% 55.2M 2s
    ## 101200K .......... .......... .......... .......... .......... 75% 19.7M 2s
    ## 101250K .......... .......... .......... .......... .......... 75% 18.5M 2s
    ## 101300K .......... .......... .......... .......... .......... 75% 43.2M 2s
    ## 101350K .......... .......... .......... .......... .......... 75% 42.7M 2s
    ## 101400K .......... .......... .......... .......... .......... 75% 15.7M 2s
    ## 101450K .......... .......... .......... .......... .......... 75% 36.3M 2s
    ## 101500K .......... .......... .......... .......... .......... 75% 6.34M 2s
    ## 101550K .......... .......... .......... .......... .......... 75% 50.1M 2s
    ## 101600K .......... .......... .......... .......... .......... 75% 56.0M 2s
    ## 101650K .......... .......... .......... .......... .......... 75% 57.8M 2s
    ## 101700K .......... .......... .......... .......... .......... 75% 49.8M 2s
    ## 101750K .......... .......... .......... .......... .......... 75% 54.0M 2s
    ## 101800K .......... .......... .......... .......... .......... 75% 53.2M 2s
    ## 101850K .......... .......... .......... .......... .......... 75% 5.37M 2s
    ## 101900K .......... .......... .......... .......... .......... 75% 11.7M 2s
    ## 101950K .......... .......... .......... .......... .......... 75% 14.6M 2s
    ## 102000K .......... .......... .......... .......... .......... 75% 44.9M 2s
    ## 102050K .......... .......... .......... .......... .......... 75% 15.6M 2s
    ## 102100K .......... .......... .......... .......... .......... 75% 13.0M 2s
    ## 102150K .......... .......... .......... .......... .......... 75% 14.2M 2s
    ## 102200K .......... .......... .......... .......... .......... 75% 33.4M 2s
    ## 102250K .......... .......... .......... .......... .......... 75% 13.5M 2s
    ## 102300K .......... .......... .......... .......... .......... 75% 13.2M 2s
    ## 102350K .......... .......... .......... .......... .......... 75% 44.0M 2s
    ## 102400K .......... .......... .......... .......... .......... 76% 15.8M 2s
    ## 102450K .......... .......... .......... .......... .......... 76% 49.2M 2s
    ## 102500K .......... .......... .......... .......... .......... 76% 13.7M 2s
    ## 102550K .......... .......... .......... .......... .......... 76% 13.6M 2s
    ## 102600K .......... .......... .......... .......... .......... 76% 13.0M 2s
    ## 102650K .......... .......... .......... .......... .......... 76% 13.9M 2s
    ## 102700K .......... .......... .......... .......... .......... 76% 46.4M 2s
    ## 102750K .......... .......... .......... .......... .......... 76% 8.50M 2s
    ## 102800K .......... .......... .......... .......... .......... 76% 45.8M 2s
    ## 102850K .......... .......... .......... .......... .......... 76% 4.86M 2s
    ## 102900K .......... .......... .......... .......... .......... 76% 48.3M 2s
    ## 102950K .......... .......... .......... .......... .......... 76% 12.5M 2s
    ## 103000K .......... .......... .......... .......... .......... 76% 49.7M 2s
    ## 103050K .......... .......... .......... .......... .......... 76% 13.7M 2s
    ## 103100K .......... .......... .......... .......... .......... 76% 50.1M 2s
    ## 103150K .......... .......... .......... .......... .......... 76% 10.4M 2s
    ## 103200K .......... .......... .......... .......... .......... 76% 50.2M 2s
    ## 103250K .......... .......... .......... .......... .......... 76% 14.4M 2s
    ## 103300K .......... .......... .......... .......... .......... 76% 49.5M 2s
    ## 103350K .......... .......... .......... .......... .......... 76% 17.4M 2s
    ## 103400K .......... .......... .......... .......... .......... 76% 47.9M 2s
    ## 103450K .......... .......... .......... .......... .......... 76% 12.0M 2s
    ## 103500K .......... .......... .......... .......... .......... 76% 14.1M 2s
    ## 103550K .......... .......... .......... .......... .......... 76% 43.2M 2s
    ## 103600K .......... .......... .......... .......... .......... 76% 16.4M 2s
    ## 103650K .......... .......... .......... .......... .......... 76% 56.8M 2s
    ## 103700K .......... .......... .......... .......... .......... 77% 14.4M 2s
    ## 103750K .......... .......... .......... .......... .......... 77% 42.3M 2s
    ## 103800K .......... .......... .......... .......... .......... 77% 17.9M 2s
    ## 103850K .......... .......... .......... .......... .......... 77% 43.7M 2s
    ## 103900K .......... .......... .......... .......... .......... 77% 13.6M 2s
    ## 103950K .......... .......... .......... .......... .......... 77% 58.6M 2s
    ## 104000K .......... .......... .......... .......... .......... 77% 7.80M 2s
    ## 104050K .......... .......... .......... .......... .......... 77% 56.7M 2s
    ## 104100K .......... .......... .......... .......... .......... 77% 15.5M 2s
    ## 104150K .......... .......... .......... .......... .......... 77% 48.7M 2s
    ## 104200K .......... .......... .......... .......... .......... 77% 14.9M 2s
    ## 104250K .......... .......... .......... .......... .......... 77% 58.2M 2s
    ## 104300K .......... .......... .......... .......... .......... 77% 18.9M 2s
    ## 104350K .......... .......... .......... .......... .......... 77% 49.8M 2s
    ## 104400K .......... .......... .......... .......... .......... 77% 16.9M 2s
    ## 104450K .......... .......... .......... .......... .......... 77% 28.5M 2s
    ## 104500K .......... .......... .......... .......... .......... 77% 13.0M 2s
    ## 104550K .......... .......... .......... .......... .......... 77% 10.5M 2s
    ## 104600K .......... .......... .......... .......... .......... 77% 41.0M 2s
    ## 104650K .......... .......... .......... .......... .......... 77% 58.3M 2s
    ## 104700K .......... .......... .......... .......... .......... 77% 16.9M 2s
    ## 104750K .......... .......... .......... .......... .......... 77% 18.7M 2s
    ## 104800K .......... .......... .......... .......... .......... 77% 40.1M 2s
    ## 104850K .......... .......... .......... .......... .......... 77% 28.7M 2s
    ## 104900K .......... .......... .......... .......... .......... 77% 20.7M 2s
    ## 104950K .......... .......... .......... .......... .......... 77% 24.2M 2s
    ## 105000K .......... .......... .......... .......... .......... 77% 22.6M 2s
    ## 105050K .......... .......... .......... .......... .......... 78% 51.4M 2s
    ## 105100K .......... .......... .......... .......... .......... 78% 23.9M 2s
    ## 105150K .......... .......... .......... .......... .......... 78% 18.8M 2s
    ## 105200K .......... .......... .......... .......... .......... 78% 29.9M 2s
    ## 105250K .......... .......... .......... .......... .......... 78% 19.1M 2s
    ## 105300K .......... .......... .......... .......... .......... 78% 24.3M 2s
    ## 105350K .......... .......... .......... .......... .......... 78% 44.0M 2s
    ## 105400K .......... .......... .......... .......... .......... 78% 14.1M 2s
    ## 105450K .......... .......... .......... .......... .......... 78% 46.5M 2s
    ## 105500K .......... .......... .......... .......... .......... 78% 18.3M 2s
    ## 105550K .......... .......... .......... .......... .......... 78% 55.1M 2s
    ## 105600K .......... .......... .......... .......... .......... 78% 16.3M 2s
    ## 105650K .......... .......... .......... .......... .......... 78% 50.9M 2s
    ## 105700K .......... .......... .......... .......... .......... 78% 14.4M 2s
    ## 105750K .......... .......... .......... .......... .......... 78% 55.8M 2s
    ## 105800K .......... .......... .......... .......... .......... 78% 68.1M 2s
    ## 105850K .......... .......... .......... .......... .......... 78% 15.2M 2s
    ## 105900K .......... .......... .......... .......... .......... 78% 55.0M 2s
    ## 105950K .......... .......... .......... .......... .......... 78% 20.9M 2s
    ## 106000K .......... .......... .......... .......... .......... 78% 68.3M 2s
    ## 106050K .......... .......... .......... .......... .......... 78% 15.8M 2s
    ## 106100K .......... .......... .......... .......... .......... 78% 70.1M 2s
    ## 106150K .......... .......... .......... .......... .......... 78% 15.2M 2s
    ## 106200K .......... .......... .......... .......... .......... 78% 16.6M 2s
    ## 106250K .......... .......... .......... .......... .......... 78% 54.7M 2s
    ## 106300K .......... .......... .......... .......... .......... 78% 64.8M 2s
    ## 106350K .......... .......... .......... .......... .......... 78% 12.7M 2s
    ## 106400K .......... .......... .......... .......... .......... 79% 65.1M 2s
    ## 106450K .......... .......... .......... .......... .......... 79% 18.3M 2s
    ## 106500K .......... .......... .......... .......... .......... 79% 53.8M 2s
    ## 106550K .......... .......... .......... .......... .......... 79% 73.4M 2s
    ## 106600K .......... .......... .......... .......... .......... 79% 11.5M 2s
    ## 106650K .......... .......... .......... .......... .......... 79% 70.1M 2s
    ## 106700K .......... .......... .......... .......... .......... 79% 18.9M 2s
    ## 106750K .......... .......... .......... .......... .......... 79% 75.0M 2s
    ## 106800K .......... .......... .......... .......... .......... 79% 24.1M 2s
    ## 106850K .......... .......... .......... .......... .......... 79% 27.4M 2s
    ## 106900K .......... .......... .......... .......... .......... 79% 30.1M 2s
    ## 106950K .......... .......... .......... .......... .......... 79% 30.4M 2s
    ## 107000K .......... .......... .......... .......... .......... 79% 33.7M 2s
    ## 107050K .......... .......... .......... .......... .......... 79% 40.0M 2s
    ## 107100K .......... .......... .......... .......... .......... 79% 20.7M 2s
    ## 107150K .......... .......... .......... .......... .......... 79% 36.7M 2s
    ## 107200K .......... .......... .......... .......... .......... 79% 21.9M 2s
    ## 107250K .......... .......... .......... .......... .......... 79% 58.3M 2s
    ## 107300K .......... .......... .......... .......... .......... 79% 58.3M 2s
    ## 107350K .......... .......... .......... .......... .......... 79% 17.7M 2s
    ## 107400K .......... .......... .......... .......... .......... 79% 83.0M 2s
    ## 107450K .......... .......... .......... .......... .......... 79% 24.6M 2s
    ## 107500K .......... .......... .......... .......... .......... 79% 42.8M 2s
    ## 107550K .......... .......... .......... .......... .......... 79% 17.7M 2s
    ## 107600K .......... .......... .......... .......... .......... 79% 50.5M 2s
    ## 107650K .......... .......... .......... .......... .......... 79% 12.7M 2s
    ## 107700K .......... .......... .......... .......... .......... 79% 53.6M 2s
    ## 107750K .......... .......... .......... .......... .......... 80% 91.2M 2s
    ## 107800K .......... .......... .......... .......... .......... 80% 18.5M 2s
    ## 107850K .......... .......... .......... .......... .......... 80% 71.8M 2s
    ## 107900K .......... .......... .......... .......... .......... 80% 27.0M 2s
    ## 107950K .......... .......... .......... .......... .......... 80% 23.2M 2s
    ## 108000K .......... .......... .......... .......... .......... 80% 80.1M 2s
    ## 108050K .......... .......... .......... .......... .......... 80% 19.9M 2s
    ## 108100K .......... .......... .......... .......... .......... 80% 58.5M 2s
    ## 108150K .......... .......... .......... .......... .......... 80% 46.9M 2s
    ## 108200K .......... .......... .......... .......... .......... 80% 18.6M 2s
    ## 108250K .......... .......... .......... .......... .......... 80% 60.2M 2s
    ## 108300K .......... .......... .......... .......... .......... 80% 17.0M 2s
    ## 108350K .......... .......... .......... .......... .......... 80% 51.4M 2s
    ## 108400K .......... .......... .......... .......... .......... 80% 86.8M 2s
    ## 108450K .......... .......... .......... .......... .......... 80% 18.0M 2s
    ## 108500K .......... .......... .......... .......... .......... 80% 52.1M 2s
    ## 108550K .......... .......... .......... .......... .......... 80% 98.3M 2s
    ## 108600K .......... .......... .......... .......... .......... 80% 18.2M 2s
    ## 108650K .......... .......... .......... .......... .......... 80% 76.6M 2s
    ## 108700K .......... .......... .......... .......... .......... 80% 19.2M 2s
    ## 108750K .......... .......... .......... .......... .......... 80% 41.5M 2s
    ## 108800K .......... .......... .......... .......... .......... 80% 21.6M 2s
    ## 108850K .......... .......... .......... .......... .......... 80% 85.8M 2s
    ## 108900K .......... .......... .......... .......... .......... 80% 26.3M 2s
    ## 108950K .......... .......... .......... .......... .......... 80% 15.4M 2s
    ## 109000K .......... .......... .......... .......... .......... 80% 40.4M 2s
    ## 109050K .......... .......... .......... .......... .......... 80% 21.0M 2s
    ## 109100K .......... .......... .......... .......... .......... 81% 41.3M 2s
    ## 109150K .......... .......... .......... .......... .......... 81% 50.1M 2s
    ## 109200K .......... .......... .......... .......... .......... 81% 35.7M 2s
    ## 109250K .......... .......... .......... .......... .......... 81% 20.6M 2s
    ## 109300K .......... .......... .......... .......... .......... 81% 44.7M 2s
    ## 109350K .......... .......... .......... .......... .......... 81% 51.5M 2s
    ## 109400K .......... .......... .......... .......... .......... 81% 21.1M 2s
    ## 109450K .......... .......... .......... .......... .......... 81% 22.4M 2s
    ## 109500K .......... .......... .......... .......... .......... 81% 16.3M 2s
    ## 109550K .......... .......... .......... .......... .......... 81% 49.1M 2s
    ## 109600K .......... .......... .......... .......... .......... 81% 49.9M 2s
    ## 109650K .......... .......... .......... .......... .......... 81% 13.4M 2s
    ## 109700K .......... .......... .......... .......... .......... 81% 49.3M 2s
    ## 109750K .......... .......... .......... .......... .......... 81% 58.4M 2s
    ## 109800K .......... .......... .......... .......... .......... 81% 22.9M 2s
    ## 109850K .......... .......... .......... .......... .......... 81% 56.2M 2s
    ## 109900K .......... .......... .......... .......... .......... 81% 24.6M 2s
    ## 109950K .......... .......... .......... .......... .......... 81% 46.4M 2s
    ## 110000K .......... .......... .......... .......... .......... 81% 51.8M 2s
    ## 110050K .......... .......... .......... .......... .......... 81% 15.2M 2s
    ## 110100K .......... .......... .......... .......... .......... 81% 49.0M 2s
    ## 110150K .......... .......... .......... .......... .......... 81% 20.2M 1s
    ## 110200K .......... .......... .......... .......... .......... 81% 44.7M 1s
    ## 110250K .......... .......... .......... .......... .......... 81% 51.9M 1s
    ## 110300K .......... .......... .......... .......... .......... 81% 22.4M 1s
    ## 110350K .......... .......... .......... .......... .......... 81% 37.2M 1s
    ## 110400K .......... .......... .......... .......... .......... 81% 20.1M 1s
    ## 110450K .......... .......... .......... .......... .......... 82% 49.4M 1s
    ## 110500K .......... .......... .......... .......... .......... 82% 49.4M 1s
    ## 110550K .......... .......... .......... .......... .......... 82% 17.4M 1s
    ## 110600K .......... .......... .......... .......... .......... 82% 53.4M 1s
    ## 110650K .......... .......... .......... .......... .......... 82% 29.8M 1s
    ## 110700K .......... .......... .......... .......... .......... 82% 37.4M 1s
    ## 110750K .......... .......... .......... .......... .......... 82% 56.4M 1s
    ## 110800K .......... .......... .......... .......... .......... 82% 18.7M 1s
    ## 110850K .......... .......... .......... .......... .......... 82% 47.7M 1s
    ## 110900K .......... .......... .......... .......... .......... 82% 21.2M 1s
    ## 110950K .......... .......... .......... .......... .......... 82% 29.6M 1s
    ## 111000K .......... .......... .......... .......... .......... 82% 46.5M 1s
    ## 111050K .......... .......... .......... .......... .......... 82% 23.5M 1s
    ## 111100K .......... .......... .......... .......... .......... 82% 41.3M 1s
    ## 111150K .......... .......... .......... .......... .......... 82% 17.7M 1s
    ## 111200K .......... .......... .......... .......... .......... 82% 42.0M 1s
    ## 111250K .......... .......... .......... .......... .......... 82% 46.5M 1s
    ## 111300K .......... .......... .......... .......... .......... 82% 27.4M 1s
    ## 111350K .......... .......... .......... .......... .......... 82% 43.1M 1s
    ## 111400K .......... .......... .......... .......... .......... 82% 19.3M 1s
    ## 111450K .......... .......... .......... .......... .......... 82% 45.4M 1s
    ## 111500K .......... .......... .......... .......... .......... 82% 48.1M 1s
    ## 111550K .......... .......... .......... .......... .......... 82% 53.5M 1s
    ## 111600K .......... .......... .......... .......... .......... 82% 32.6M 1s
    ## 111650K .......... .......... .......... .......... .......... 82% 38.9M 1s
    ## 111700K .......... .......... .......... .......... .......... 82% 33.3M 1s
    ## 111750K .......... .......... .......... .......... .......... 82% 33.8M 1s
    ## 111800K .......... .......... .......... .......... .......... 83% 32.4M 1s
    ## 111850K .......... .......... .......... .......... .......... 83% 43.5M 1s
    ## 111900K .......... .......... .......... .......... .......... 83% 14.5M 1s
    ## 111950K .......... .......... .......... .......... .......... 83% 47.0M 1s
    ## 112000K .......... .......... .......... .......... .......... 83% 6.37M 1s
    ## 112050K .......... .......... .......... .......... .......... 83% 39.7M 1s
    ## 112100K .......... .......... .......... .......... .......... 83% 44.6M 1s
    ## 112150K .......... .......... .......... .......... .......... 83% 50.4M 1s
    ## 112200K .......... .......... .......... .......... .......... 83% 43.8M 1s
    ## 112250K .......... .......... .......... .......... .......... 83% 43.2M 1s
    ## 112300K .......... .......... .......... .......... .......... 83% 47.7M 1s
    ## 112350K .......... .......... .......... .......... .......... 83% 43.6M 1s
    ## 112400K .......... .......... .......... .......... .......... 83% 44.6M 1s
    ## 112450K .......... .......... .......... .......... .......... 83% 58.9M 1s
    ## 112500K .......... .......... .......... .......... .......... 83% 29.3M 1s
    ## 112550K .......... .......... .......... .......... .......... 83% 27.9M 1s
    ## 112600K .......... .......... .......... .......... .......... 83% 37.6M 1s
    ## 112650K .......... .......... .......... .......... .......... 83% 17.7M 1s
    ## 112700K .......... .......... .......... .......... .......... 83% 41.5M 1s
    ## 112750K .......... .......... .......... .......... .......... 83% 20.3M 1s
    ## 112800K .......... .......... .......... .......... .......... 83% 15.0M 1s
    ## 112850K .......... .......... .......... .......... .......... 83% 52.2M 1s
    ## 112900K .......... .......... .......... .......... .......... 83% 54.7M 1s
    ## 112950K .......... .......... .......... .......... .......... 83% 16.5M 1s
    ## 113000K .......... .......... .......... .......... .......... 83% 13.8M 1s
    ## 113050K .......... .......... .......... .......... .......... 83% 51.9M 1s
    ## 113100K .......... .......... .......... .......... .......... 83% 19.8M 1s
    ## 113150K .......... .......... .......... .......... .......... 84% 50.5M 1s
    ## 113200K .......... .......... .......... .......... .......... 84% 19.5M 1s
    ## 113250K .......... .......... .......... .......... .......... 84% 31.8M 1s
    ## 113300K .......... .......... .......... .......... .......... 84% 14.2M 1s
    ## 113350K .......... .......... .......... .......... .......... 84% 47.4M 1s
    ## 113400K .......... .......... .......... .......... .......... 84% 51.1M 1s
    ## 113450K .......... .......... .......... .......... .......... 84% 16.9M 1s
    ## 113500K .......... .......... .......... .......... .......... 84% 7.66M 1s
    ## 113550K .......... .......... .......... .......... .......... 84% 50.5M 1s
    ## 113600K .......... .......... .......... .......... .......... 84% 55.1M 1s
    ## 113650K .......... .......... .......... .......... .......... 84% 13.6M 1s
    ## 113700K .......... .......... .......... .......... .......... 84% 17.7M 1s
    ## 113750K .......... .......... .......... .......... .......... 84% 56.1M 1s
    ## 113800K .......... .......... .......... .......... .......... 84% 23.0M 1s
    ## 113850K .......... .......... .......... .......... .......... 84% 25.1M 1s
    ## 113900K .......... .......... .......... .......... .......... 84% 56.1M 1s
    ## 113950K .......... .......... .......... .......... .......... 84% 18.2M 1s
    ## 114000K .......... .......... .......... .......... .......... 84% 55.3M 1s
    ## 114050K .......... .......... .......... .......... .......... 84% 13.8M 1s
    ## 114100K .......... .......... .......... .......... .......... 84% 55.1M 1s
    ## 114150K .......... .......... .......... .......... .......... 84% 12.5M 1s
    ## 114200K .......... .......... .......... .......... .......... 84% 58.7M 1s
    ## 114250K .......... .......... .......... .......... .......... 84% 15.1M 1s
    ## 114300K .......... .......... .......... .......... .......... 84% 61.6M 1s
    ## 114350K .......... .......... .......... .......... .......... 84% 17.9M 1s
    ## 114400K .......... .......... .......... .......... .......... 84% 54.0M 1s
    ## 114450K .......... .......... .......... .......... .......... 84% 10.5M 1s
    ## 114500K .......... .......... .......... .......... .......... 85% 56.3M 1s
    ## 114550K .......... .......... .......... .......... .......... 85% 15.2M 1s
    ## 114600K .......... .......... .......... .......... .......... 85% 53.0M 1s
    ## 114650K .......... .......... .......... .......... .......... 85% 12.3M 1s
    ## 114700K .......... .......... .......... .......... .......... 85% 58.2M 1s
    ## 114750K .......... .......... .......... .......... .......... 85% 14.5M 1s
    ## 114800K .......... .......... .......... .......... .......... 85% 48.9M 1s
    ## 114850K .......... .......... .......... .......... .......... 85% 36.3M 1s
    ## 114900K .......... .......... .......... .......... .......... 85% 21.0M 1s
    ## 114950K .......... .......... .......... .......... .......... 85% 48.0M 1s
    ## 115000K .......... .......... .......... .......... .......... 85% 16.0M 1s
    ## 115050K .......... .......... .......... .......... .......... 85% 25.3M 1s
    ## 115100K .......... .......... .......... .......... .......... 85% 18.2M 1s
    ## 115150K .......... .......... .......... .......... .......... 85% 45.2M 1s
    ## 115200K .......... .......... .......... .......... .......... 85% 10.4M 1s
    ## 115250K .......... .......... .......... .......... .......... 85% 44.1M 1s
    ## 115300K .......... .......... .......... .......... .......... 85% 17.4M 1s
    ## 115350K .......... .......... .......... .......... .......... 85% 50.9M 1s
    ## 115400K .......... .......... .......... .......... .......... 85% 9.99M 1s
    ## 115450K .......... .......... .......... .......... .......... 85% 53.9M 1s
    ## 115500K .......... .......... .......... .......... .......... 85% 55.9M 1s
    ## 115550K .......... .......... .......... .......... .......... 85% 17.7M 1s
    ## 115600K .......... .......... .......... .......... .......... 85% 40.6M 1s
    ## 115650K .......... .......... .......... .......... .......... 85% 17.8M 1s
    ## 115700K .......... .......... .......... .......... .......... 85% 38.4M 1s
    ## 115750K .......... .......... .......... .......... .......... 85% 60.9M 1s
    ## 115800K .......... .......... .......... .......... .......... 85% 17.5M 1s
    ## 115850K .......... .......... .......... .......... .......... 86% 56.6M 1s
    ## 115900K .......... .......... .......... .......... .......... 86% 8.24M 1s
    ## 115950K .......... .......... .......... .......... .......... 86% 49.1M 1s
    ## 116000K .......... .......... .......... .......... .......... 86% 32.5M 1s
    ## 116050K .......... .......... .......... .......... .......... 86% 34.5M 1s
    ## 116100K .......... .......... .......... .......... .......... 86% 19.8M 1s
    ## 116150K .......... .......... .......... .......... .......... 86% 44.5M 1s
    ## 116200K .......... .......... .......... .......... .......... 86% 33.7M 1s
    ## 116250K .......... .......... .......... .......... .......... 86% 28.4M 1s
    ## 116300K .......... .......... .......... .......... .......... 86% 34.1M 1s
    ## 116350K .......... .......... .......... .......... .......... 86% 23.2M 1s
    ## 116400K .......... .......... .......... .......... .......... 86% 42.2M 1s
    ## 116450K .......... .......... .......... .......... .......... 86% 45.2M 1s
    ## 116500K .......... .......... .......... .......... .......... 86% 23.1M 1s
    ## 116550K .......... .......... .......... .......... .......... 86% 33.1M 1s
    ## 116600K .......... .......... .......... .......... .......... 86% 47.2M 1s
    ## 116650K .......... .......... .......... .......... .......... 86% 46.6M 1s
    ## 116700K .......... .......... .......... .......... .......... 86% 21.0M 1s
    ## 116750K .......... .......... .......... .......... .......... 86% 49.7M 1s
    ## 116800K .......... .......... .......... .......... .......... 86% 17.8M 1s
    ## 116850K .......... .......... .......... .......... .......... 86% 29.8M 1s
    ## 116900K .......... .......... .......... .......... .......... 86% 42.5M 1s
    ## 116950K .......... .......... .......... .......... .......... 86% 19.7M 1s
    ## 117000K .......... .......... .......... .......... .......... 86% 24.4M 1s
    ## 117050K .......... .......... .......... .......... .......... 86% 26.3M 1s
    ## 117100K .......... .......... .......... .......... .......... 86% 31.2M 1s
    ## 117150K .......... .......... .......... .......... .......... 86% 48.6M 1s
    ## 117200K .......... .......... .......... .......... .......... 87% 40.9M 1s
    ## 117250K .......... .......... .......... .......... .......... 87% 17.2M 1s
    ## 117300K .......... .......... .......... .......... .......... 87% 30.5M 1s
    ## 117350K .......... .......... .......... .......... .......... 87% 49.2M 1s
    ## 117400K .......... .......... .......... .......... .......... 87% 13.6M 1s
    ## 117450K .......... .......... .......... .......... .......... 87% 15.9M 1s
    ## 117500K .......... .......... .......... .......... .......... 87% 39.3M 1s
    ## 117550K .......... .......... .......... .......... .......... 87% 34.0M 1s
    ## 117600K .......... .......... .......... .......... .......... 87% 28.7M 1s
    ## 117650K .......... .......... .......... .......... .......... 87% 89.7M 1s
    ## 117700K .......... .......... .......... .......... .......... 87% 14.5M 1s
    ## 117750K .......... .......... .......... .......... .......... 87% 76.7M 1s
    ## 117800K .......... .......... .......... .......... .......... 87% 17.3M 1s
    ## 117850K .......... .......... .......... .......... .......... 87% 77.3M 1s
    ## 117900K .......... .......... .......... .......... .......... 87% 58.4M 1s
    ## 117950K .......... .......... .......... .......... .......... 87% 15.4M 1s
    ## 118000K .......... .......... .......... .......... .......... 87% 74.9M 1s
    ## 118050K .......... .......... .......... .......... .......... 87% 88.3M 1s
    ## 118100K .......... .......... .......... .......... .......... 87% 17.1M 1s
    ## 118150K .......... .......... .......... .......... .......... 87% 85.6M 1s
    ## 118200K .......... .......... .......... .......... .......... 87% 25.5M 1s
    ## 118250K .......... .......... .......... .......... .......... 87% 64.1M 1s
    ## 118300K .......... .......... .......... .......... .......... 87% 32.9M 1s
    ## 118350K .......... .......... .......... .......... .......... 87% 24.8M 1s
    ## 118400K .......... .......... .......... .......... .......... 87% 28.6M 1s
    ## 118450K .......... .......... .......... .......... .......... 87% 27.4M 1s
    ## 118500K .......... .......... .......... .......... .......... 87% 79.1M 1s
    ## 118550K .......... .......... .......... .......... .......... 88% 33.8M 1s
    ## 118600K .......... .......... .......... .......... .......... 88% 32.9M 1s
    ## 118650K .......... .......... .......... .......... .......... 88% 23.8M 1s
    ## 118700K .......... .......... .......... .......... .......... 88% 69.8M 1s
    ## 118750K .......... .......... .......... .......... .......... 88% 49.3M 1s
    ## 118800K .......... .......... .......... .......... .......... 88% 22.4M 1s
    ## 118850K .......... .......... .......... .......... .......... 88% 32.5M 1s
    ## 118900K .......... .......... .......... .......... .......... 88% 22.1M 1s
    ## 118950K .......... .......... .......... .......... .......... 88%  103M 1s
    ## 119000K .......... .......... .......... .......... .......... 88% 42.6M 1s
    ## 119050K .......... .......... .......... .......... .......... 88% 25.6M 1s
    ## 119100K .......... .......... .......... .......... .......... 88% 34.2M 1s
    ## 119150K .......... .......... .......... .......... .......... 88% 22.6M 1s
    ## 119200K .......... .......... .......... .......... .......... 88% 31.0M 1s
    ## 119250K .......... .......... .......... .......... .......... 88% 83.0M 1s
    ## 119300K .......... .......... .......... .......... .......... 88% 22.3M 1s
    ## 119350K .......... .......... .......... .......... .......... 88% 34.7M 1s
    ## 119400K .......... .......... .......... .......... .......... 88% 24.0M 1s
    ## 119450K .......... .......... .......... .......... .......... 88% 68.3M 1s
    ## 119500K .......... .......... .......... .......... .......... 88% 29.1M 1s
    ## 119550K .......... .......... .......... .......... .......... 88% 25.2M 1s
    ## 119600K .......... .......... .......... .......... .......... 88% 33.9M 1s
    ## 119650K .......... .......... .......... .......... .......... 88%  102M 1s
    ## 119700K .......... .......... .......... .......... .......... 88% 25.0M 1s
    ## 119750K .......... .......... .......... .......... .......... 88% 25.0M 1s
    ## 119800K .......... .......... .......... .......... .......... 88% 31.4M 1s
    ## 119850K .......... .......... .......... .......... .......... 88% 22.8M 1s
    ## 119900K .......... .......... .......... .......... .......... 89% 32.0M 1s
    ## 119950K .......... .......... .......... .......... .......... 89% 84.3M 1s
    ## 120000K .......... .......... .......... .......... .......... 89% 25.2M 1s
    ## 120050K .......... .......... .......... .......... .......... 89% 35.5M 1s
    ## 120100K .......... .......... .......... .......... .......... 89% 24.2M 1s
    ## 120150K .......... .......... .......... .......... .......... 89% 82.9M 1s
    ## 120200K .......... .......... .......... .......... .......... 89% 44.7M 1s
    ## 120250K .......... .......... .......... .......... .......... 89% 13.1M 1s
    ## 120300K .......... .......... .......... .......... .......... 89% 86.3M 1s
    ## 120350K .......... .......... .......... .......... .......... 89% 39.5M 1s
    ## 120400K .......... .......... .......... .......... .......... 89% 20.1M 1s
    ## 120450K .......... .......... .......... .......... .......... 89% 76.9M 1s
    ## 120500K .......... .......... .......... .......... .......... 89% 26.8M 1s
    ## 120550K .......... .......... .......... .......... .......... 89% 40.7M 1s
    ## 120600K .......... .......... .......... .......... .......... 89% 63.6M 1s
    ## 120650K .......... .......... .......... .......... .......... 89% 24.3M 1s
    ## 120700K .......... .......... .......... .......... .......... 89% 24.7M 1s
    ## 120750K .......... .......... .......... .......... .......... 89% 79.3M 1s
    ## 120800K .......... .......... .......... .......... .......... 89% 29.1M 1s
    ## 120850K .......... .......... .......... .......... .......... 89% 20.2M 1s
    ## 120900K .......... .......... .......... .......... .......... 89% 83.5M 1s
    ## 120950K .......... .......... .......... .......... .......... 89% 36.8M 1s
    ## 121000K .......... .......... .......... .......... .......... 89% 10.3M 1s
    ## 121050K .......... .......... .......... .......... .......... 89% 91.9M 1s
    ## 121100K .......... .......... .......... .......... .......... 89%  105M 1s
    ## 121150K .......... .......... .......... .......... .......... 89% 18.4M 1s
    ## 121200K .......... .......... .......... .......... .......... 89% 88.3M 1s
    ## 121250K .......... .......... .......... .......... .......... 90% 84.7M 1s
    ## 121300K .......... .......... .......... .......... .......... 90% 19.6M 1s
    ## 121350K .......... .......... .......... .......... .......... 90% 81.6M 1s
    ## 121400K .......... .......... .......... .......... .......... 90% 68.8M 1s
    ## 121450K .......... .......... .......... .......... .......... 90% 18.6M 1s
    ## 121500K .......... .......... .......... .......... .......... 90% 61.0M 1s
    ## 121550K .......... .......... .......... .......... .......... 90% 17.4M 1s
    ## 121600K .......... .......... .......... .......... .......... 90% 44.1M 1s
    ## 121650K .......... .......... .......... .......... .......... 90% 93.4M 1s
    ## 121700K .......... .......... .......... .......... .......... 90% 20.2M 1s
    ## 121750K .......... .......... .......... .......... .......... 90% 56.8M 1s
    ## 121800K .......... .......... .......... .......... .......... 90%  100M 1s
    ## 121850K .......... .......... .......... .......... .......... 90% 23.6M 1s
    ## 121900K .......... .......... .......... .......... .......... 90% 45.4M 1s
    ## 121950K .......... .......... .......... .......... .......... 90% 32.4M 1s
    ## 122000K .......... .......... .......... .......... .......... 90% 29.2M 1s
    ## 122050K .......... .......... .......... .......... .......... 90% 70.6M 1s
    ## 122100K .......... .......... .......... .......... .......... 90% 25.4M 1s
    ## 122150K .......... .......... .......... .......... .......... 90% 46.5M 1s
    ## 122200K .......... .......... .......... .......... .......... 90% 19.3M 1s
    ## 122250K .......... .......... .......... .......... .......... 90% 85.8M 1s
    ## 122300K .......... .......... .......... .......... .......... 90% 49.1M 1s
    ## 122350K .......... .......... .......... .......... .......... 90% 21.9M 1s
    ## 122400K .......... .......... .......... .......... .......... 90% 53.5M 1s
    ## 122450K .......... .......... .......... .......... .......... 90% 21.4M 1s
    ## 122500K .......... .......... .......... .......... .......... 90% 64.3M 1s
    ## 122550K .......... .......... .......... .......... .......... 90% 71.0M 1s
    ## 122600K .......... .......... .......... .......... .......... 91% 20.8M 1s
    ## 122650K .......... .......... .......... .......... .......... 91% 39.2M 1s
    ## 122700K .......... .......... .......... .......... .......... 91% 22.1M 1s
    ## 122750K .......... .......... .......... .......... .......... 91% 30.4M 1s
    ## 122800K .......... .......... .......... .......... .......... 91% 94.3M 1s
    ## 122850K .......... .......... .......... .......... .......... 91% 53.1M 1s
    ## 122900K .......... .......... .......... .......... .......... 91% 19.6M 1s
    ## 122950K .......... .......... .......... .......... .......... 91% 67.2M 1s
    ## 123000K .......... .......... .......... .......... .......... 91% 80.1M 1s
    ## 123050K .......... .......... .......... .......... .......... 91% 13.1M 1s
    ## 123100K .......... .......... .......... .......... .......... 91%  113M 1s
    ## 123150K .......... .......... .......... .......... .......... 91% 13.5M 1s
    ## 123200K .......... .......... .......... .......... .......... 91%  115M 1s
    ## 123250K .......... .......... .......... .......... .......... 91%  109M 1s
    ## 123300K .......... .......... .......... .......... .......... 91% 15.0M 1s
    ## 123350K .......... .......... .......... .......... .......... 91%  120M 1s
    ## 123400K .......... .......... .......... .......... .......... 91% 14.2M 1s
    ## 123450K .......... .......... .......... .......... .......... 91% 37.6M 1s
    ## 123500K .......... .......... .......... .......... .......... 91%  110M 1s
    ## 123550K .......... .......... .......... .......... .......... 91% 28.7M 1s
    ## 123600K .......... .......... .......... .......... .......... 91% 49.8M 1s
    ## 123650K .......... .......... .......... .......... .......... 91% 44.5M 1s
    ## 123700K .......... .......... .......... .......... .......... 91% 99.5M 1s
    ## 123750K .......... .......... .......... .......... .......... 91% 38.1M 1s
    ## 123800K .......... .......... .......... .......... .......... 91% 20.9M 1s
    ## 123850K .......... .......... .......... .......... .......... 91% 35.5M 1s
    ## 123900K .......... .......... .......... .......... .......... 91% 21.1M 1s
    ## 123950K .......... .......... .......... .......... .......... 92% 86.5M 1s
    ## 124000K .......... .......... .......... .......... .......... 92% 53.2M 1s
    ## 124050K .......... .......... .......... .......... .......... 92% 20.6M 1s
    ## 124100K .......... .......... .......... .......... .......... 92% 54.1M 1s
    ## 124150K .......... .......... .......... .......... .......... 92% 19.0M 1s
    ## 124200K .......... .......... .......... .......... .......... 92% 97.5M 1s
    ## 124250K .......... .......... .......... .......... .......... 92% 37.7M 1s
    ## 124300K .......... .......... .......... .......... .......... 92% 21.6M 1s
    ## 124350K .......... .......... .......... .......... .......... 92% 57.2M 1s
    ## 124400K .......... .......... .......... .......... .......... 92% 96.9M 1s
    ## 124450K .......... .......... .......... .......... .......... 92% 12.5M 1s
    ## 124500K .......... .......... .......... .......... .......... 92% 31.6M 1s
    ## 124550K .......... .......... .......... .......... .......... 92% 21.2M 1s
    ## 124600K .......... .......... .......... .......... .......... 92% 57.6M 1s
    ## 124650K .......... .......... .......... .......... .......... 92% 66.5M 1s
    ## 124700K .......... .......... .......... .......... .......... 92% 20.6M 1s
    ## 124750K .......... .......... .......... .......... .......... 92% 49.1M 1s
    ## 124800K .......... .......... .......... .......... .......... 92% 19.1M 1s
    ## 124850K .......... .......... .......... .......... .......... 92% 43.7M 1s
    ## 124900K .......... .......... .......... .......... .......... 92% 95.4M 1s
    ## 124950K .......... .......... .......... .......... .......... 92% 16.0M 1s
    ## 125000K .......... .......... .......... .......... .......... 92% 39.2M 1s
    ## 125050K .......... .......... .......... .......... .......... 92% 49.3M 1s
    ## 125100K .......... .......... .......... .......... .......... 92% 30.6M 1s
    ## 125150K .......... .......... .......... .......... .......... 92% 51.1M 1s
    ## 125200K .......... .......... .......... .......... .......... 92% 55.7M 1s
    ## 125250K .......... .......... .......... .......... .......... 92% 17.1M 1s
    ## 125300K .......... .......... .......... .......... .......... 93% 59.3M 1s
    ## 125350K .......... .......... .......... .......... .......... 93% 65.0M 1s
    ## 125400K .......... .......... .......... .......... .......... 93% 14.1M 1s
    ## 125450K .......... .......... .......... .......... .......... 93% 52.7M 1s
    ## 125500K .......... .......... .......... .......... .......... 93% 17.5M 1s
    ## 125550K .......... .......... .......... .......... .......... 93% 58.0M 1s
    ## 125600K .......... .......... .......... .......... .......... 93% 47.0M 1s
    ## 125650K .......... .......... .......... .......... .......... 93% 18.4M 1s
    ## 125700K .......... .......... .......... .......... .......... 93% 41.1M 1s
    ## 125750K .......... .......... .......... .......... .......... 93% 45.3M 1s
    ## 125800K .......... .......... .......... .......... .......... 93% 20.2M 1s
    ## 125850K .......... .......... .......... .......... .......... 93% 27.2M 1s
    ## 125900K .......... .......... .......... .......... .......... 93% 41.9M 1s
    ## 125950K .......... .......... .......... .......... .......... 93% 40.6M 1s
    ## 126000K .......... .......... .......... .......... .......... 93% 24.0M 1s
    ## 126050K .......... .......... .......... .......... .......... 93% 48.6M 0s
    ## 126100K .......... .......... .......... .......... .......... 93% 36.0M 0s
    ## 126150K .......... .......... .......... .......... .......... 93% 41.7M 0s
    ## 126200K .......... .......... .......... .......... .......... 93% 50.7M 0s
    ## 126250K .......... .......... .......... .......... .......... 93% 26.9M 0s
    ## 126300K .......... .......... .......... .......... .......... 93% 41.8M 0s
    ## 126350K .......... .......... .......... .......... .......... 93% 51.5M 0s
    ## 126400K .......... .......... .......... .......... .......... 93% 33.3M 0s
    ## 126450K .......... .......... .......... .......... .......... 93% 39.9M 0s
    ## 126500K .......... .......... .......... .......... .......... 93% 45.4M 0s
    ## 126550K .......... .......... .......... .......... .......... 93% 32.2M 0s
    ## 126600K .......... .......... .......... .......... .......... 93% 24.0M 0s
    ## 126650K .......... .......... .......... .......... .......... 94% 31.7M 0s
    ## 126700K .......... .......... .......... .......... .......... 94% 30.2M 0s
    ## 126750K .......... .......... .......... .......... .......... 94% 63.2M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 29.2M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 26.5M 0s
    ## 126900K .......... .......... .......... .......... .......... 94% 28.0M 0s
    ## 126950K .......... .......... .......... .......... .......... 94% 21.8M 0s
    ## 127000K .......... .......... .......... .......... .......... 94% 42.3M 0s
    ## 127050K .......... .......... .......... .......... .......... 94% 50.8M 0s
    ## 127100K .......... .......... .......... .......... .......... 94% 24.9M 0s
    ## 127150K .......... .......... .......... .......... .......... 94% 43.1M 0s
    ## 127200K .......... .......... .......... .......... .......... 94% 42.5M 0s
    ## 127250K .......... .......... .......... .......... .......... 94% 50.1M 0s
    ## 127300K .......... .......... .......... .......... .......... 94% 26.4M 0s
    ## 127350K .......... .......... .......... .......... .......... 94% 27.5M 0s
    ## 127400K .......... .......... .......... .......... .......... 94% 33.9M 0s
    ## 127450K .......... .......... .......... .......... .......... 94% 53.9M 0s
    ## 127500K .......... .......... .......... .......... .......... 94% 28.2M 0s
    ## 127550K .......... .......... .......... .......... .......... 94% 41.2M 0s
    ## 127600K .......... .......... .......... .......... .......... 94% 23.9M 0s
    ## 127650K .......... .......... .......... .......... .......... 94% 37.9M 0s
    ## 127700K .......... .......... .......... .......... .......... 94% 46.7M 0s
    ## 127750K .......... .......... .......... .......... .......... 94% 31.2M 0s
    ## 127800K .......... .......... .......... .......... .......... 94% 34.6M 0s
    ## 127850K .......... .......... .......... .......... .......... 94% 26.0M 0s
    ## 127900K .......... .......... .......... .......... .......... 94% 34.2M 0s
    ## 127950K .......... .......... .......... .......... .......... 94% 66.9M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 36.4M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 21.9M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 24.8M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 62.2M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 49.3M 0s
    ## 128250K .......... .......... .......... .......... .......... 95% 12.8M 0s
    ## 128300K .......... .......... .......... .......... .......... 95% 42.9M 0s
    ## 128350K .......... .......... .......... .......... .......... 95% 52.6M 0s
    ## 128400K .......... .......... .......... .......... .......... 95% 15.8M 0s
    ## 128450K .......... .......... .......... .......... .......... 95% 47.3M 0s
    ## 128500K .......... .......... .......... .......... .......... 95% 35.3M 0s
    ## 128550K .......... .......... .......... .......... .......... 95% 57.8M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 19.4M 0s
    ## 128650K .......... .......... .......... .......... .......... 95% 44.1M 0s
    ## 128700K .......... .......... .......... .......... .......... 95% 47.7M 0s
    ## 128750K .......... .......... .......... .......... .......... 95% 22.5M 0s
    ## 128800K .......... .......... .......... .......... .......... 95% 36.1M 0s
    ## 128850K .......... .......... .......... .......... .......... 95% 54.4M 0s
    ## 128900K .......... .......... .......... .......... .......... 95% 35.6M 0s
    ## 128950K .......... .......... .......... .......... .......... 95% 41.2M 0s
    ## 129000K .......... .......... .......... .......... .......... 95% 31.6M 0s
    ## 129050K .......... .......... .......... .......... .......... 95% 27.6M 0s
    ## 129100K .......... .......... .......... .......... .......... 95% 29.8M 0s
    ## 129150K .......... .......... .......... .......... .......... 95% 26.9M 0s
    ## 129200K .......... .......... .......... .......... .......... 95% 41.8M 0s
    ## 129250K .......... .......... .......... .......... .......... 95% 55.0M 0s
    ## 129300K .......... .......... .......... .......... .......... 95% 24.0M 0s
    ## 129350K .......... .......... .......... .......... .......... 96% 42.0M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 32.7M 0s
    ## 129450K .......... .......... .......... .......... .......... 96% 41.9M 0s
    ## 129500K .......... .......... .......... .......... .......... 96% 37.3M 0s
    ## 129550K .......... .......... .......... .......... .......... 96% 39.5M 0s
    ## 129600K .......... .......... .......... .......... .......... 96% 41.6M 0s
    ## 129650K .......... .......... .......... .......... .......... 96% 37.3M 0s
    ## 129700K .......... .......... .......... .......... .......... 96% 38.2M 0s
    ## 129750K .......... .......... .......... .......... .......... 96% 41.0M 0s
    ## 129800K .......... .......... .......... .......... .......... 96% 35.4M 0s
    ## 129850K .......... .......... .......... .......... .......... 96% 36.6M 0s
    ## 129900K .......... .......... .......... .......... .......... 96% 33.9M 0s
    ## 129950K .......... .......... .......... .......... .......... 96% 48.5M 0s
    ## 130000K .......... .......... .......... .......... .......... 96% 34.7M 0s
    ## 130050K .......... .......... .......... .......... .......... 96% 38.7M 0s
    ## 130100K .......... .......... .......... .......... .......... 96% 44.2M 0s
    ## 130150K .......... .......... .......... .......... .......... 96% 46.1M 0s
    ## 130200K .......... .......... .......... .......... .......... 96% 54.7M 0s
    ## 130250K .......... .......... .......... .......... .......... 96% 36.2M 0s
    ## 130300K .......... .......... .......... .......... .......... 96% 42.5M 0s
    ## 130350K .......... .......... .......... .......... .......... 96% 52.4M 0s
    ## 130400K .......... .......... .......... .......... .......... 96% 59.0M 0s
    ## 130450K .......... .......... .......... .......... .......... 96% 63.4M 0s
    ## 130500K .......... .......... .......... .......... .......... 96% 46.8M 0s
    ## 130550K .......... .......... .......... .......... .......... 96% 64.0M 0s
    ## 130600K .......... .......... .......... .......... .......... 96% 34.1M 0s
    ## 130650K .......... .......... .......... .......... .......... 97% 49.9M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 47.1M 0s
    ## 130750K .......... .......... .......... .......... .......... 97% 54.4M 0s
    ## 130800K .......... .......... .......... .......... .......... 97% 35.8M 0s
    ## 130850K .......... .......... .......... .......... .......... 97% 53.1M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 59.3M 0s
    ## 130950K .......... .......... .......... .......... .......... 97% 54.5M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 34.1M 0s
    ## 131050K .......... .......... .......... .......... .......... 97% 58.7M 0s
    ## 131100K .......... .......... .......... .......... .......... 97% 38.3M 0s
    ## 131150K .......... .......... .......... .......... .......... 97% 37.2M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 62.1M 0s
    ## 131250K .......... .......... .......... .......... .......... 97% 55.6M 0s
    ## 131300K .......... .......... .......... .......... .......... 97% 40.0M 0s
    ## 131350K .......... .......... .......... .......... .......... 97% 66.8M 0s
    ## 131400K .......... .......... .......... .......... .......... 97% 43.2M 0s
    ## 131450K .......... .......... .......... .......... .......... 97% 53.8M 0s
    ## 131500K .......... .......... .......... .......... .......... 97% 62.2M 0s
    ## 131550K .......... .......... .......... .......... .......... 97% 71.4M 0s
    ## 131600K .......... .......... .......... .......... .......... 97% 68.5M 0s
    ## 131650K .......... .......... .......... .......... .......... 97% 39.6M 0s
    ## 131700K .......... .......... .......... .......... .......... 97% 41.1M 0s
    ## 131750K .......... .......... .......... .......... .......... 97% 64.1M 0s
    ## 131800K .......... .......... .......... .......... .......... 97% 44.6M 0s
    ## 131850K .......... .......... .......... .......... .......... 97% 77.6M 0s
    ## 131900K .......... .......... .......... .......... .......... 97% 25.2M 0s
    ## 131950K .......... .......... .......... .......... .......... 97% 40.9M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 62.3M 0s
    ## 132050K .......... .......... .......... .......... .......... 98% 67.3M 0s
    ## 132100K .......... .......... .......... .......... .......... 98% 51.6M 0s
    ## 132150K .......... .......... .......... .......... .......... 98% 60.1M 0s
    ## 132200K .......... .......... .......... .......... .......... 98% 51.2M 0s
    ## 132250K .......... .......... .......... .......... .......... 98% 79.3M 0s
    ## 132300K .......... .......... .......... .......... .......... 98% 38.5M 0s
    ## 132350K .......... .......... .......... .......... .......... 98% 45.4M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 44.5M 0s
    ## 132450K .......... .......... .......... .......... .......... 98% 68.0M 0s
    ## 132500K .......... .......... .......... .......... .......... 98% 45.8M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 54.3M 0s
    ## 132600K .......... .......... .......... .......... .......... 98% 46.5M 0s
    ## 132650K .......... .......... .......... .......... .......... 98% 62.1M 0s
    ## 132700K .......... .......... .......... .......... .......... 98% 40.3M 0s
    ## 132750K .......... .......... .......... .......... .......... 98% 73.3M 0s
    ## 132800K .......... .......... .......... .......... .......... 98% 39.4M 0s
    ## 132850K .......... .......... .......... .......... .......... 98% 58.6M 0s
    ## 132900K .......... .......... .......... .......... .......... 98% 60.7M 0s
    ## 132950K .......... .......... .......... .......... .......... 98% 54.4M 0s
    ## 133000K .......... .......... .......... .......... .......... 98% 30.7M 0s
    ## 133050K .......... .......... .......... .......... .......... 98% 45.8M 0s
    ## 133100K .......... .......... .......... .......... .......... 98% 57.8M 0s
    ## 133150K .......... .......... .......... .......... .......... 98% 64.1M 0s
    ## 133200K .......... .......... .......... .......... .......... 98% 57.3M 0s
    ## 133250K .......... .......... .......... .......... .......... 98% 33.6M 0s
    ## 133300K .......... .......... .......... .......... .......... 98% 57.3M 0s
    ## 133350K .......... .......... .......... .......... .......... 99% 71.6M 0s
    ## 133400K .......... .......... .......... .......... .......... 99% 40.9M 0s
    ## 133450K .......... .......... .......... .......... .......... 99% 58.8M 0s
    ## 133500K .......... .......... .......... .......... .......... 99% 21.8M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 60.4M 0s
    ## 133600K .......... .......... .......... .......... .......... 99% 34.0M 0s
    ## 133650K .......... .......... .......... .......... .......... 99% 78.9M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 52.6M 0s
    ## 133750K .......... .......... .......... .......... .......... 99% 24.3M 0s
    ## 133800K .......... .......... .......... .......... .......... 99% 49.7M 0s
    ## 133850K .......... .......... .......... .......... .......... 99% 67.0M 0s
    ## 133900K .......... .......... .......... .......... .......... 99% 52.1M 0s
    ## 133950K .......... .......... .......... .......... .......... 99% 30.8M 0s
    ## 134000K .......... .......... .......... .......... .......... 99% 50.8M 0s
    ## 134050K .......... .......... .......... .......... .......... 99% 54.5M 0s
    ## 134100K .......... .......... .......... .......... .......... 99% 70.9M 0s
    ## 134150K .......... .......... .......... .......... .......... 99% 19.8M 0s
    ## 134200K .......... .......... .......... .......... .......... 99% 55.9M 0s
    ## 134250K .......... .......... .......... .......... .......... 99% 77.7M 0s
    ## 134300K .......... .......... .......... .......... .......... 99% 25.4M 0s
    ## 134350K .......... .......... .......... .......... .......... 99% 33.1M 0s
    ## 134400K .......... .......... .......... .......... .......... 99% 43.9M 0s
    ## 134450K .......... .......... .......... .......... .......... 99% 52.2M 0s
    ## 134500K .......... .......... .......... .......... .......... 99% 67.9M 0s
    ## 134550K .......... .......... .......... .......... .......... 99% 60.9M 0s
    ## 134600K .......... .......... .......... .......... .......... 99% 48.5M 0s
    ## 134650K .......... .......... .......... .......... .......... 99% 73.6M 0s
    ## 134700K .......... .......... .......... ..........           100% 22.6M=7.5s
    ## 
    ## 2020-12-04 14:23:08 (17.6 MB/s) - ‘silva_nr99_v138_train_set.fa.gz.8’ saved [137973851/137973851]

Ici on assigne une taxonomie a nos séquences via la fonction
“assignTaxonomy” dans les différents échantillons.C’est une première
façon de visualiser la taxonomie.

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

Ici, nous avons une seconde façon de visualiser la taxonomie. La
fonction “addSpecies” affecte des binômes genre-espèce aux séquences.
Puis il y a une fusion avec le tableau taxinomique. Seul ce qui est
cohérent va être inclus dans le tableau de retour.

``` r
taxa <- addSpecies(taxa, "~/silva_species_assignment_v138.fa.gz")
```

Ici on peut voir les assignements taxonomiques. On supprime les noms de
séquence pour permettre seulement l’affichage.

``` r
taxa.print <- taxa 
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum         Class         Order           Family          
    ## [1,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidota" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      Genus         Species
    ## [1,] NA            NA     
    ## [2,] NA            NA     
    ## [3,] NA            NA     
    ## [4,] NA            NA     
    ## [5,] "Bacteroides" NA     
    ## [6,] NA            NA

On peut voir que les Bacteroidetes sont un des taxons les plus
abondants. Ceci est en corrélation avec le fait qu’il s’agisse
d’échantillons fécaux.

# Evaluer la précision :

## Evalution de la précision de DADA2 sur la communauté estimée :

Ici la fonction “sort” va ordonner des variables. La fonction “cat”
permet de sortir les objets en concaténant les représentions. Ici DADA2
a déduit 20 séquences d’échantillons dans la communauté estimée.

``` r
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
```

    ## DADA2 inferred 20 sample sequences present in the Mock community.

Ici, la communauté estimée contenait 20 souches bactériennes. Comme
DADA2 a identifié 20 ASV le taux d’erreur résiduel est de 0%.

``` r
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
```

    ## Of those, 20 were exact matches to the expected reference sequences.

# Bonus : Transfert à phyloseq :

Installation du package phyloseq.

``` r
BiocManager::install("phyloseq")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.3 (2020-10-10)

    ## Installing package(s) 'phyloseq'

    ## Installation path not writeable, unable to update packages: codetools,
    ##   KernSmooth, nlme

Avec “library” on vérifie que l’on a bien le package phyloseq.

``` r
library(phyloseq); packageVersion("phyloseq")
```

    ## [1] '1.32.0'

A nouveau avec “library” on vérifie que l’on a bien le package
Biostrings".

``` r
library(Biostrings); packageVersion("Biostrings")
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## [1] '2.56.0'

On vérifie que l’on a bien le package ggplot2.

``` r
library(ggplot2); packageVersion("ggplot2")
```

    ## [1] '3.3.2'

La fonction “theme\_set” modifie le theme pour la session R. “theme\_bw”
donne une couleur de fond blanche et une grille grise.

``` r
theme_set(theme_bw())
```

Ici on construit un exemple à partir des informations dans les fichiers.

``` r
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
```

La fonction “phyloseq” permet de construire l’objet ps à partir des
données de DADA2.

``` r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
```

“DNAStringSet” permet de repésenter les chaines de séquences d’ADN. Ceci
permet de donner les objets de classe phyloseq. On a la table des OTU,
des exemples de données, la table de taxonomie et les séquences de
références.

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 232 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 232 taxa by 7 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 232 reference sequences ]

## Visualisez l’alpha-diversité.

La fonction“plot\_richness” estime des indices de diversité alpha grâce
aux indice de Shannon et Simpson et renvoie un graphique.

``` r
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

Sur l’axe des ordonnées on peut voir la mesure de l’alpha diversité et
le nombre de jours sur l’axe des abscisses. On remarque qu’il n’y a pas
de différence entre les échantillons précoces et tardifs.

## Création d’un graphique d’ordination :

La fonction “transform\_sample\_counts” va transformer les comptages
d’échantillons en fonction des OTU. La fonction “ordinate” va
effectuer une ordination sur les données Phyloseq.

``` r
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
```

    ## Run 0 stress 0.08574537 
    ## Run 1 stress 0.08002299 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04283623  max resid 0.1433489 
    ## Run 2 stress 0.08942859 
    ## Run 3 stress 0.09421601 
    ## Run 4 stress 0.09421602 
    ## Run 5 stress 0.08942864 
    ## Run 6 stress 0.1331726 
    ## Run 7 stress 0.08942864 
    ## Run 8 stress 0.08002299 
    ## ... Procrustes: rmse 1.802237e-06  max resid 5.836014e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.09421601 
    ## Run 10 stress 0.08574537 
    ## Run 11 stress 0.09421601 
    ## Run 12 stress 0.08574537 
    ## Run 13 stress 0.08002299 
    ## ... Procrustes: rmse 3.764295e-06  max resid 1.038775e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.08574537 
    ## Run 15 stress 0.1233096 
    ## Run 16 stress 0.09421602 
    ## Run 17 stress 0.08002299 
    ## ... Procrustes: rmse 3.464814e-06  max resid 9.676464e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.08574537 
    ## Run 19 stress 0.08574537 
    ## Run 20 stress 0.08002299 
    ## ... Procrustes: rmse 3.635033e-06  max resid 1.054972e-05 
    ## ... Similar to previous best
    ## *** Solution reached

La fonction “plot\_ordination” va tracer une ordination phyloseq sous
forme de graphique ggplot2.

``` r
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

Ici on peut voir sur le graphique d’ordination qu’il y a une séparation
très distincte entre les échantillons précoces et tardifs.

# Histogrammes :

Ici la fonction “plot\_bar” permet de créer un histogramme, c’est à dire
des graphiques qui résument les différences d’abondance des taxons entre
les échantillons.

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Day", fill="Family") + facet_wrap(~When, scales="free_x")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

Dans ce graphique on peut voir une mince différence dans la distribution
taxonomique pour les échantillons précoces et tardifs. Dans les deux cas
on peut voir une abondance des Muribaculaceae. On peut également voir
globalement une très faible diminution de l’abondance pour les
échantillons précoces. D’autres analyses sont nécessaires pour aller
plus loin.

# Tutoriel phyloseq

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

``` r
library(DECIPHER)
```

    ## Loading required package: RSQLite

``` r
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
```

    ## negative edges length changed to 0!

``` r
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

# Import des données phyloseq du tutoriel phyloseq :

Ici les résultat de traitement de dada2 sont organisés en un objet
phyloseq. La fonction “url” permet ici d’ouvrir l’URL dans l’objet
ps\_connect. La fonction “readRDS” permet d’écrire un seul objet R dans
un fichier et de le restaurer.

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds")
ps = readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

# Filtration taxonomique :

La fonction “rank\_names” permet de déterminer les rangs taxonomiques.

``` r
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

Ici on voit bient les rangs taxonomiques allant du kingdom au genre.

La fonction “tax\_table” permet de construire et accéder à une table
avec des noms taxonomiques sous forme de colonne. La fonction “table” va
créer un tableau à chaque combinaison. En somme, le code va créer un
tableau, avec le nombre de fonctionnalités pour chaque phyla.

``` r
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

On peut voir ici quelques phylums avec le nombre de caractéristiques qui
est observés : pour le phyla Actinobacteria il y a 13 caractéristiques
observées.

La fonction “subset\_taxa” va supprimer les OTU qui ont un résultat de
valeur manquante ou NA.

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
```

La fonction “apply” permet d’appliquer une fonction à chaque colonne
d’un tableau de données. La fonction “taxa\_are\_rows” permet
d’accéder à l’emplacement à partir des objets otu\_table. La fonction
“data.frame” va créer des collections de données. En somme, ici on
regarde les caractéristiques de l’ensemble des données.

``` r
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
```

“plyr” est un package dont dépend la fonction “ddply”. “ddply” va
appliquer, pour chaque sous-ensemble, la fonction et combiner les
résultats. Ici on va rechercher les prévalences totales et les moyennes
et les sommes des caractéristiques de chaque phyla.

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

On peut voir la moyenne de la prévalence dans la colonne 1 et la somme
de la prévalence dans la colonne 2.

Grâce à la fonction “c” on va combiner des arguments pour les filtrer
grâce à la fonction “subset\_taxa”. Ceci va retirer les entrées avec un
phyla non identifié.

``` r
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

# Filtration de prévalence :

Dans la première étape grâce à la fonction on renvoie les sous ensemble
de phyla qu’ils restent après la filtration. Puis avec la fonction
“ggplot” va entrer les données pour un graphique. Dans la seconde
étape on va inclure une estimation des paramètres grâce à la fonction
“aes”.

``` r
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->
L’ensemble des figures donne la prévalence des taxas par rapport aux
dénombrement totaux. Chaque figure répresente un taxa différent. On
peut voir que les Firmicutes ont une abondance plus importante que les
autres taxa. On peut également voir une abondance moins importante pour
Candidatus\_Saccharibacteria, Tenericutes et Verrucomicrobia.

Ici on définit le seuil de prévalence à 5% du total des échantillons.

``` r
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

La fonction “prune\_taxa” permet de retirer les taxons qui n’ont pas le
seuil de prévalence de 5%.

``` r
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

# Taxons agglomérés :

L’agglomération va permettre de limiter les redondance de nombreuses
espèces ou sous-espèces d’une communauté microbienne. On va chercher à
combiner les descendants du même genre. Ici on va chercher combien on a
de genres après la filtration. La fonction “get\_taxa\_unique” va donner
les taxons observés à un rang particulier, ici, le genre. La fonction
“length” permet de donner la longueur, qui est ici 49. On a donc 49
genres présents après la filtration.

``` r
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

La fonction “tax\_glom” permet d’agglomérer les taxons du même rang,
ici, le genre.

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
```

Ici on définit un distance de 0.4. Puis la fonction “tip\_glom” va
réunir les taxons qui sont liés et qui ont une distance inférieur à
0.4.

``` r
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

La fonction “plot\_tree” permet de tracer un arbre phylogénétique. Ici
on va créer 3 arbres phylogénétique. Dans un premier temps on va tracer
l’abre original. Puis on trace l’arbre en ayant réunit par genre et
enfin le troisième arbre sera celui de la distance de 0.4 que l’on avait
fait pour la dernière agglomération.

``` r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```

Tout d’abord il faut faire appel à la library “gridExtra” pour avoir
plusieurs tracés. Grâce à la fonction “grid.arrange” on va afficher
plusieurs tracés des trois arbres.

``` r
library("gridExtra")
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-52-1.png)<!-- -->
Globalement on peut voir que sur l’arbre regroupé en agglomération par
genre contient moins de branches que celui avant agglomération, mais
également que celui après regroupement sur une distance de 0.4. Il
semble difficile d’extraire des analyses de ces arbres car nous n’avons
pas de noms d’espèces au bout des branches.

# Transformation de la valeur d’abondance :

Ici “subset\_taxa” va permettre de filtrer un sous ensemble basé sur les
phylum. On cherche à obtenir un graphique d’abondance relative. La
fonction aes\_string permet de générer des mappages esthétiques.

``` r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

La fonction “transform\_sample\_counts” permet de transformer les
données d’abondance et ce, échantillon par échantillon dans l’objet
ps3ra.

``` r
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

La fonction “plot\_abundance” trace l’abondance relative des
échantillons. Grâce à la fonction “grid.arrange” on peut avoir les
graphiques des abondances avant transformation et après. Cela permet
d’avoir plusieurs graphiques dans un seul graphique.

``` r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->
Les 4 premiers graphiques montrent les abondances d’origines et les 4
derniers montrent les abondances après transformation des données. On
peut voir que les Clostridiales ont une abondance bien plus imporante et
ce avant et après la modificiation des données, même si on pourrait
détecter une plus grande abondance avant. On voit également pour que
les Bacillales ont une abondance moindre. On peut également dénoter deux
modules d’abondances chez les Erysipelotrichales et les Lactobacillales.
Il y a également les données en fonction du sexe de l’hôte, qui ne
semblent pas, ici, avoir un impact majeur.

# Sous-ensemble par taxonomie :

Comme nous avions remarquer une abondance bimodale chez les
Lactobacillales on va regarder graphiquement uniquement cet ordre, grâce
aux fonctions que l’on a utilisé précédemment.

``` r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-56-1.png)<!-- -->
Ici on a un graphique des abondances relatives de l’ordre des
Lactobacillales en fonction du sexe de l’hôte. Il apparait clairement
que la distibution bimodale provient d’une abondance des Lactobacillus
et des Streptococcus. Notons que les Lactobacillus ont une plus grande
abondance que les Streptococcus.

# Prétraitement :

La fonction “qplot” va permettre de tracer rapidement un histogramme,
comme on lui indique vouloir ce type de graphique, en fonction de l’âge
des sujets de l’étude.

``` r
library("ggplot2")
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-57-1.png)<!-- -->
On peut voir que trois groupes d’âges se distingues sur ce graphique.
Nous avons un premier groupe qui a un âge entre 0 et100 jours, un autre
entre 100 et 200 jours et un autre entre 300 et 400 jours, mais qui est
moins important quantitativement.

La fonction “log10” calcule le logarithme de la valeur en base 10. La
fonction “rowSums” permet de faire la somme des valeurs de ps. Cela
donne un histogramme qui va donner les profondeurs de lectures.

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->
D’après le tutoriel phyloseq cet histogramme suggère que la
transformation pourrait être suffisante pour normaliser les données
d’abondance dans les analyses exploratrices.

Ici une analyse PCoA va être faite en regroupant les âges avec une
couleur spécifique. La fonction “cut” va diviser les âges en plusieurs
sections. La fonction “list” va construire une liste en fonction des
âges. On aura donc les jeunes, les âges du milieu et les âgées. La
fonction “ordinate” va ordiner avec la méthode MDS (effectue une analyse
des coordonnées principales à l’échelle multidimensionnelle). Puis la
fonction “plot\_ordination” va tracer les résultats d’ordination.

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTGAAAGCCCCCGGCTTAACCGGGGAGGGTCATTGGAAACTGGAAGACTGGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-59-1.png)<!-- -->
Ce graphique montre l’analyse d’ordination avec les log des abondances
en fonction des groupes d’âges. On peut voir quelques valeurs
abérrantes, principalement dans le groupe d’âge du milieu. Globalement
on peut distinguer une ressemblance au niveau de l’abondance entre les
différents groupes d’âges même si les âges entre 100 et 200 jours
semblent tirer vers le haut et les jeunes vers le bas.

Ici on veut regarder, via un histogramme, les échantillons abérrants,
qui montrent une abondance relative élevée.

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->
On peut voir que les échantillons aberrants sont représentés par un seul
ASV (amplicon sequences variant).

# Différents projections d’ordination :

Ces étapes servent à calculer les ordinations après avoir supprimé les
valeurs abberrantes.

Ic, grâce à la fonction “prune\_samples” on définit un sous ensemble que
l’on souhaite conserver. La fonction “c”, précédent cette étape indique
que l’on va retirer les séquences abbérentes.

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)
```

Cette étape va permettre de déterminer quels sont les échantillosn qui
ont moins de 1000 reads.

``` r
which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

On peut voir que l’échantillon F5D145 a 69 read. L’échantillon M2D149 a
185 reads.

La fonction “prune\_samples” va ici permettre de retirer les
échantillons que l’on a déterminé précédemment, c’est à dire, ceux qui
ont moins de 1000 reads. Puis on débombre chaque échantillon
individuellement à l’aide de la fonction “transform\_sample\_counts”.

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
```

Ici on cherche à effectuer une PCoA avec la distance de Bray-Curtis.
Cela sera réalisé avec les mêmes groupes d’âges que précédemment, mais
également un regroupement en fonction des environnements.

``` r
out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->
Globalement sur la PCoA on peut voir que l’âge joue un rôle important en
fonction des échantillons. On peut voir que pour les groupes d’âges
moyens le score semble être plus importants que pour les jeunes. En
revanche il n’y a pas réellement de distinction entre les portées, elles
semblent être entremêlés. On peut également voir que les souris âgées se
situent entre les jeunes et les medium.

Ici on va créer un graphique DPCoA avec l’identifiant de l’échantillon.
La fonction “coord\_fixed” permet d’imposer un rapport précis pour la
représentation des unités sur les axes.

``` r
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->
On peut voir un allongement sur l’axe des abscisses, montrant 75% de la
diversité. La DPCoA donne des indications sur l’ordination
phylogénétique des échantillons. Nous pouvons voir à nouveau que les
âges médiants se situent à un score plus élevés que les jeunes.

On décide de réaliser une PCoA en fonction des espèces ou du phylum.

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-66-1.png)<!-- -->
On peut voir que les Bacteroidetes dénote par rapport aux autres
phhylum. En revanche nous retrouvons beaucoup de Firmicutes sur ce
graphe tandis que les autres phylum semblent moins présent.

Ici on réalise une PCoA mais avec la méthode PCoA et la distance Unifrac
pondéré.

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTCCGTAGGCGGACTTATAAGTCAGTGGTGAAAGCCTGTCGCTTAACGATAGAACTGCCATTGATACTGTAAGTCTTGAGTATATTTGAGGTAGCTGGAATAAGTAGTGTAGCGGTGAAATGCATAGATATTACTTAGAACACCAATTGCGAAGGCAGGTTACCAAGATATAACTGACGCTGAGGGACGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->
On voit une position des échantillons avec les âges moyens au dessus des
jeunes. On semble distinguer que les échantillons âgées se retouvent
entre les positions des deux autres groupes d’âges. A nouveau, il ne
semble pas y avoir de tendance en ce qui concerne les portées.

# PCA sur les rangs :

On va créer une nouvelle matrice en représentant les abondances par
leurs rangs.

``` r
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
```

Ici, tous les rangs dont l’abondance est inférieru à 329 vont être liés
à 1 pour éviter une différence de rang trop importante.

``` r
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1
```

Tout d’abord on charge les lybrary dplyr et reshape2. Puis on exécute ce
code afin de donne visualiser le rang en fonction de l’abondance.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-70-1.png)<!-- -->
On peut voir le rang en fonction de l’abondance chez différents
échantillons. On peut distinguer une corrélation entre l’abondance des
échantillons et le rang. On peut également voir que c’est dans
l’échantillon F7D125 qu’il y a une abondance plus importante et F3D65
qui en a une plus petite.

Nous allons maintenant réaliser une PCA. La fonction “dudi.pca” effectue
une analyse des données et les renvoie sous forme d’objets de classe.
Globalement, ce bloc permet d’annoter la figure.

``` r
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
col_scores <- col_scores %>%
  left_join(tax)
```

    ## Joining, by = "seq"

Puis on exécute le code suivant pour donner la PCA en fonction des
différents ordres déterminés par la légende précédente.

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->
Ici nous pouvons voir 3 graphiques de PCA, en fonction des groupes
d’âges. Il semble que l’ordre des Clostridiales soit prépondérant
chez chaque groupe d’âge. L’analyse de la PCA, avec des rangs, démontre
des similitudes avec l’analyse des PCoA. On peut voir qu’il y a une plus
grande abondance chez les jeunes et cela décline en avancant avec l’âge.

# Analyse de correspondance canonique :

Il s’agit d’une méthode qui va générer une ordination d’une espèce par
table d’échantillons en intégrant des informations supplémentaires sur
les échantillons. On peut utiliser la fonction “ordination” et on peut
affilier des échantillons supplémentaires avec d’autres
caractéristiques.

``` r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)
```

Ici on génère la CCpnA.

``` r
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

    ## Warning in class(x) <- c(setdiff(subclass, tibble_class), tibble_class): Setting
    ## class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an
    ## S4 object

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->
L’analyse de correspondance canonique est réalisé entre les portées 1
et 2. Dans les deux cas on peut distinguer une abondances de l’ordre des
Clostridiales. On ne distingue que peu des autres ordre de bactéries. Il
y aurait une légère baisse du score chez la portée 2.

# Enseignement supervisé :

Ici on va essayer de prédire l’âge à partir de la composition du
microbiome. On prend 8 souris au hasard.

``` r
library(caret)
```

    ## Loading required package: lattice

``` r
sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))
trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
plsFit <- train(age ~ ., data = training,
                method = "pls", preProc = "center")
```

Ici on va faire une prédiction de modèle grâce à la fonction “predict”.

``` r
plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$age)
```

    ##            
    ## plsClasses  (0,100] (100,400]
    ##   (0,100]        62         1
    ##   (100,400]       3        43

On peut voir ici qu’il y aurait selon ce modèle 65 souris qui auraient
entre 0 et 100, et 48 entre 100 et 400. Il y a 8 souris entre 100 et
400.

Il s’agit d’un autre exemple avec la library randomForest. Cela va
implémenter un algorithme permettant une classification.

``` r
library(randomForest)
```

    ## randomForest 4.6-14

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)
```

    ##            
    ## rfClasses   (0,100] (100,400]
    ##   (0,100]        62         1
    ##   (100,400]       3        43

On peut voir qu’avec cet algorithme les valeurs ont quelque peu
changées. On voit une augmentation des souris dans l’âge de 0 à 100 et
une diminution pour l’âge de 100 à 400.

Ici le code suivant va générer extraire les coordonnées et va permettr
d’inclure des annotations sur le biplot PLS.

``` r
library(vegan)
```

    ## Loading required package: permute

    ## This is vegan 2.5-7

    ## 
    ## Attaching package: 'vegan'

    ## The following object is masked from 'package:caret':
    ## 
    ##     tolerance

``` r
pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"

pls_biplot$scores <- data.frame(sample_data(pslog)[inTrain, ],
                                pls_biplot$scores)

tax <- tax_table(ps)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)
```

Le code suivant va permettre de générer le PLS en permettant de générer
deux graphes en fonction pour séparer les échantillons; On sépare les
échantillons pour le groupe de l’âge de 0 à 100 et un autre groupe de
100 à 400.

``` r
ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ age2) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-79-1.png)<!-- -->
Globalement on peut voir pour les deux groupes d’âges différents qu’il y
a une abondance des Clostridiales. En revanche il semble y avoir une
plus petite abondance bactérienne dans le cas des souris plus âgées.

Ici on va utiliser “randomForest” pour évaluer les proximités entre les
points des données. L’algorithme va calculer une distance entre les
échantillons. Si une paire d’échantillon se produit fréquemment, elle
aura une faible distance. On classe tout ceci dans une PCoA.

``` r
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])

ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 1, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-80-1.png)<!-- -->
On peut voir ici qu’il y a bien une séparation des classes d’âges pour
le groupe des jeunes (de 0 à 100 jours). En revanche on peut voir que
les souris plus âgées se regroupent (de 100 à 200 jours et plus de 200
jours).

Ici, on va rechercher la famille et le genre de la bactérie ayant le
plus d’influence dans la prédiction du modèle random forest.

``` r
as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
```

    ## [1] "Lachnospiraceae" "Roseburia"

Il s’agit d’une bactérie de la famille des Lachnospiracées qui a le plus
d’influence dans la prédiction.

On va étudier l’abondance de cette famille bactérienne dans
l’échantillon en fonction des jours de vie. Cette commande va générer
un histogramme permettant de voir quand il y a une abondance.

``` r
impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund)) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-82-1.png)<!-- -->
Ici on peut voir que l’abondance de la famille des Lachnospiracées est
plus importantes pour les tranches d’âges allant de 100 à 400 jours. En
effet on peut voir que l’abondance tourne autour de 2 et 3, tandis qu’on
a une abondance proche de 0 pour les souris entre 0 et 100 jours.

# Analyses basées sur des graphiques :

## Créer et tracer des graphiques :

On charge les packages nécessaires aux graphiques. “ggnetwork” trace un
réseaux en se basant sur une matrice de distance.

``` r
library("igraph")
```

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:vegan':
    ## 
    ##     diversity

    ## The following object is masked from 'package:permute':
    ## 
    ##     permute

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:ape':
    ## 
    ##     edges, mst, ring

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     union

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     union

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     normalize, path, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library("phyloseqGraphTest")
library("igraph")
library("ggnetwork")
```

La fonction “make\_network” va créer un réseau de microbiome par
échantillon basée sur une distance de 0.35. Le réseau est crée par la
matrice de dissilarité de Jaccard. On assigne un attribut pour savoir de
quelle souris il s’agit mais aussi de la portée. Il faut ajouter
net\_graph \<- ggnetwork(net) et changer net\_graph après.

``` r
net <- make_network(ps, max.dist=0.35)
sampledata <- data.frame(sample_data(ps))
V(net)$id <- sampledata[names(V(net)), "host_subject_id"]
V(net)$litter <- sampledata[names(V(net)), "family_relationship"]
net_graph <- ggnetwork(net)
```

Ici on trace le réseau en ajoutant des couleurs aux souris mais
également une forme pour la portée.

``` r
ggplot(net_graph, aes(x = x, y = y, xend = xend, yend = yend), layout = "fruchtermanreingold") +
  geom_edges(color = "darkgray") +
  geom_nodes(aes(color = id, shape = litter),  size = 3 ) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.key.height = unit(0.5,"line")) +
  guides(col = guide_legend(override.aes = list(size = .5)))
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-85-1.png)<!-- -->
Globalement ici nous pouvons déjà voir un regroupement en fonction de la
portée. En effet, la portée 2 symbolisée par des triangles présente de
nombreux regroupement, sur les extrémités du réseau. On peut voir un
regroupement pour la portée 1 également. De plus, on peut distinguer un
regroupement en fonction des souris également. On peut voir un
regroupement donc, en fonction des souris mais également de la portée.

## Tests à deux échantillons basée sur des graphiques :

## Minimum Spanning Tree (MST) :

Ici on cherche à faire un arbre couvrant le minimum qui sera basé sur
les distances entre les échantillons. On prend la distance de Jaccard.
Le but ici est de savoir si les deux portées proviennent de la même
distribution. La fonction “graph\_perm\_test” effectue des tests de
permutation basés sur des graphes.

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "mst")
gt$pval
```

    ## [1] 0.004

Ici on va tracer un graphique et un histogramme de permutation à partir
du MST. La fonction “plot\_test\_network” permet de créer un graphique
en réponse à la fonction “graph\_perm\_test”. Il y aura des noeuds
colorés par type d’échantillons, ici pour la portée 1 et pour la portée
2. On a aussi des arêtes marquées comme purs ou mixtes. Ensuite, la
fonction “plot\_permutations” trace un histogramme de la distribution de
permutation du nombre d’arêtes pures ainsi qu’une marque montrant le
nombre observé d’arêtes pures.

``` r
plotNet1=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm1=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet1, plotPerm1)
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-87-1.png)<!-- -->
On peut voir que le graphique du MST que les échantillons sont regroupés
par portée. Globalement on peut voir que les arêtes sont le plus souvent
pures. Ceci est confirmé par l’histogramme qui dénombre beaucoup
d’arêtes pures. On peut conclure que les échantillons proviennent de
deux distributions différentes.

## Voisins les plus proches :

On cherche à obtenir un graphe des voisins les plus proches en plaçant
une arête entre deux échantillons chaque fois que l’un des est dans un
ensemble d’une k proximité. Ici on définit le nombre de voisins les plus
proches à 1.

``` r
gt <- graph_perm_test(ps, "family_relationship", grouping = "host_subject_id",
                      distance = "jaccard", type = "knn", knn = 1)
```

Ici on trace à nouveau un graphique et un histogramme de permutation en
tenant compte de la proximité.

``` r
plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-89-1.png)<!-- -->
Ici on peut voir que lorsque des échantillons partagent des bords, il
sont bien souvent de la même portée. En effet, on constate beaucoup
d’échantillons avec des arêtes communes qui sont de la même portée.
Encore une fois, grâce à l’histogramme on retrouve beaucoup d’arêtes
pures.

## Modélisation linéaire :

Le but ici est de décrire comment un environnement va agir sur la
structure globale de la communauté. La fonction “estimate\_richness”
permet de résumer la diversité alpha, et renvoie les résultats sous
forme de data.frame. On utilise pour ça la mesure de diversité de
Shannon.

``` r
library("nlme")
```

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     collapse

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     collapse

``` r
library("reshape2")
ps_alpha_div <- estimate_richness(ps, split = TRUE, measure = "Shannon")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>%
  as.factor()
ps_samp <- sample_data(ps) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ps_alpha_div, by = "SampleID") %>%
  melt(measure.vars = "Shannon",
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")

# reorder's facet from lowest to highest diversity
diversity_means <- ps_samp %>%
  group_by(host_subject_id) %>%
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
ps_samp$host_subject_id <- factor(ps_samp$host_subject_id)
#                                  diversity_means$host_subject_id)
```

``` r
alpha_div_model <- lme(fixed = alpha_diversity ~ age_binned, data = ps_samp,
                       random = ~ 1 | host_subject_id)
```

``` r
new_data <- expand.grid(host_subject_id = levels(ps_samp$host_subject_id),
                        age_binned = levels(ps_samp$age_binned))
new_data$pred <- predict(alpha_div_model, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model$sigma ^ 2
```

``` r
# fitted values, with error bars
ggplot(ps_samp %>% left_join(new_data)) +
  geom_errorbar(aes(x = age_binned, ymin = pred - 2 * sqrt(pred_var),
                    ymax = pred + 2 * sqrt(pred_var)),
                col = "#858585", size = .1) +
  geom_point(aes(x = age_binned, y = alpha_diversity,
                 col = family_relationship), size = 0.8) +
  facet_wrap(~host_subject_id) +
  scale_y_continuous(limits = c(2.4, 4.6), breaks = seq(0, 5, .5)) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Binned Age", y = "Shannon Diversity", color = "Litter") +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)),
        axis.text.x = element_text(angle = -90, size = 6),
        axis.text.y = element_text(size = 6))
```

    ## Joining, by = c("host_subject_id", "age_binned")

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-93-1.png)<!-- -->
Ici on peut voir qu’il y a une diversité plus importante chez toutes les
souris à l’âge de 0 à 100 jours. La diversité diminue et on ne distingue
qu’une très faible diversité de 200 à 400 jours. On distingue un profil
comparable peu importe la portée. La diversité décline avec l’âge. Il ne
semble pas y avoir de tendance en fonction du sexe non plus.

## Tests multiples hierarchiques :

Ici on va chercher à savoir si des bactéries abondantes peuvent être
liées à l’âge. Pour se faire on utilise un test hierarchique, c’est à
dire que nous n’allons tester des groupes taxonomiques que si les
niveaux plus élevés sont associés. On va générer un histogramme avec
DESeq.

``` r
library("reshape2")
library("DESeq2")
```

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     sampleNames

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:igraph':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
#New version of DESeq2 needs special levels
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship = gsub(" ", "", sample_data(ps)$family_relationship)
ps_dds <- phyloseq_to_deseq2(ps, design = ~ age_binned + family_relationship)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
# geometric mean, set to zero when all coordinates are zero
geo_mean_protected <- function(x) {
  if (all(x == 0)) {
    return (0)
  }
  exp(mean(log(x[x != 0])))
}

geoMeans <- apply(counts(ps_dds), 1, geo_mean_protected)
ps_dds <- estimateSizeFactors(ps_dds, geoMeans = geoMeans)
ps_dds <- estimateDispersions(ps_dds)
```

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

``` r
abund <- getVarianceStabilizedData(ps_dds)
```

``` r
short_names <- substr(rownames(abund), 1, 5)%>%
  make.names(unique = TRUE)
rownames(abund) <- short_names
```

``` r
abund_sums <- rbind(data.frame(sum = colSums(abund),
                               sample = colnames(abund),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(pslog)),
                               sample = rownames(otu_table(pslog)),
                               type = "log(1 + x)"))

ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 20) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-96-1.png)<!-- -->
La figure ci dessous représente l’abondance de la transformation DEseq.
Le premier histogramme donne l’abondance totale de DESeq2 dans chaque
échantillon. Le second histogramme a été repris pour faire une
comparaison. On peut voir qu’il y a une plus grande abondance sur
l’histogramme DESeq2, même si elles semblent plus étalées sur
l’histogramme du bas.

Ces tests nécessitent de faire des tests univariés pour chaque groupe
taxnonomique, et cela se fait grâce à la fonction “treePValues”.

``` r
library("structSSI")
el <- phy_tree(pslog)$edge
el0 <- el
el0 <- el0[nrow(el):1, ]
el_names <- c(short_names, seq_len(phy_tree(pslog)$Nnode))
el[, 1] <- el_names[el0[, 1]]
el[, 2] <- el_names[as.numeric(el0[, 2])]
unadj_p <- treePValues(el, abund, sample_data(pslog)$age_binned)
```

``` r
hfdr_res <- hFDR.adjust(unadj_p, el, .75)
summary(hfdr_res)
```

    ## Number of hypotheses: 764 
    ## Number of tree discoveries: 579 
    ## Estimated tree FDR: 1 
    ## Number of tip discoveries: 280 
    ## Estimated tips FDR: 1 
    ## 
    ##  hFDR adjusted p-values: 
    ##                 unadjp         adjp adj.significance
    ## GCAAG.95  1.861873e-82 3.723745e-82              ***
    ## GCAAG.70  1.131975e-75 2.263950e-75              ***
    ## GCAAG.187 5.148758e-59 1.029752e-58              ***
    ## GCAAG.251 3.519276e-50 7.038553e-50              ***
    ## GCAAG.148 1.274481e-49 2.548962e-49              ***
    ## GCAAG.30  9.925218e-49 1.985044e-48              ***
    ## GCGAG.76  1.722591e-46 3.445183e-46              ***
    ## GCAAG.167 6.249050e-43 1.249810e-42              ***
    ## 255       8.785479e-40 1.757096e-39              ***
    ## GCAAG.64  2.727610e-36 5.455219e-36              ***
    ## [only 10 most significant hypotheses shown] 
    ## --- 
    ## Signif. codes:  0 '***' 0.015 '**' 0.15 '*' 0.75 '.' 1.5 '-' 1

``` r
#interactive part: not run
plot(hfdr_res, height = 5000) # opens in a browser
```

Ici nous devrions avoir une image qui s’ouvre dans une fenêtre internet.
La commande ne fonctionne pas. Cela renvoie un sous arbre avec des
bactéries qui ont différentes abondances. Il en ressort qu’une
association entre le groupe d’âge et l’abondance bactérienne n’est pas
une tendance majeure. Cela apparait pour quelques groupes taxonomiques.

Ici, on va regarder l’identité taxonomique des bactéries qui donnent une
corrélation entre abondance bactérienne et âge.

``` r
tax <- tax_table(pslog)[, c("Family", "Genus")] %>%
  data.frame()
tax$seq <- short_names
```

``` r
options(digits=3)
hfdr_res@p.vals$seq <- rownames(hfdr_res@p.vals)
tax %>%
  left_join(hfdr_res@p.vals) %>%
  arrange(adjp) %>% head(10)
```

    ## Joining, by = "seq"

    ##             Family            Genus       seq   unadjp     adjp
    ## 1  Lachnospiraceae             <NA>  GCAAG.95 1.86e-82 3.72e-82
    ## 2  Lachnospiraceae        Roseburia  GCAAG.70 1.13e-75 2.26e-75
    ## 3  Lachnospiraceae Clostridium_XlVa GCAAG.187 5.15e-59 1.03e-58
    ## 4  Lachnospiraceae             <NA> GCAAG.251 3.52e-50 7.04e-50
    ## 5  Lachnospiraceae Clostridium_XlVa GCAAG.148 1.27e-49 2.55e-49
    ## 6  Lachnospiraceae             <NA>  GCAAG.30 9.93e-49 1.99e-48
    ## 7  Ruminococcaceae     Ruminococcus  GCGAG.76 1.72e-46 3.45e-46
    ## 8  Lachnospiraceae Clostridium_XlVa GCAAG.167 6.25e-43 1.25e-42
    ## 9  Lachnospiraceae        Roseburia  GCAAG.64 2.73e-36 5.46e-36
    ## 10            <NA>             <NA>   GCAAG.1 5.22e-35 1.04e-34
    ##    adj.significance
    ## 1               ***
    ## 2               ***
    ## 3               ***
    ## 4               ***
    ## 5               ***
    ## 6               ***
    ## 7               ***
    ## 8               ***
    ## 9               ***
    ## 10              ***

On peut voir que la majorité des bactéries qui montrent une corrélation
entre l’abondance et le groupe d’âge sont des bactéries de la famille
des Lachnospiraceae. Rappelons nous, que, plus haut nous avions déjà
retrouvé cette famille, en effet elle avait déjà le plus d’influence
dans le modèle avec prédiction.

# Multitable techniques :

Ici une corrélation canonique clairsemée est réalisée CCA. Il s agit
d’une méthode pour comparer les échantillons et identifier des
caractéristiques qui présentent des variatiosn intéressantes. Ici nous
utilisons un nouvel ensembke de données : avec deux tableaux, un pour
les bactéries et un avec les métabolites.

``` r
metab <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/metabolites.csv",row.names = 1)
microbe_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/microbe.rda")
load(microbe_connect)
microbe
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 20609 taxa and 12 samples ]
    ## tax_table()   Taxonomy Table:    [ 20609 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 20609 tips and 20607 internal nodes ]

Ici on a filtré les bactéries et les métabolites que l’on souhaite
étudiés. On va supprimer les métabolites que l’on ne retrouve pas dans
nos échantillons.

``` r
library("genefilter")
```

    ## 
    ## Attaching package: 'genefilter'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     rowSds, rowVars

``` r
keep_ix <- rowSums(metab == 0) <= 3
metab <- metab[keep_ix, ]
microbe <- prune_taxa(taxa_sums(microbe) > 4, microbe)
microbe <- filter_taxa(microbe, filterfun(kOverA(3, 2)), TRUE)
metab <- log(1 + metab, base = 10)
X <- otu_table(microbe)
X[X > 50] <- 50
dim(X)
```

    ## [1] 174  12

``` r
dim(metab)
```

    ## [1] 405  12

On réalise ici une CCA clairsemée. Cela permet de comparer des ensembles
de données dans des tables qui ont des dimensions très grandes. La
technique va sélectionner des sous-ensemble qui vont refléter le plus de
covariance. Ceci signifie que les signaux sont un reflet correct.

``` r
library(PMA)
cca_res <- CCA(t(X),  t(metab), penaltyx = .15, penaltyz = .15)
```

    ## 123456789101112131415

``` r
cca_res
```

    ## Call: CCA(x = t(X), z = t(metab), penaltyx = 0.15, penaltyz = 0.15)
    ## 
    ## 
    ## Num non-zeros u's:  5 
    ## Num non-zeros v's:  15 
    ## Type of x:  standard 
    ## Type of z:  standard 
    ## Penalty for x: L1 bound is  0.15 
    ## Penalty for z: L1 bound is  0.15 
    ## Cor(Xu,Zv):  0.974

Avec les paramètres que l’on vient de définir, on a sélectionnés 6
bactéries et 11 métabolites. Cela engendre une corrélation à 0.976
entre les deux tableaux ce qui est significatif. On peut donc dire que
les données des bactéries et des métabolites reflètent des signaux
similaires.

Ici, une PCA est réalisée pour relier les métabolites et les bactéries
aux échantillons.

``` r
combined <- cbind(t(X[cca_res$u != 0, ]),
                  t(metab[cca_res$v != 0, ]))
pca_res <- dudi.pca(combined, scannf = F, nf = 3)
```

``` r
genotype <- substr(rownames(pca_res$li), 1, 2)
sample_type <- substr(rownames(pca_res$l1), 3, 4)
feature_type <- grepl("\\.", colnames(combined))
feature_type <- ifelse(feature_type, "Metabolite", "OTU")
sample_info <- data.frame(pca_res$li, genotype, sample_type)
feature_info <- data.frame(pca_res$c1,
                           feature = substr(colnames(combined), 1, 6))
```

``` r
ggplot() +  geom_point(data = sample_info,
            aes(x = Axis1, y = Axis2, col = sample_type, shape = genotype), size = 3) + 
  geom_label_repel(data = feature_info,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = feature_type),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = feature_info,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = feature_type),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed(sqrt(pca_res$eig[2] / pca_res$eig[2])) +
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pca_res$eig[1] / sum(pca_res$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pca_res$eig[2] / sum(pca_res$eig), 2)),
       fill = "Feature Type", col = "Sample Type")
```

![](Tutoriel-phyloseq_files/figure-gfm/unnamed-chunk-109-1.png)<!-- -->
Ici nous pouvons voir les variations entre les échantillons. Ce graphe
donne un lien entre le génotype : type sauvage ou knockout. PD et ST
correspondent aux différents régimes. Les plus grandes variations sont
entre les régimes alimentaires PD et ST.

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or
by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output
will be saved alongside it (click the *Preview* button or press
*Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the
editor. Consequently, unlike *Knit*, *Preview* does not run any R code
chunks. Instead, the output of the chunk when it was last run in the
editor is displayed.
