Evaluating feature tables with R
================
Bernice Waweru
Tue 22, Jun 2021

-   [Analysis in R](#analysis-in-r)
    -   [Install the packages and prepare the input
        data.](#install-the-packages-and-prepare-the-input-data)
    -   [Use the functions](#use-the-functions)
        -   [Function `*otuReport*` to summarise the two
            datasets](#function-otureport-to-summarise-the-two-datasets)
        -   [Calculate Alpha diversity](#calculate-alpha-diversity)
    -   [Working with the merged
        dataset](#working-with-the-merged-dataset)

## Analysis in R

[**otuSummary**](https://github.com/cam315/otuSummary/) is an r package
that provide functions to summarize various aspects of on otu table
generated from either *qiime* or *mothur*. We use this package to
compare the fwd and reverse tables generated for the data from qiime2.

### Install the packages and prepare the input data.

``` r
#check example data from qiime

data("otuqiime")
dim(otuqiime)
names(otuqiime)[1:10]
rownames(otuqiime)[1:10]
head(otuqiime$taxonomy) #the OTU and taxonomy table are combined into one dataframe, taxa as rows and the samples as columns
str(otuqiime)

#load the data from qiime from forward and reverse and combine the feature and taxonomy tables by the feature ID
#this will generate the one data frame that we need to use with the package

#============== forward ================

fwd_otu <- read.table("C:/Users/BWaweru/OneDrive - CGIAR/Documents/Fellows/Yves_Armel/jan_2021/exported_files/yves_fwd_deblur_marjm_265_table.txt", header = T)

# reading the taxonomy table with read.table gives and error that some lines are missing elements.
# using read.csv however works

fwd_tax <- read.csv("C:/Users/BWaweru/OneDrive - CGIAR/Documents/Fellows/Yves_Armel/jan_2021/exported_files/yves_fwd_deblur_marjm_265_taxonomy.txt", sep = '\t', header = T)

colnames(fwd_tax) <- c("OTU_ID", "taxonomy")

#========== merge 

merge(fwd_otu, fwd_tax, by = "OTU_ID") -> fwd

#the OTU_ID column needs to be the rownames

rownames(fwd) <- fwd$OTU_ID
rownames(fwd)[1:10]
fwd$OTU_ID <- NULL
names(fwd)
#============= reverse ==============

rev_otu <- read.table("C:/Users/BWaweru/OneDrive - CGIAR/Documents/Fellows/Yves_Armel/jan_2021/exported_files/yves_rev_deblur_marjm_240_table.txt", header = T)

rev_tax <- read.csv("C:/Users/BWaweru/OneDrive - CGIAR/Documents/Fellows/Yves_Armel/jan_2021/exported_files/yves_rev_deblur_marjm_240_repseqs_taxonomy.tsv", header = T, sep = '\t')

#a bit of cleaning before merging
rev_tax$Confidence <- NULL
colnames(rev_tax) <- c("OTU_ID", "taxonomy")
colnames(rev_otu) <- c("OTU_ID","S10","S12","S2","S3","S4","S6","S7")
names(rev_otu)

#merge
merge(rev_otu, rev_tax, by="OTU_ID") -> rev

#assgn the rownames

rownames(rev) <- rev$OTU_ID
rownames(rev)[1:10]
rev$OTU_ID <- NULL
names(rev)
```

### Use the functions

#### Function `*otuReport*` to summarise the two datasets

``` r
otuReport(fwd, siteInCol = T, taxhead = "taxonomy", platform = "qiime", pattern = ";", percent = F, taxlevel = "family") -> fwd_otu_report

otuReport(rev, siteInCol = T, taxhead = "taxonomy", platform = "qiime", pattern = ";", percent = F, taxlevel = "family") -> rev_otu_report

#Showing the reads from each family present per sample

kableExtra::kable(as.data.frame(fwd_otu_report[[4]]), caption = "Families present in FORWARD reads with reads per sample")
```

<table>
<caption>
Families present in FORWARD reads with reads per sample
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
S10
</th>
<th style="text-align:right;">
S12
</th>
<th style="text-align:right;">
S2
</th>
<th style="text-align:right;">
S3
</th>
<th style="text-align:right;">
S4
</th>
<th style="text-align:right;">
S6
</th>
<th style="text-align:right;">
S7
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
1316
</td>
<td style="text-align:right;">
192
</td>
<td style="text-align:right;">
691
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
386
</td>
<td style="text-align:right;">
11330
</td>
<td style="text-align:right;">
2682
</td>
<td style="text-align:right;">
2775
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
185
</td>
<td style="text-align:right;">
141
</td>
<td style="text-align:right;">
142
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
449
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1331
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
425
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
33813
</td>
<td style="text-align:right;">
35344
</td>
<td style="text-align:right;">
2897
</td>
<td style="text-align:right;">
2734
</td>
<td style="text-align:right;">
18511
</td>
<td style="text-align:right;">
402
</td>
<td style="text-align:right;">
1835
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
125050
</td>
<td style="text-align:right;">
129188
</td>
<td style="text-align:right;">
83106
</td>
<td style="text-align:right;">
104039
</td>
<td style="text-align:right;">
198943
</td>
<td style="text-align:right;">
26898
</td>
<td style="text-align:right;">
52395
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
441
</td>
<td style="text-align:right;">
183
</td>
<td style="text-align:right;">
370
</td>
<td style="text-align:right;">
168
</td>
<td style="text-align:right;">
1168
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
402
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
2650
</td>
<td style="text-align:right;">
1271
</td>
<td style="text-align:right;">
802
</td>
<td style="text-align:right;">
1641
</td>
<td style="text-align:right;">
6416
</td>
<td style="text-align:right;">
344
</td>
<td style="text-align:right;">
1840
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
356
</td>
<td style="text-align:right;">
220
</td>
<td style="text-align:right;">
53
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
273
</td>
<td style="text-align:right;">
289
</td>
<td style="text-align:right;">
2380
</td>
<td style="text-align:right;">
6088
</td>
<td style="text-align:right;">
1708
</td>
<td style="text-align:right;">
1303
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
kableExtra::kable(as.data.frame(rev_otu_report[[4]]), caption = "Families present in REVERSE reads with reads per sample")
```

<table>
<caption>
Families present in REVERSE reads with reads per sample
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
S10
</th>
<th style="text-align:right;">
S12
</th>
<th style="text-align:right;">
S2
</th>
<th style="text-align:right;">
S3
</th>
<th style="text-align:right;">
S4
</th>
<th style="text-align:right;">
S6
</th>
<th style="text-align:right;">
S7
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
581
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
617
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
141
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
386
</td>
<td style="text-align:right;">
15252
</td>
<td style="text-align:right;">
3886
</td>
<td style="text-align:right;">
2679
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
368
</td>
<td style="text-align:right;">
139
</td>
<td style="text-align:right;">
184
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
2130
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1273
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
190
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
41050
</td>
<td style="text-align:right;">
42453
</td>
<td style="text-align:right;">
3303
</td>
<td style="text-align:right;">
2362
</td>
<td style="text-align:right;">
22967
</td>
<td style="text-align:right;">
232
</td>
<td style="text-align:right;">
1678
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
173693
</td>
<td style="text-align:right;">
176749
</td>
<td style="text-align:right;">
92366
</td>
<td style="text-align:right;">
127660
</td>
<td style="text-align:right;">
271817
</td>
<td style="text-align:right;">
30382
</td>
<td style="text-align:right;">
52585
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
5182
</td>
<td style="text-align:right;">
3273
</td>
<td style="text-align:right;">
849
</td>
<td style="text-align:right;">
267
</td>
<td style="text-align:right;">
3242
</td>
<td style="text-align:right;">
209
</td>
<td style="text-align:right;">
443
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
2734
</td>
<td style="text-align:right;">
1704
</td>
<td style="text-align:right;">
1350
</td>
<td style="text-align:right;">
2235
</td>
<td style="text-align:right;">
12304
</td>
<td style="text-align:right;">
703
</td>
<td style="text-align:right;">
2919
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
536
</td>
<td style="text-align:right;">
1134
</td>
<td style="text-align:right;">
165
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
287
</td>
<td style="text-align:right;">
198
</td>
<td style="text-align:right;">
1576
</td>
<td style="text-align:right;">
4994
</td>
<td style="text-align:right;">
1549
</td>
<td style="text-align:right;">
815
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
#Checking what taxa are present in both reveals 12 families in both the forward and reverse, very similar

kableExtra::kable(fwd_otu_report[[1]], caption = "Taxa down to family level present in FORWARD reads") #12 families
```

<table>
<caption>
Taxa down to family level present in FORWARD reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
</tr>
</tbody>
</table>

``` r
kable(rev_otu_report[[1]], caption = "Taxa down to family level present in REVERSE reads") #12 families 
```

<table>
<caption>
Taxa down to family level present in REVERSE reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
</tr>
</tbody>
</table>

``` r
#Count abundance of the present taxa

kable(as.data.frame(fwd_otu_report[[2]] %>% sort(decreasing = T)), caption = "Count abundance of taxa down to family present FORWARD reads")
```

<table>
<caption>
Count abundance of taxa down to family present FORWARD reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
Var1
</th>
<th style="text-align:right;">
Freq
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
8961
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
464
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
276
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
228
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
220
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
78
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
58
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
42
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
23
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>

``` r
kable(as.data.frame(rev_otu_report[[2]] %>% sort(decreasing = T)), caption = "Count abundance of taxa down to family present in REVERSE reads")
```

<table>
<caption>
Count abundance of taxa down to family present in REVERSE reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
Var1
</th>
<th style="text-align:right;">
Freq
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
10868
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
537
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
397
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
301
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
207
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
101
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
54
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
48
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
40
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
3
</td>
</tr>
</tbody>
</table>

``` r
#Readsum of taxa i.e reads observed for each taxa at family level
kableExtra::kable(fwd_otu_report[[5]] %>% sort(decreasing = T), caption = "Reads observed per taxa at family level in the FORWARD reads")
```

<table>
<caption>
Reads observed per taxa at family level in the FORWARD reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
719619
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
95536
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
17226
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
14964
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
12041
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
2763
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
2241
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
1756
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
708
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
495
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
470
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
73
</td>
</tr>
</tbody>
</table>

``` r
kable(rev_otu_report[[5]] %>% sort(decreasing = T), caption = "Reads observed per taxa at family level in the REVERSE reads")
```

<table>
<caption>
Reads observed per taxa at family level in the REVERSE reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
925252
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
114045
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
23949
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
22374
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
13465
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
9419
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
2235
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
2089
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
1493
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
1463
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
691
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
50
</td>
</tr>
</tbody>
</table>

``` r
#Relative abundance of present taxa per sample
kable(fwd_otu_report[[8]] %>% as.data.frame(), caption = "Relative abundance of observed taxa in FORWARD reads per sample" )
```

<table>
<caption>
Relative abundance of observed taxa in FORWARD reads per sample
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
S10
</th>
<th style="text-align:right;">
S12
</th>
<th style="text-align:right;">
S2
</th>
<th style="text-align:right;">
S3
</th>
<th style="text-align:right;">
S4
</th>
<th style="text-align:right;">
S6
</th>
<th style="text-align:right;">
S7
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
0.0012322176
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0524816893
</td>
<td style="text-align:right;">
0.0013052729
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
0.0055449790
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0343174698
</td>
<td style="text-align:right;">
1.0156750457
</td>
<td style="text-align:right;">
0.0835374635
</td>
<td style="text-align:right;">
2.0883072925
</td>
<td style="text-align:right;">
0.0035414527
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
0.0289571127
</td>
<td style="text-align:right;">
0.0036083497
</td>
<td style="text-align:right;">
0.4273078498
</td>
<td style="text-align:right;">
8.7443755837
</td>
<td style="text-align:right;">
1.1669139434
</td>
<td style="text-align:right;">
8.3864728460
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0022140303
</td>
<td style="text-align:right;">
0.1427810665
</td>
<td style="text-align:right;">
0.0613478248
</td>
<td style="text-align:right;">
0.4291456375
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0033210455
</td>
<td style="text-align:right;">
0.3465335072
</td>
<td style="text-align:right;">
0.0087018191
</td>
<td style="text-align:right;">
0.0695095047
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
1.0272518890
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
1.2844147602
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
20.8324861838
</td>
<td style="text-align:right;">
21.2555854247
</td>
<td style="text-align:right;">
3.2070229041
</td>
<td style="text-align:right;">
2.1100726254
</td>
<td style="text-align:right;">
8.0539686822
</td>
<td style="text-align:right;">
1.2149052555
</td>
<td style="text-align:right;">
3.2492828558
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
77.0444029598
</td>
<td style="text-align:right;">
77.6925806316
</td>
<td style="text-align:right;">
91.9996014745
</td>
<td style="text-align:right;">
80.2962128287
</td>
<td style="text-align:right;">
86.5583000126
</td>
<td style="text-align:right;">
81.2898546345
</td>
<td style="text-align:right;">
92.7772072104
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
0.2717039721
</td>
<td style="text-align:right;">
0.1100546665
</td>
<td style="text-align:right;">
0.4095956074
</td>
<td style="text-align:right;">
0.1296606441
</td>
<td style="text-align:right;">
0.5081862363
</td>
<td style="text-align:right;">
0.0936867237
</td>
<td style="text-align:right;">
0.7118319935
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
1.6326882674
</td>
<td style="text-align:right;">
0.7643687493
</td>
<td style="text-align:right;">
0.8878261543
</td>
<td style="text-align:right;">
1.2665066490
</td>
<td style="text-align:right;">
2.7915435722
</td>
<td style="text-align:right;">
1.0396204177
</td>
<td style="text-align:right;">
3.2581364876
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
0.0147866107
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.3940973952
</td>
<td style="text-align:right;">
0.1697937007
</td>
<td style="text-align:right;">
0.0230598207
</td>
<td style="text-align:right;">
0.1662183807
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
0.1681976970
</td>
<td style="text-align:right;">
0.1738021782
</td>
<td style="text-align:right;">
2.6346960690
</td>
<td style="text-align:right;">
4.6986547708
</td>
<td style="text-align:right;">
0.7431353524
</td>
<td style="text-align:right;">
3.9378645471
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
</tbody>
</table>

``` r
kable(rev_otu_report[[8]], caption = "Relative abundance of observed taxa in REVERSE reads per sample")
```

<table>
<caption>
Relative abundance of observed taxa in REVERSE reads per sample
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
S10
</th>
<th style="text-align:right;">
S12
</th>
<th style="text-align:right;">
S2
</th>
<th style="text-align:right;">
S3
</th>
<th style="text-align:right;">
S4
</th>
<th style="text-align:right;">
S6
</th>
<th style="text-align:right;">
S7
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0315843998
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
0.0094111320
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0438181547
</td>
<td style="text-align:right;">
0.3670107261
</td>
<td style="text-align:right;">
0.0727118682
</td>
<td style="text-align:right;">
1.6992096059
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
0.0631890293
</td>
<td style="text-align:right;">
0.0133683286
</td>
<td style="text-align:right;">
0.3844047204
</td>
<td style="text-align:right;">
9.6345053251
</td>
<td style="text-align:right;">
1.2285144333
</td>
<td style="text-align:right;">
7.3779295530
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.2324611828
</td>
<td style="text-align:right;">
0.0439432595
</td>
<td style="text-align:right;">
0.5067334967
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
0.0053777897
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0049793358
</td>
<td style="text-align:right;">
1.3454954329
</td>
<td style="text-align:right;">
0.0056904940
</td>
<td style="text-align:right;">
0.1927790477
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.8041388198
</td>
<td style="text-align:right;">
0.0000000000
</td>
<td style="text-align:right;">
0.5232574151
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
18.3965223626
</td>
<td style="text-align:right;">
18.9175218684
</td>
<td style="text-align:right;">
3.2893492008
</td>
<td style="text-align:right;">
1.4920470481
</td>
<td style="text-align:right;">
7.2607542434
</td>
<td style="text-align:right;">
0.6389248437
</td>
<td style="text-align:right;">
2.9119305857
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
77.8403692749
</td>
<td style="text-align:right;">
78.7612906676
</td>
<td style="text-align:right;">
91.9842652990
</td>
<td style="text-align:right;">
80.6412896542
</td>
<td style="text-align:right;">
85.9318342043
</td>
<td style="text-align:right;">
83.6716146622
</td>
<td style="text-align:right;">
91.2537960954
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
2.3223088644
</td>
<td style="text-align:right;">
1.4584846554
</td>
<td style="text-align:right;">
0.8454912115
</td>
<td style="text-align:right;">
0.1686606951
</td>
<td style="text-align:right;">
1.0249212025
</td>
<td style="text-align:right;">
0.5755831566
</td>
<td style="text-align:right;">
0.7687635575
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
1.2252397598
</td>
<td style="text-align:right;">
0.7593210671
</td>
<td style="text-align:right;">
1.3444206543
</td>
<td style="text-align:right;">
1.4118226725
</td>
<td style="text-align:right;">
3.8897688079
</td>
<td style="text-align:right;">
1.9360524359
</td>
<td style="text-align:right;">
5.0655097614
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
0.0089629829
</td>
<td style="text-align:right;">
0.0017824438
</td>
<td style="text-align:right;">
0.5337847931
</td>
<td style="text-align:right;">
0.7163341882
</td>
<td style="text-align:right;">
0.0521628619
</td>
<td style="text-align:right;">
0.6334168709
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
0.1286188043
</td>
<td style="text-align:right;">
0.0882309691
</td>
<td style="text-align:right;">
1.5694866305
</td>
<td style="text-align:right;">
3.1546498553
</td>
<td style="text-align:right;">
0.4896986251
</td>
<td style="text-align:right;">
2.2444989122
</td>
<td style="text-align:right;">
0.0000000000
</td>
</tr>
</tbody>
</table>

``` r
#save as a csv the relative abudance table

write.csv(fwd_otu_report[[8]] %>% as.data.frame(), file = "fwd_relative_abund_df.csv", quote = F)
```

> At family level above we do not see a difference in the taxa
> identified. Only with the abundance, with the most abundant families
> i.e **Glomeraceae**, **Claroideoglomeraceae**, and
> **Diversisporaceae** in both the forward and reverse reads being more
> abundant within the revere reads. This is observed consistently in
> terms of *`reads per sample`*, *`count abundance`*, *`read sums`* and
> *`relative abundance per sample`*.

Can we see a difference at genus level?

Let’s find out.

``` r
otuReport(fwd, siteInCol = T, taxhead = "taxonomy", platform = "qiime", pattern = ";", percent = F, taxlevel = "genus") -> fwd_genus_otu_report

otuReport(rev, siteInCol = T, taxhead = "taxonomy", platform = "qiime", pattern = ";", percent = F, taxlevel = "genus") -> rev_genus_otu_report
```

First we look at the how many genus we have from each datas set

``` r
nrow(as.data.frame(fwd_genus_otu_report[[1]])) #129 genus
```

    ## [1] 129

``` r
nrow(as.data.frame(rev_genus_otu_report[[1]])) #136 genus
```

    ## [1] 136

> At **genus** level is where we see a difference between the two
> datasets, at family level they were the same, both had 12 families.

> Which are this genera that are different between the two?

We manipulate a little the result of above so we can compare the
sixth(genus)column of the forward and reverse dataframes (dfs) to find
the differences

We split based on *-&gt;*, then put the results into a dataframe

``` r
fwd_df <- data.frame(do.call(rbind, strsplit(as.vector(fwd_genus_otu_report[[1]]), "->", fixed=TRUE)))
rev_df <- data.frame(do.call(rbind, strsplit(as.vector(rev_genus_otu_report[[1]]), "->", fixed=TRUE)))
```

Now we do an `anti_join`, of the two.This function returns rows from the
first df that are not common between the two based on the column of
comparison, i.e unique to fwd\_df.

``` r
kable(as.data.frame(anti_join(fwd_df, rev_df, by="X6")), caption = "Genera present only in taxa identified from the FORWARD reads")
```

<table>
<caption>
Genera present only in taxa identified from the FORWARD reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
X1
</th>
<th style="text-align:left;">
X2
</th>
<th style="text-align:left;">
X3
</th>
<th style="text-align:left;">
X4
</th>
<th style="text-align:left;">
X5
</th>
<th style="text-align:left;">
X6
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Toljander08-Glo8\_VTX00114
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Diversisporales
</td>
<td style="text-align:left;">
Acaulosporaceae
</td>
<td style="text-align:left;">
Acaulospora\_sp.\_VTX00014
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00304
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Franke A\_VTX00204
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Paraglomeromycetes
</td>
<td style="text-align:left;">
Paraglomerales
</td>
<td style="text-align:left;">
Paraglomeraceae
</td>
<td style="text-align:left;">
Paraglomus\_Para2-OTU13\_VTX00281
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Diversisporales
</td>
<td style="text-align:left;">
Diversisporaceae
</td>
<td style="text-align:left;">
Diversispora\_Torrecillas12b Div1\_VTX00054
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Alguacil11d Glo G17\_VTX00172
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00156
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Diversisporales
</td>
<td style="text-align:left;">
Gigasporaceae
</td>
<td style="text-align:left;">
Gigaspora\_sp.\_VTX00039
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00125
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_NF14\_VTX00167
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00140
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Glo-F\_VTX00167
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_NF08\_VTX00154
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Diversisporales
</td>
<td style="text-align:left;">
Pacisporaceae
</td>
<td style="text-align:left;">
Pacispora\_Schechter08 Paci1\_VTX00011
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Alguacil11d Glo G9\_VTX00172
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Alguacil11d Glo G1\_VTX00113
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Franke A3\_VTX00092
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Archaeosporomycetes
</td>
<td style="text-align:left;">
Archaeosporales
</td>
<td style="text-align:left;">
Archaeosporaceae
</td>
<td style="text-align:left;">
Archaeospora\_sp.\_VTX00005
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Toljander08-Glo4\_VTX00143
</td>
</tr>
</tbody>
</table>

> We have **20** genera that are present in the fwd dataset not present
> in the reverse.

Let’s do the other way,

``` r
kable(as.data.frame(anti_join(rev_df, fwd_df, by="X6")), caption= "Genera present only in taxa identified from the REVERSE reads" )
```

<table>
<caption>
Genera present only in taxa identified from the REVERSE reads
</caption>
<thead>
<tr>
<th style="text-align:left;">
X1
</th>
<th style="text-align:left;">
X2
</th>
<th style="text-align:left;">
X3
</th>
<th style="text-align:left;">
X4
</th>
<th style="text-align:left;">
X5
</th>
<th style="text-align:left;">
X6
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Alguacil11d Glo G12\_VTX00064
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Diversisporales
</td>
<td style="text-align:left;">
Acaulosporaceae
</td>
<td style="text-align:left;">
Acaulospora\_MO-A5\_VTX00026
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_geosporum
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Torrecillas12b Glo G13\_VTX00409
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Claroideoglomeraceae
</td>
<td style="text-align:left;">
Claroideoglomus\_claroideum
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00063
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00414
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_NF05\_VTX00322
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Paraglomeromycetes
</td>
<td style="text-align:left;">
Paraglomerales
</td>
<td style="text-align:left;">
Paraglomeraceae
</td>
<td style="text-align:left;">
Paraglomus\_Alguacil12b PARA1\_VTX00350
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Claroideoglomeraceae
</td>
<td style="text-align:left;">
Claroideoglomus\_GlBb1.3\_VTX00055
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Glo72\_VTX00344
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Diversisporales
</td>
<td style="text-align:left;">
Diversisporaceae
</td>
<td style="text-align:left;">
Diversispora\_MO-GC1\_VTX00060
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_iranicum\_VTX00155
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00212
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00195
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00085
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_NF13\_VTX00419
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_PODO16.1
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_sp.\_VTX00053
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Glo G4\_VTX00419
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_NES19\_VTX00086
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Claroideoglomeraceae
</td>
<td style="text-align:left;">
Claroideoglomus\_Glo G8\_VTX00057
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_JP3\_VTX00128
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Glo G2\_VTX00156
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Paraglomeromycetes
</td>
<td style="text-align:left;">
Paraglomerales
</td>
<td style="text-align:left;">
Paraglomeraceae
</td>
<td style="text-align:left;">
Paraglomus\_Alguacil11d Pa2\_VTX00239
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Diversisporales
</td>
<td style="text-align:left;">
Acaulosporaceae
</td>
<td style="text-align:left;">
Acaulospora\_SG07 Aca unk1\_VTX00230
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:left;">
Glomeromycota
</td>
<td style="text-align:left;">
Glomeromycetes
</td>
<td style="text-align:left;">
Glomerales
</td>
<td style="text-align:left;">
Glomeraceae
</td>
<td style="text-align:left;">
Glomus\_Kupea martinetugei symbiont\_VTX00204
</td>
</tr>
</tbody>
</table>

> **27** genera that are present in the reverse data and not in the
> forward.

Now the questions;

1.  Based on this is it enough that we continue the analyses henceforth
    separately for the forward and reverse reads?
2.  Can we merge the data sets ensuring we keep the taxa unique from
    both in the final dataset?
3.  Maybe we should investigate further the alpha and better diversities
    of the whole/unique datasets from the forward and reverse datasets?

#### Calculate Alpha diversity

The relative abundance is calculated in percentage, hence the total of
relative abundance per sample equals 100. In calculating the alpha
diversity, we need to set a threshold for rare biosphere/taxa, by
default its set at one as per the manual. We will go with that.

The `alphadiversity` function

``` r
alphaDiversity(temp1, siteInCol = F, taxhead = NULL, threshold = 1, percent = F, write = T)

temp1 <- temp[1:7, ]
as.integer(temp1[,2:4])

ncol(temp1)
cols <- c(1:10391)
temp1[cols] <- lapply(temp1[cols], as.numeric) #this worked to transform the data from character to numeric

temp1[1:7,1:10]

class(temp1)
str(temp1) %>% head()
```

### Working with the merged dataset

Using qiime2 we merged the forward and reverse feature tables and the
taxonomies as below;

        #====== Objective ========
        #merge the feature tables and taxonomies generated from the reverse and forward reads separately into one feature table and taxonomy file

        #====== merge the feature tables =========================

        qiime feature-table merge \
         --i-tables ${fwd}/yves_fwd_deblur_marjm_265_table.qza \
         --i-tables ${rev}/yves_rev_deblur_marjm_240_table.qza \
         --p-overlap-method error_on_overlapping_feature \
         --output-dir ${out}/fwd_rev_mgd_feature_table.qza


        #=============== merge the taxonomies ================

        qiime feature-table merge-taxa \
         --i-data ${fwd}/yves_fwd_deblur_marjm_265_repseqs_taxonomy.qza \
         --i-data ${rev}/yves_rev_deblur_marjm_240_repseqs_taxonomy.qza \
         --output-dir ${out}/fwd_rev_mgd_taxonomy.qza

        #============ export the merged files to use in R ================
        #=================================================================

        qiime tools export \
         --input-path ${out}/fwd_rev_mgd_feature_table_feb23.qza  \
         --output-path ${out}/fwd_rev_mgd_feature_table_fe23

        qiime tools export \
         --input-path ${out}/fwd_rev_mgd_taxonomy_feb23.qza  \
         --output-path ${out}/fwd_rev_mgd_taxonomy_feb23

        #======== convert exported files to tsv format =======

        biom convert \
         --input-fp ${out}/fwd_rev_mgd_feature_table_fe23/feature-table.biom \
        --output-fp ${out}/fwd_rev_mgd_feature_table_feb23.tsv \
        --to-tsv

Working with this data here to look at the what features it contains and
if check that all feature from the fwd and reverse tables are present in
the merged file.

``` r
#========== load the data

mgd_table <- read.table("C:/Users/BWaweru/OneDrive - CGIAR/Documents/Fellows/Yves_Armel/jan_2021/fwd_rev_mgd_feature_table_feb23.tsv", header = T)
mgd_taxo <- read.table("C:/Users/BWaweru/OneDrive - CGIAR/Documents/Fellows/Yves_Armel/jan_2021/fwd_rev_mgd_taxonomy_feb23.tsv", header = T, sep = '\t')
mgd_taxo$Confidence <- NULL #remove the colum with confidence values

#========= merge the two table

mgd_data <- merge(mgd_table, mgd_taxo, by="OTU_ID")
names(mgd_data)
rownames(mgd_data) <- mgd_data$OTU_ID #make the OTU_ID column the rownames
mgd_data$OTU_ID <- NULL #drop the OTU_ID as a column as they are now the rownames
names(mgd_data)
```

Now we use the *otuSUMMARY* function `otureport` to generate a report
from the merged dataset at genus level

``` r
#======= first at family level to see how many familes are reported

otuReport(mgd_data, siteInCol = T, taxhead = "taxonomy", platform = "qiime", pattern = ";", percent = F, taxlevel = "family") -> mgd_otu_report_fam

length(mgd_otu_report_fam[[1]]) # 12 families reported for the mgd data set, similar to the fwd and rev data set
```

    ## [1] 12

``` r
kable(mgd_otu_report_fam[[1]], caption = "Familes observed in the merged data set")
```

<table>
<caption>
Familes observed in the merged data set
</caption>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
</tr>
</tbody>
</table>

``` r
kable(mgd_otu_report_fam[[5]] %>% sort(decreasing = T), caption = "Familes present in the merged data and their count abundance")
```

<table>
<caption>
Familes present in the merged data and their count abundance
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae
</td>
<td style="text-align:right;">
1644871
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae
</td>
<td style="text-align:right;">
209581
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae
</td>
<td style="text-align:right;">
39600
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
38913
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae
</td>
<td style="text-align:right;">
21460
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA
</td>
<td style="text-align:right;">
16228
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Acaulosporaceae
</td>
<td style="text-align:right;">
3734
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Pacisporaceae
</td>
<td style="text-align:right;">
3219
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
2797
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;NA
</td>
<td style="text-align:right;">
2730
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Gigasporaceae
</td>
<td style="text-align:right;">
1161
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Archaeosporomycetes-&gt;Archaeosporales-&gt;Archaeosporaceae
</td>
<td style="text-align:right;">
123
</td>
</tr>
</tbody>
</table>

> The three familes **Glomeraceae**, **Claroideoglomeraceae**, and
> **Diversisporaceae** are still the most abundant in that order

``` r
#====== then at genus level to do the comparison

otuReport(mgd_data, siteInCol = T, taxhead = "taxonomy", platform = "qiime", pattern = ";", percent = F, taxlevel = "genus") -> mgd_otu_report

length(mgd_otu_report[[1]]) #156 genera reported in the merged data set
```

    ## [1] 156

``` r
#a look a the 10 most abundant genera in the merged data set

kable(mgd_otu_report[[5]] %>% sort(decreasing = T) %>% head(n=10L), caption = "10 most abundant genera in the merged data set")
```

<table>
<caption>
10 most abundant genera in the merged data set
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae-&gt;Glomus\_coremioides\_VTX00268
</td>
<td style="text-align:right;">
780285
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae-&gt;NA
</td>
<td style="text-align:right;">
756000
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae-&gt;NA
</td>
<td style="text-align:right;">
166217
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;NA-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
38913
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Diversisporales-&gt;Diversisporaceae-&gt;NA
</td>
<td style="text-align:right;">
36493
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae-&gt;Claroideoglomus\_SG07
Glo unk6\_VTX00193
</td>
<td style="text-align:right;">
22681
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;NA-&gt;NA
</td>
<td style="text-align:right;">
16228
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Claroideoglomeraceae-&gt;Claroideoglomus\_claroideum
</td>
<td style="text-align:right;">
13829
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Paraglomeromycetes-&gt;Paraglomerales-&gt;Paraglomeraceae-&gt;NA
</td>
<td style="text-align:right;">
13309
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi-&gt;Glomeromycota-&gt;Glomeromycetes-&gt;Glomerales-&gt;Glomeraceae-&gt;Glomus\_Alguacil11d
Glo G13\_VTX00067
</td>
<td style="text-align:right;">
10972
</td>
</tr>
</tbody>
</table>

> Now we look a to see if all the genera from the forward and reverse
> reads are present in the merged dataset

We generate a dataframe we can use to do the intersection

``` r
mgd_df <- data.frame(do.call(rbind, strsplit(as.vector(mgd_otu_report[[1]]), "->", fixed=TRUE)))

#intersection between the fwd and the mgd

length(intersect(fwd_df$X6, mgd_df$X6)) #118
```

    ## [1] 118

``` r
length(intersect(mgd_df$X6, rev_df$X6)) #125
```

    ## [1] 125

``` r
# we don't expect any differences, all the genera found in the fwd/rev data sets should also be in the merged data

length(setdiff(fwd_df$X6, mgd_df$X6)) #zero, all genera found in forward are also in the merged data
```

    ## [1] 0

``` r
length(setdiff(rev_df$X6, mgd_df$X6)) #zero, all genera found in reverse are also in the merged data
```

    ## [1] 0

Essentially the merged data set is a concatenation of the individual
feature and taxa tables hence the merged data set is simply a data
object with two combined hence the result.

Ideally however, we want to be able to analyze our data set such that a
**feature/OTU** is identified or is a result of the *forward and reverse
read pair*. In qiime as of now we cannot trace which read gave rise to
which feature.

We need to understand this data better, can we trace which read pair
generates which taxa?
