Analysis of raw reads with BWA and BLAST
================
Bernice Waweru
Tue 22, Jun 2021

-   [Taxonomic assignemnts of reads with BWA and
    BLAST](#taxonomic-assignemnts-of-reads-with-bwa-and-blast)
    -   [March 11 2021](#march-11-2021)
    -   [April 8 2021](#april-8-2021)
        -   [Read length distribution of database
            reads](#read-length-distribution-of-database-reads)
-   [Session Information](#session-information)

## Taxonomic assignemnts of reads with BWA and BLAST

##### March 11 2021

Set up for mapping with *bwa-mem*, first with the whole set of files
which was taking too long, then repeat the mapping with sequences
filtered that above 300 bp, not very many per sample. Moreover even for
these that has a good number of sequences, most were not properly paired
Maybe we should set minimum read lengths to those we used to generate
feature in qiime2, 240 as the minimum as that was the minimum for the
reverse reads.

The **maarjAM.5.fasta** file has *5,934* seq in total, *5,551* of them
with lengths greater than 300bp.

In terms of tracing, we can filter as below; 1. Filter reads that mapped
2. Filter reads that are properly paired in mapping 3. Filter reads
where mates mapped to different positions 4. Find the taxonomic
classification of db hits from the above two filters and how to move on
with those. 5. Check on the mapping quality of the the reads that
mapped, 6. Make a db of the hits and use those for blast, check on the
percentage identities and the query coverage, OR is there a way to
control for this during the mapping?

    i. *-k* seed length?
    ii. *-c* discard a MEM if it has more than INT hits on the genomes
    iii. *-A* Matching score 
    iv. *T* don't output with alignment score lower than INT

We re-do the mapping with reads filtered at minimum length of 240bp. We
also use samtools to look at the statistics after mapping to see how
well it worked, below is the are the commandswe used;

    # ========== define path variables ==================================

    ref='/home/bngina/Fellows/Yves_Tchiechoua/maarJAM_db/maarjAM.5.fasta'

    fastq_dir='/home/bngina/Fellows/Yves_Tchiechoua/orig_data'


    bwa_out='/home/bngina/Fellows/Yves_Tchiechoua/bwa_out'


    # ========= first step is to index the reference db =================

    bwa index ${ref}


    # =========== using a loop to align the files =======================

    for R1 in ${fastq_dir}/S*R1_001.fastq.gz;\
    do echo ${R1};\
    R2=$(echo $R1 | sed 's/R1/R2/g');\
    out_name=$(echo $R1 | cut -f 7 -d "/" | cut -f 1 -d "." | sed 's/_R1_001//g');\
    echo $R2 ;\
    echo $out_name ;\

    #============ filter for reads with minimum length 300bp ============
    zcat ${R1} | seqkit seq -m 240 > ${fastq_dir}/${out_name}_R1_240.fq ;\
    zcat ${R2} | seqkit seq -m 240 > ${fastq_dir}/${out_name}_R2_240.fq ;\

    # ======== do the mapping ==================
    bwa mem -t 4 ${ref} ${fastq_dir}/${out_name}_R1_240.fq ${fastq_dir}/${out_name}_R2_240.fq  >  ${bwa_out}/${out_name}_240.sam ;\

    # ==== use samtools to sort by name then save the output as bam file ========
    samtools sort -n -o ${bwa_out}/${out_name}_srtd_240.bam -O bam -@ 4 ${bwa_out}/${out_name}_240.sam ;\

    #========== filter for reads that mapped; get stats before and after in a text file ===================
    echo ${out_name} >> ${bwa_out}/${out_name}_stats_240.txt ;\

    samtools flagstat ${bwa_out}/${out_name}_srtd_240.bam >> ${bwa_out}/${out_name}_stats_240.txt ;\

    # ======== filtering that reads that mapped ====================================
    samtools view -b -F 4 ${bwa_out}/${out_name}_srtd_240.bam > ${bwa_out}/${out_name}_srtd_mpd_240.bam ;\

    echo ${bwa_out}'/'${out_name}'_srtd_mpd.bam' >> ${bwa_out}/${out_name}_stats_240.txt ;\

    samtools flagstat ${bwa_out}/${out_name}_srtd_mpd_240.bam >> ${bwa_out}/${out_name}_stats_240.txt ;\

    done

After above, we look at one of the text files with the stats, for
example we look at *S4* that had the highest number of reads.

    less -S bwa_out/S4_S3_L001_stats_240.txt

Gives us a look at the file;

    S4_S3_L001

    3730788 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    3591614 + 0 mapped (96.27% : N/A)
    3730788 + 0 paired in sequencing
    1864318 + 0 read1
    1866470 + 0 read2
    0 + 0 properly paired (0.00% : N/A)
    3455894 + 0 with itself and mate mapped
    135720 + 0 singletons (3.64% : N/A)
    3411860 + 0 with mate mapped to a different chr
    444469 + 0 with mate mapped to a different chr (mapQ>=5)

    /home/bngina/Fellows/Yves_Tchiechoua/bwa_out/S4_S3_L001_srtd_mpd.bam

    3591614 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 secondary
    0 + 0 supplementary
    0 + 0 duplicates
    3591614 + 0 mapped (100.00% : N/A)
    3591614 + 0 paired in sequencing
    1804268 + 0 read1
    1787346 + 0 read2
    0 + 0 properly paired (0.00% : N/A)
    3455894 + 0 with itself and mate mapped
    135720 + 0 singletons (3.78% : N/A)
    3411860 + 0 with mate mapped to a different chr
    444469 + 0 with mate mapped to a different chr (mapQ>=5)

WE see that we had **3,730,788** reads in total that passed QC,
**3,591,614 (96.27% of total)** mapped to our reference. *1,804,268* are
Read 1 and *1,787,346* are Read 2. *3,455,894* had mate pairs and both
mapped. 3,411,860 with mate mapped to a different chromosome or
position, in our case to a difference taxonomic sequence. Further only
444,469 with the mate mapped to a different sequence with a mapping
quality higher than 5(&gt;5).

We see that a majority of the read pairs did not map to the same
taxonomic sequence, looking at this from the bwa output with the below;

    samtools view bwa_out/S4_S3_L001_srtd_240.bam | less -S

We see the output;

    QNAME                                           FLAG    RNAME           POS     MAPQ    CIAGR               MRNM            MPOS    ISIZE
    M03021:27:000000000-C8YMF:1:1101:1914:12854     113     HQ323483        82      0       281M20S             AB546438        52      0       ACGGGGAGGTAGTGACAATAAATAA
    M03021:27:000000000-C8YMF:1:1101:1914:12854     129     AJ276079        276     0       300M                AB365857        237     0       GGGAATTAGGGTTCGATTCCGGAGA
    M03021:27:000000000-C8YMF:1:1101:1930:13801     81      HE798939        10      0       300M                AJ563863        232     0       GGTAATTCCAGCTACAATAGCGTAT
    M03021:27:000000000-C8YMF:1:1101:1930:13801     129     EU340308        220     13      301M                AF485877        192     0       GGAATGAGTACAATTTAAATCTCTT
    M03021:27:000000000-C8YMF:1:1101:1944:12381     113     EU340308        198     26      293M7S              DQ336523        64      0       GGCTCTTTCGGGTTTAGTAATTGGA
    M03021:27:000000000-C8YMF:1:1101:1944:12381     161     AB555671        295     0       227M74S             HE798937        164     0       GGATAGAGGCCTACCATGGTGGTAA
    M03021:27:000000000-C8YMF:1:1101:1956:11727     145     JN090188        388     0       203S98M             AB365857        158     0       TAGTCTATCAACTAGTACGACAGCA
    M03021:27:000000000-C8YMF:1:1101:1971:13875     97      HE798808        111     0       47M1D3M1I178M72S    JN131598        89      0       TGTGTCACT
    M03021:27:000000000-C8YMF:1:1101:1971:13875     145     HE615037        403     0       107S107M1D87M       FJ831568        212     0       TTTATTCGCCGTTACTT
    M03021:27:000000000-C8YMF:1:1101:1980:11096     65      AB546438        19      0       294M7S              FR821529        150     0       CGGGTAACGGGGTGTTAGGGCACGA
    M03021:27:000000000-C8YMF:1:1101:1980:11096     177     EU340308        218     6       300M                FJ831602        511     0       TTGGAAGAAGTACAATTTAAGGTTC
    M03021:27:000000000-C8YMF:1:1101:1985:13603     113     FN869855        224     0       93S208M             FR773152        478     0       TACCTTCTTAGTCTTTCTTTTGTTT
    M03021:27:000000000-C8YMF:1:1101:1987:12515     97      FJ831568        152     17      301M                AJ496107        174     0       ACCCAATCCCGACACGGGGAGGTAG
    M03021:27:000000000-C8YMF:1:1101:1987:12515     145     FN263141        155     0       300M                EU340308        359     0       AAAACCGCTATGTCATTAATTTGGT
    M03021:27:000000000-C8YMF:1:1101:1989:12133     97      FR750204        336     0       143M157S            FJ009602        445     0       GGTTTTAACGGGTAACG
    M03021:27:000000000-C8YMF:1:1101:1991:10789     129     AB555671        305     0       176M125S            DQ164810        139     0       CTACCATGGTGGTAACG
    M03021:27:000000000-C8YMF:1:1101:1992:14275     113     JF683559        84      0       150S150M            HE615079        40      0       ACGTCTGATCTTGCATA
    M03021:27:000000000-C8YMF:1:1101:1993:14109     81      FR821553        57      0       301M                DQ396751        122     0       TTAGGGCACGACACCGGAGAGGGAG
    M03021:27:000000000-C8YMF:1:1101:1993:14109     161     AY493667        261     0       301M                FR847991        162     0       CTACCATGGTGGTAACGGGTAACGG
    M03021:27:000000000-C8YMF:1:1101:1995:10596     81      GU322408        208     0       43S257M             EU340308        359     0       GTTCAGGCGTCTCGTGGGCTCGGAG
    M03021:27:000000000-C8YMF:1:1101:1995:10596     161     DQ085198        479     0       254M47S             FR847991        167     0       GTTTAAAGCAGGCACACGCTTGAAT
    M03021:27:000000000-C8YMF:1:1101:1995:14821     113     HE798939        77      0       77S177M1I46M        AY129612        67      0       TTTTGTTTTGTTCAAGC
    M03021:27:000000000-C8YMF:1:1101:1997:11165     81      FR773152        498     0       270M                EU417623        251     0       AATGAGTACAATTTAAAGCCCTTAA
    M03021:27:000000000-C8YMF:1:1101:1997:11165     161     EU340318        222     0       270M                FJ831568        82      0       AATGAGTACAATTTAAATCCCTTAA
    M03021:27:000000000-C8YMF:1:1101:2008:13400     81      DQ396756        126     0       301M                EU340308        202     0       CAAGGAAGGCAGCAGGCGCGCAAAT
    M03021:27:000000000-C8YMF:1:1101:2008:13400     129     DQ396756        56      0       281M20S             AB365822        147     0       TGGTTTTAACGGGTAACGGGGTGTT
    M03021:27:000000000-C8YMF:1:1101:2017:12071     113     HE615079        24      0       111S190M            FR728614        160     0       CCAGAAGCACATATCCT
    M03021:27:000000000-C8YMF:1:1101:2024:12862     81      HE798939        56      0       209M1D92M           EU340293        184     0       AAGCTCGGAGTTGAATT
    M03021:27:000000000-C8YMF:1:1101:2024:12862     137     EU340308        310     0       232M1D61M8S         =               310     0       ATATTAAAGTTGTTGCAGTTAAAAA
    M03021:27:000000000-C8YMF:1:1101:2026:13382     65      FN869785        188     0       8S235M              AM746143        158     0       CTTTCGGATCTCGTAATTGGAATGA
    M03021:27:000000000-C8YMF:1:1101:2026:13382     145     FN869790        190     0       10S233M             DQ164810        120     0       CTGTCGGATCTCGTAATTGGAATGA
    M03021:27:000000000-C8YMF:1:1101:2030:10684     81      AJ699064        371     0       244M                FN646000        63      0       CCGGAGAGGGAGCCTGAGAAACGGC
    M03021:27:000000000-C8YMF:1:1101:2030:10684     161     DQ336444        244     0       246M55S             DQ371679        156     0       CCGGAGAGGGAGCCTGAGAAACGGC
    M03021:27:000000000-C8YMF:1:1101:2031:11716     113     DQ336486        237     0       289M12S             HE798939        139     0       CACGACACCGGAGAGGGAGCCTGAG
    M03021:27:000000000-C8YMF:1:1101:2031:11716     161     AB546401        7       0       301M                HE614962        245     0       CCATGGTGGTAACGGGTAACGGGGT
    M03021:27:000000000-C8YMF:1:1101:2031:14432     97      AB546438        24      0       301M                HQ656911        111     0       AACGGGGTGTTAGGGCACGACACCG
    M03021:27:000000000-C8YMF:1:1101:2031:14432     145     EU340308        111     8       301M                DQ396726        115     0       AACGGCTGCCACATCCAAGGATGGC
    M03021:27:000000000-C8YMF:1:1101:2042:14721     97      EU417585        413     0       188M75S             HE614949        44      0       CTTTCGGGTTTCGTAATTGGGATGA

The headers from the output are explained in the [bwa
manual](http://bio-bwa.sourceforge.net/bwa.shtml) as below;

    Col   Field   Description
    1       QNAME     Query (pair) NAME
    2       FLAG      bitwise FLAG
    3       RNAME     Reference sequence NAME
    4       POS     1-based leftmost POSition/coordinate of clipped sequence
    5       MAPQ      MAPping Quality (Phred-scaled)
    6       CIAGR     extended CIGAR string
    7       MRNM      Mate Reference sequence NaMe (‘=’ if same as RNAME)
    8       MPOS      1-based Mate POSistion
    9       ISIZE     Inferred insert SIZE
    10    SEQ       query SEQuence on the same strand as the reference
    11    QUAL    query QUALity (ASCII-33 gives the Phred base quality)
    12    OPT       variable OPTional fields in the format TAG:VTYPE:VALUE

Looking at some informative columns, **column 3** shows where the **read
mapped**, with the **mapping quality** at **column 5** and information
about where the **mate mapped** in **column 7**. We see for a majority
of the reads in our view, **the mapping quality is 0, and the mate
paired to a different sequence**.

Arguably, the sequences don’t differ very much. it could be just a one
or two base pair substitution that differentiates the reference
taxonomic sequences, so we need to look at this further. Maybe we can do
a blast, so that we look at the read 1 and read 2 two alignments and
their classification to see if there is a big difference.

Pursuing blast, below are the command line steps undertaken;

    # ===== first we use the files where we filterd for min seq length of 240 bp ==========
    # ===== we then interleave thos files, the tweak the read names


    # ===== using  reformat.sh from bbmap suite of tools to interleave
    module load  bbmap/38.67

    # ====== dir to store interleaved fasta files

    mkdir -p /home/bngina/Fellows/Yves_Tchiechoua/orig_data/fasta_inter

    fasta_dir='/home/bngina/Fellows/Yves_Tchiechoua/orig_data/fasta_inter'


    for read in ${fastq_dir}/S10*R1_001.fastq.gz;\
    do out_name=$(echo $read | cut -f 7 -d "/" | cut -f 1 -d "." | sed 's/_R1_001//g');\
       R1=${fastq_dir}/${out_name}_R1_240.fq ;\
       R2=${fastq_dir}/${out_name}_R2_240.fq ;\
       echo $R2 ;\
       echo $R1 ;\
       echo $out_name ;\
       reformat.sh in1=${R1} in2=${R2}  out=${fasta_dir}/${out_name}_240.fasta ;\
    # ===== tweaking the sequnce headers so as to track the read 1 and 2 in the blast output
       cat ${fasta_dir}/${out_name}_240.fasta | sed 's/ 1:N:/_1:N:/g' | sed 's/ 2:N:/_2:N:/g' > ${fasta_dir}/${out_name}_240_twkd.fasta ;\
    done

    #index the maarjam reference database for blast

    makeblastdb \
     -in ${ref} \
     -dbtype nucl


    # ===== run the blast with the indexed maarjAM db ======


    # ===== the output dir for the blast

    mkdir -p /home/bngina/Fellows/Yves_Tchiechoua/blast_out

    out_dir='/home/bngina/Fellows/Yves_Tchiechoua/blast_out'

    for file in ${fasta_dir}/*twkd.fasta ;\
    do echo ${file};\
    out=$(echo ${file} | cut -f 8 -d "/" | sed 's/_twkd.fasta//g' ) ;\
    echo ${out};\
    blastn\
     -num_threads 4 \
     -query ${file}\
     -task blastn\
     -db ${ref}\
     # ===== because we needed extrax information, i.e the query cover percentages, we run blast twice
     # ===== the first to get the normal blast output format 6, then the next to get the query cover
     -outfmt "6 qcovs" \
     -evalue 0.001\
     -max_target_seqs 2\
    # -out ${out_dir}/${out}_qcovs.out ;\
    # ===== we concantenate the blast out puts, the normal and the querry cover
    pr -mts ${out_dir}/${out}.out ${out_dir}/${out}_qcovs.out > ${out_dir}/${out}_blast_f.out;\
    # ===== we follow a few steps in order to extract the taxonmic information then add it to our final file

    # ==== sort the blast output based on the second column with the tax sequence IDs
    sort -k 2,2 ${out_dir}/${out}_blast_f.out > ${out_dir}/${out}_blast_f_srtd_tax_seq_ids.out;\

    # ===== extract the sorted Ids and store them in a text file
    awk '{ print $2 }' ${out_dir}/${out}_blast_f_srtd_tax_seq_ids.out > ${out_dir}/${out}_tax_seq_ids.txt ;\

    # ===== use join to get the corresponding ids in the text file from the maarjam taxonomy file and store them in a text file
    join ${out_dir}/${out}_tax_seq_ids.txt /home/bngina/Fellows/Yves_Tchiechoua/maarjAM.id_to_taxonomy.5_rev.txt > ${out_dir}/${out}_srtd_tax_ids_with_taxonomy.txt;\

    # ===== join the sorted blast output with the out put from join with the taxonmoy information and then sort them based on the query id before saving the final output in a text file
    pr -mts ${out_dir}/${out}_blast_f_srtd_tax_seq_ids.out ${out_dir}/${out}_srtd_tax_ids_with_taxonomy.txt | sort -k 1,1 > ${out_dir}/${out}_blast_f_with_taxonomy_srtd_query_id.txt;\
    done

And now a look at one of the sample file final outputs

    qseqid                                                  sseqid        pident   length    mismatch    gapopen     qstart  qend    sstart  send    evalue  bitscore   qcovs   salltitles
    M03021:27:000000000-C8YMF:1:1101:10000:2562_1:N:0:2     FJ009603        98.333  300     5       0       1       300     59      358     6.44e-148       518     99      FJ009603 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_intraradices_VTX00114
    M03021:27:000000000-C8YMF:1:1101:10000:2562_1:N:0:2     HE576805        98.662  299     4       0       1       299     34      332     5.29e-149       522     99      HE576805 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Claroideoglomeraceae;Claroideoglomus_Alguacil12b GLO G1_VTX00056
    M03021:27:000000000-C8YMF:1:1101:10000:2562_2:N:0:2     DQ396709        96.949  295     9       0       6       300     493     199     8.92e-140       491     98      DQ396709 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_PF14_VTX00083
    M03021:27:000000000-C8YMF:1:1101:10000:2562_2:N:0:2     FN869738        96.610  295     10      0       6       300     469     175     1.09e-138       488     98      FN869738 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_Glo G2_VTX00083
    M03021:27:000000000-C8YMF:1:1101:10015:5254_1:N:0:2     DQ396700        95.973  298     11      1       1       298     396     100     1.61e-136       480     99      DQ396700 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_PF22_VTX00117
    M03021:27:000000000-C8YMF:1:1101:10015:5254_1:N:0:2     DQ396756        95.638  298     13      0       1       298     397     100     5.63e-136       479     99      DQ396756 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_PF9_VTX00117
    M03021:27:000000000-C8YMF:1:1101:10015:5254_2:N:0:2     DQ396756        91.864  295     24      0       1       295     56      350     1.09e-119       425     98      DQ396756 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_PF9_VTX00117
    M03021:27:000000000-C8YMF:1:1101:10015:5254_2:N:0:2     DQ396774        91.864  295     24      0       1       295     56      350     1.09e-119       425     98      DQ396774 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_PF21_VTX00117
    M03021:27:000000000-C8YMF:1:1101:10024:5136_1:N:0:2     HE576898        100.000 296     0       0       1       296     35      330     8.23e-153       535     100     HE576898 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_Alguacil12b GLO G14_VTX00191
    M03021:27:000000000-C8YMF:1:1101:10024:5136_1:N:0:2     X86687  99.662  296     1       0       1       296     336     631     3.50e-151       529     100     X86687 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_epigaea_VTX00061
    M03021:27:000000000-C8YMF:1:1101:10024:5136_2:N:0:2     AJ276088        97.973  296     6       0       1       296     631     336     1.14e-144       508     100     AJ276088 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_epigaea_VTX00061
    M03021:27:000000000-C8YMF:1:1101:10024:5136_2:N:0:2     HE576898        97.973  296     6       0       1       296     330     35      1.14e-144       508     100     HE576898 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_Alguacil12b GLO G14_VTX00191
    M03021:27:000000000-C8YMF:1:1101:10029:2552_1:N:0:2     HE615082        98.456  259     4       0       3       261     182     440     2.35e-127       450     99      HE615082 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_Torrecillas12b Div1_VTX00062
    M03021:27:000000000-C8YMF:1:1101:10029:2552_1:N:0:2     HQ657105        98.456  259     4       0       3       261     120     378     2.35e-127       450     99      HQ657105 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_Diversispora1_VTX00062
    M03021:27:000000000-C8YMF:1:1101:10029:2552_2:N:0:2     HE615082        98.069  259     5       0       1       259     440     182     9.98e-126       444     99      HE615082 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_Torrecillas12b Div1_VTX00062
    M03021:27:000000000-C8YMF:1:1101:10029:2552_2:N:0:2     HQ657105        98.069  259     5       0       1       259     378     120     9.98e-126       444     99      HQ657105 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_Diversispora1_VTX00062
    M03021:27:000000000-C8YMF:1:1101:10030:2280_1:N:0:2     FJ009619        98.851  261     3       0       1       261     532     792     1.58e-129       457     100     FJ009619 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_geosporum_VTX00065
    M03021:27:000000000-C8YMF:1:1101:10030:2280_1:N:0:2     Z14007  98.851  261     3       0       1       261     730     990     1.58e-129       457     100     Z14007 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_mosseae_VTX00067
    M03021:27:000000000-C8YMF:1:1101:10030:2280_2:N:0:2     FJ009619        98.851  261     3       0       1       261     792     532     1.58e-129       457     100     FJ009619 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_geosporum_VTX00065
    M03021:27:000000000-C8YMF:1:1101:10030:2280_2:N:0:2     Z14007  98.851  261     3       0       1       261     990     730     1.58e-129       457     100     Z14007 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_mosseae_VTX00067
    M03021:27:000000000-C8YMF:1:1101:10031:3880_1:N:0:2     FJ831643        97.674  301     7       0       1       301     410     110     9.56e-146       511     100     FJ831643 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_NF29_VTX00062
    M03021:27:000000000-C8YMF:1:1101:10031:3880_1:N:0:2     JQ811206        97.674  301     7       0       1       301     635     335     9.56e-146       511     100     JQ811206 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Diversisporaceae;Diversispora_sp. 8479.a.1_VTX00062
    M03021:27:000000000-C8YMF:1:1101:10031:3880_2:N:0:2     JN131587        96.656  299     10      0       1       299     34      332     7.33e-141       495     99      JN131587 Fungi;Glomeromycota;Paraglomeromycetes;Paraglomerales;Paraglomeraceae;Paraglomus_majewskii_VTX00335
    M03021:27:000000000-C8YMF:1:1101:10031:3880_2:N:0:2     JN131588        96.656  299     10      0       1       299     34      332     7.33e-141       495     99      JN131588 Fungi;Glomeromycota;Paraglomeromycetes;Paraglomerales;Paraglomeraceae;Paraglomus_majewskii_VTX00335
    M03021:27:000000000-C8YMF:1:1101:10036:3813_2:N:0:2     JN131594        99.142  233     2       0       6       238     547     779     6.87e-116       412     78      JN131594 Fungi;Glomeromycota;Paraglomeromycetes;Paraglomerales;Paraglomeraceae;Paraglomus_majewskii_VTX00335
    M03021:27:000000000-C8YMF:1:1101:10036:3813_2:N:0:2     JN131595        99.142  233     2       0       6       238     547     779     6.87e-116       412     78      JN131595 Fungi;Glomeromycota;Paraglomeromycetes;Paraglomerales;Paraglomeraceae;Paraglomus_majewskii_VTX00335
    M03021:27:000000000-C8YMF:1:1101:10045:4540_1:N:0:2     JF414171        95.695  302     11      1       1       300     601     300     1.32e-137       484     100     JF414171 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Acaulosporaceae;Acaulospora_sp._VTX00328
    M03021:27:000000000-C8YMF:1:1101:10045:4540_1:N:0:2     JF414178        95.695  302     11      1       1       300     602     301     1.32e-137       484     100     JF414178 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Acaulosporaceae;Acaulospora_sp._VTX00328
    M03021:27:000000000-C8YMF:1:1101:10045:4540_2:N:0:2     GQ140613        87.000  300     39      0       1       300     35      334     8.97e-102       365     99      GQ140613 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Acaulosporaceae;Acaulospora_sp._VTX00249
    M03021:27:000000000-C8YMF:1:1101:10045:4540_2:N:0:2     JN252439        87.000  300     39      0       1       300     35      334     8.97e-102       365     99      JN252439 Fungi;Glomeromycota;Glomeromycetes;Diversisporales;Acaulosporaceae;Acaulospora_Early-19_VTX00026
    M03021:27:000000000-C8YMF:1:1101:10047:2661_2:N:0:2     FN869837        97.449  196     5       0       8       203     188     383     1.84e-91        331     65      FN869837 Fungi;Glomeromycota;Paraglomeromycetes;Paraglomerales;Paraglomeraceae;Paraglomus_Para1_VTX00335
    M03021:27:000000000-C8YMF:1:1101:10047:2661_2:N:0:2     HE615075        96.939  196     6       0       8       203     186     381     2.25e-90        327     65      HE615075 Fungi;Glomeromycota;Paraglomeromycetes;Paraglomerales;Paraglomeraceae;Paraglomus_Torrecillas12b Para1_VTX00335
    M03021:27:000000000-C8YMF:1:1101:10053:3919_1:N:0:2     FJ831624        100.000 301     0       0       1       301     452     152     1.62e-155       544     100     FJ831624 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_NF21_VTX00064
    M03021:27:000000000-C8YMF:1:1101:10053:3919_1:N:0:2     FJ831626        100.000 301     0       0       1       301     452     152     1.62e-155       544     100     FJ831626 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_NF21_VTX00064
    M03021:27:000000000-C8YMF:1:1101:10053:3919_2:N:0:2     HE615054        96.000  300     12      0       1       300     34      333     1.09e-138       488     99      HE615054 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_Torrecillas12b Glo G9_VTX00064
    M03021:27:000000000-C8YMF:1:1101:10053:3919_2:N:0:2     HE615055        96.000  300     12      0       1       300     34      333     1.09e-138       488     99      HE615055 Fungi;Glomeromycota;Glomeromycetes;Glomerales;Glomeraceae;Glomus_Torrecillas12b Glo G9_VTX00064

We sorted the above final file based on the first column, the query/read
ID, so we can see tht read 1 and read 2 are sequential in the list. We
filtered the blast output to the top two best hits, thats why we only
see two alignments per read. In [qiime2 consensus blast feature
classifier](https://docs.qiime2.org/2021.2/plugins/available/feature-classifier/classify-consensus-blast/),
the algorithm works by first discarding anything that has a less than
80% percentage identity, then a *qcov\_hsp\_perc* parameter of not less
than 80%, and for consensus, 50% of the assignments need to match the
top hit for it to be considered a consensus assignment, otherwise its
unassigned.

Next is to think of a way to filter our file so we can extract sequences
where the taxonomic classification is the same for read 1 and 2,
although for most it is different as we have observed.

#### April 8 2021

From the above results, its clear we need to understand the data
structure behind the [*maarjAM
database*](https://maarjam.botany.ut.ee/?action=about). Further,we need
to understand how the virtual taxa(VTX) are assigned. We notice that for
some assignments, the read one and 2 are assigned to the same final
virtual taxa, but the values of assignments beforer the VTX number os
different and vice versa. How is this so?

We need to read the [paper mentioned on the
webpage](https://nph.onlinelibrary.wiley.com/doi/epdf/10.1111/j.1469-8137.2010.03334.x)
that explains the database carefully, and understand the data structure.

Also we agreed to look into other possible fungal databases that can be
used to correctly classify the AMF sequences from our study. In
particular databases that are regularly maintained and recently updated.
In this regard, came across two that can be considered;

1.  The [Protist Ribosomal Reference database i.e the PR2
    databse](https://github.com/pr2database/pr2database). Its aim is to
    provide a reference database of carefully annotated 18S rRNA
    sequences using eight unique taxonomic fields (from kingdom to
    species).Data files for use can be downloaded from
    [here](https://github.com/pr2database/pr2database/releases). The
    latest version, 4.13.0, was updated just last month.
2.  The [UNITE](https://unite.ut.ee/).It is a database and sequence
    management environment centered on the eukaryotic nuclear ribosomal
    ITS region. Data files for use can be downloaded from
    [here](https://unite.ut.ee/repository.php). Latest release was in
    early 2020.

##### Read length distribution of database reads

We also needed to check the length distribution of the reads within the
maarjAM database fasta file we have. To do that, we first count the
lengths if the reds within the unix command line environment, then
transfer the file to plot a histogram with R.

We calculate the read lengths using the below code

    #!/bin/bash

    ref_dir='/home/bngina/Fellows/Yves_Tchiechoua/maarJAM_db'

    # ===== count the number of reads and their lengths

    for file in ${ref_dir}/*.fasta ;\
    do echo ${file};\
    echo ${file} 'total number of scaffolds' $(grep -c '^>' ${file}) $(awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ${file}) | sed "s/>/\n/g" >> ref_scaffolds_stats.txt ;\
    done

We upload the file into R;

``` r
# ===== we read the table skipping the first line that gives info on the file and number of reads (5,934)

read.table(file = "data/ref_scaffolds_stats.txt") -> ref_read_lens

# ===== checking it

head(ref_read_lens)

# ===== we rename the columns

colnames(ref_read_lens) <- c("read_name", "length")
names(ref_read_lens)

# ==== now to plot a basic histogram
hist(ref_read_lens$length)
```

![](Analysis-of-raw-reads-with-BWA-and-BLAST_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

We use ggplot to get a better plot;

``` r
# ===== using ggplot2

require(ggplot2)

ggplot(ref_read_lens, aes(x=length)) + geom_histogram(binwidth = 100,color="black", fill="lightblue") + theme_grey() + 
  labs(title = " Distribution of read lengths within the maarjAM database fasta file", x="Read length", y="Number of reads")
```

![](Analysis-of-raw-reads-with-BWA-and-BLAST_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

We see that most of the reads are 500bp in length, not the full \~795 bp
of the variable region target by the primer. Its also interesting though
that we do have a few sequences longer than 1000bp.

## Session Information

``` r
devtools::session_info()
```

    ## - Session info ---------------------------------------------------------------
    ##  setting  value                       
    ##  version  R version 4.0.3 (2020-10-10)
    ##  os       Windows 10 x64              
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United States.1252  
    ##  ctype    English_United States.1252  
    ##  tz       Africa/Nairobi              
    ##  date     2021-06-22                  
    ## 
    ## - Packages -------------------------------------------------------------------
    ##  package     * version date       lib source        
    ##  assertthat    0.2.1   2019-03-21 [2] CRAN (R 4.0.3)
    ##  cachem        1.0.5   2021-05-15 [2] CRAN (R 4.0.5)
    ##  callr         3.7.0   2021-04-20 [2] CRAN (R 4.0.5)
    ##  cli           2.5.0   2021-04-26 [2] CRAN (R 4.0.5)
    ##  colorspace    2.0-1   2021-05-04 [2] CRAN (R 4.0.5)
    ##  crayon        1.4.1   2021-02-08 [2] CRAN (R 4.0.5)
    ##  DBI           1.1.1   2021-01-15 [2] CRAN (R 4.0.3)
    ##  desc          1.3.0   2021-03-05 [2] CRAN (R 4.0.5)
    ##  devtools      2.4.2   2021-06-07 [2] CRAN (R 4.0.3)
    ##  digest        0.6.27  2020-10-24 [1] CRAN (R 4.0.3)
    ##  dplyr         1.0.6   2021-05-05 [2] CRAN (R 4.0.5)
    ##  ellipsis      0.3.2   2021-04-29 [2] CRAN (R 4.0.5)
    ##  evaluate      0.14    2019-05-28 [2] CRAN (R 4.0.3)
    ##  fansi         0.5.0   2021-05-25 [2] CRAN (R 4.0.5)
    ##  farver        2.1.0   2021-02-28 [2] CRAN (R 4.0.5)
    ##  fastmap       1.1.0   2021-01-25 [2] CRAN (R 4.0.5)
    ##  fs            1.5.0   2020-07-31 [2] CRAN (R 4.0.3)
    ##  generics      0.1.0   2020-10-31 [2] CRAN (R 4.0.3)
    ##  ggplot2     * 3.3.4   2021-06-16 [2] CRAN (R 4.0.3)
    ##  glue          1.4.2   2020-08-27 [2] CRAN (R 4.0.3)
    ##  gtable        0.3.0   2019-03-25 [2] CRAN (R 4.0.3)
    ##  highr         0.9     2021-04-16 [2] CRAN (R 4.0.5)
    ##  htmltools     0.5.1.1 2021-01-22 [2] CRAN (R 4.0.5)
    ##  knitr         1.33    2021-04-24 [2] CRAN (R 4.0.5)
    ##  labeling      0.4.2   2020-10-20 [2] CRAN (R 4.0.3)
    ##  lifecycle     1.0.0   2021-02-15 [2] CRAN (R 4.0.5)
    ##  magrittr      2.0.1   2020-11-17 [2] CRAN (R 4.0.3)
    ##  memoise       2.0.0   2021-01-26 [2] CRAN (R 4.0.5)
    ##  munsell       0.5.0   2018-06-12 [2] CRAN (R 4.0.3)
    ##  pillar        1.6.1   2021-05-16 [2] CRAN (R 4.0.5)
    ##  pkgbuild      1.2.0   2020-12-15 [2] CRAN (R 4.0.3)
    ##  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.0.3)
    ##  pkgload       1.2.1   2021-04-06 [2] CRAN (R 4.0.5)
    ##  prettyunits   1.1.1   2020-01-24 [2] CRAN (R 4.0.3)
    ##  processx      3.5.2   2021-04-30 [2] CRAN (R 4.0.5)
    ##  ps            1.6.0   2021-02-28 [2] CRAN (R 4.0.5)
    ##  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.0.3)
    ##  R6            2.5.0   2020-10-28 [2] CRAN (R 4.0.3)
    ##  remotes       2.4.0   2021-06-02 [2] CRAN (R 4.0.5)
    ##  rlang         0.4.11  2021-04-30 [2] CRAN (R 4.0.5)
    ##  rmarkdown     2.9     2021-06-15 [2] CRAN (R 4.0.3)
    ##  rprojroot     2.0.2   2020-11-15 [2] CRAN (R 4.0.3)
    ##  scales        1.1.1   2020-05-11 [2] CRAN (R 4.0.3)
    ##  sessioninfo   1.1.1   2018-11-05 [2] CRAN (R 4.0.3)
    ##  stringi       1.6.2   2021-05-17 [2] CRAN (R 4.0.3)
    ##  stringr       1.4.0   2019-02-10 [2] CRAN (R 4.0.3)
    ##  testthat      3.0.2   2021-02-14 [2] CRAN (R 4.0.5)
    ##  tibble        3.1.2   2021-05-16 [2] CRAN (R 4.0.5)
    ##  tidyselect    1.1.1   2021-04-30 [2] CRAN (R 4.0.5)
    ##  usethis       2.0.1   2021-02-10 [2] CRAN (R 4.0.5)
    ##  utf8          1.2.1   2021-03-12 [2] CRAN (R 4.0.5)
    ##  vctrs         0.3.8   2021-04-29 [2] CRAN (R 4.0.5)
    ##  withr         2.4.2   2021-04-18 [2] CRAN (R 4.0.5)
    ##  xfun          0.24    2021-06-15 [2] CRAN (R 4.0.3)
    ##  yaml          2.2.1   2020-02-01 [2] CRAN (R 4.0.3)
    ## 
    ## [1] C:/Users/BWaweru/OneDrive - CGIAR/Documents/R/win-library/4.0
    ## [2] C:/R/R-4.0.3/library
