# UTR_Checker
Assess for U3-R-U5 presence in user-provided HIV sequences.

Problem: HIV sequences are often plagued with sequence defects and/or sequencing/assembly artifacts.

Aim: develop a tool to screen for HIV-1 LTR sequences, including U3 region, R region, U5 region at the 5' and 3' ends. 

Background: HIV DNA has LTR-HIV-LTR structure, mosre specifically U3-R-U5-HIV-U3-R-U5. 5' U3-R-U5 should be identical to 3' U3-R-U5. HIV sequences with different sequences denote possible sequencing artifacts or synthetic constructs. E.g., NL4-3 made from two HIV-1 cDNAs spliced together during molecular cloning, as opposed to during an infective cycle. Note, the current historical HIV reference genome GenBank:K03455.1 has sequencing artifacts in its U3, R, and U5:

```
CLUSTAL format alignment by MAFFT (v7.511)


HXB2_5'_LTR_U3_region     tggaagggctaattcactcccaacgaagacaagatatccttgatctgtggatctaccaca
HXB2_3'_LTR_U3_region     tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccaca
                          *********************** ************************************

HXB2_5'_LTR_U3_region     cacaaggctacttccctgattagcagaactacacaccagggccagggatcagatatccac
HXB2_3'_LTR_U3_region     cacaaggctacttccctgattagcagaactacacaccagggccaggggtcagatatccac
                          ***********************************************.************

HXB2_5'_LTR_U3_region     tgacctttggatggtgctacaagctagtaccagttgagccagagaagttagaagaagcca
HXB2_3'_LTR_U3_region     tgacctttggatggtgctacaagctagtaccagttgagccagataagatagaagaggcca
                          ******************************************* *** *******.****

HXB2_5'_LTR_U3_region     acaaaggagagaacaccagcttgttacaccctgtgagcctgcatggaatggatgacccgg
HXB2_3'_LTR_U3_region     ataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccgg
                          *.********************************************.*************

HXB2_5'_LTR_U3_region     agagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacatggcccgag
HXB2_3'_LTR_U3_region     agagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacgtggcccgag
                          **************************************************.*********

HXB2_5'_LTR_U3_region     agctgcatccggagtacttcaagaactgctgacatcgagcttgctacaagggactttccg
HXB2_3'_LTR_U3_region     agctgcatccggagtacttcaagaactgctgacatcgagcttgctacaagggactttccg
                          ************************************************************

HXB2_5'_LTR_U3_region     ctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagat
HXB2_3'_LTR_U3_region     ctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagat
                          ************************************************************

HXB2_5'_LTR_U3_region     cctgcatataagcagctgctttttgcctgtactgg
HXB2_3'_LTR_U3_region     cctgcatataagcagctgctttttgcctgtactgg
                          ***********************************
```
```
CLUSTAL format alignment by MAFFT (v7.511)


HXB2_3'_LTR_R_repeat      gtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccact
HXB2_5'_LTR_R_repeat      gtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccact
                          ************************************************************

HXB2_3'_LTR_R_repeat      gcttaagcctcaataaagcttgccttgagtgcttca
HXB2_5'_LTR_R_repeat      gcttaagcctcaataaagcttgccttgagtgcttca
                          ************************************
```

```
CLUSTAL format alignment by MAFFT (v7.511)


HXB2_3'_LTR_U5_region     agtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagaccctttta
HXB2_5'_LTR_U5_region     agtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagaccctttta
                          ************************************************************

HXB2_3'_LTR_U5_region     gtcagtgtggaaaatctctagca
HXB2_5'_LTR_U5_region     gtcagtgtggaaaatctctagca
                          ***********************
```
```
CLUSTAL format alignment by MAFFT (v7.511)


MW079479.1:1-634          tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccaca
MW079479.1:9086-9719      tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccaca
                          ************************************************************

MW079479.1:1-634          cacaaggctacttccctgattagcagaactacacaccagggccaggggtcagatatccac
MW079479.1:9086-9719      cacaaggctacttccctgattagcagaactacacaccagggccaggggtcagatatccac
                          ************************************************************

MW079479.1:1-634          tgacctttggatggtgctacaagctagtaccagttgagccagataaggtagaagaggcca
MW079479.1:9086-9719      tgacctttggatggtgctacaagctagtaccagttgagccagataaggtagaagaggcca
                          ************************************************************

MW079479.1:1-634          ataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccgg
MW079479.1:9086-9719      ataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccgg
                          ************************************************************

MW079479.1:1-634          agagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacgtggcccgag
MW079479.1:9086-9719      agagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacgtggcccgag
                          ************************************************************

MW079479.1:1-634          agctgcatccggagtacttcaagaactgctgatatcgagcttgctacaagggactttccg
MW079479.1:9086-9719      agctgcatccggagtacttcaagaactgctgatatcgagcttgctacaagggactttccg
                          ************************************************************

MW079479.1:1-634          ctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagat
MW079479.1:9086-9719      ctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagat
                          ************************************************************

MW079479.1:1-634          cctgcatataagcagctgctttttgcctgtactgggtctctctggttagaccagatctga
MW079479.1:9086-9719      cctgcatataagcagctgctttttgcctgtactgggtctctctggttagaccagatctga
                          ************************************************************

MW079479.1:1-634          gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
MW079479.1:9086-9719      gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
                          ************************************************************

MW079479.1:1-634          tgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
MW079479.1:9086-9719      tgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
                          ************************************************************

MW079479.1:1-634          agacccttttagtcagtgtggaaaatctctagca
MW079479.1:9086-9719      agacccttttagtcagtgtggaaaatctctagca
                          **********************************
```
[MW079479.1](https://www.ncbi.nlm.nih.gov/nuccore/MW079479) is first described in version 1 of Gener, 2019[^1].
[^1]: Alejandro R. Gener. "Full-coverage sequencing of HIV-1 provirus from a reference plasmid" bioRxiv 611848; doi: https://doi.org/10.1101/611848.

The most current work is Gener et al. 2021[^2].

[^2]: Alejandro R. Gener, Wei Zou, Brian T. Foley, Deborah P. Hyink, Paul E. Klotman. "Reference plasmid pHXB2_D is an HIV-1 molecular clone that exhibits identical LTRs and a single integration site indicative of an HIV provirus" bioRxiv 611848; doi: https://doi.org/10.1101/611848.

(Pairwise alignment done with MAFFT online server. Method FFT-NS-i (Standard). Command: mafft --reorder --auto input.)

References:
   Katoh et al. (2002) describes FFT-NS-1, FFT-NS-2 and FFT-NS-i.
   
   Kuraku et al. (2013) outlines this web service.
   
MAFFT home:
   https://mafft.cbrc.jp/alignment/software/
