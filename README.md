# UTR_Checker
Assess for U3-R-U5 presence in user-provided HIV sequences.

Problem: HIV sequences are often plagued with sequence defects and/or sequencing/assembly artifacts.

Aim: develop a tool to screen for HIV-1 LTR sequences, including U3 region, R region, U5 region at the 5' and 3' ends. 

Background: HIV DNA has LTR-HIV-LTR structure, mosre specifically U3-R-U5-HIV-U3-R-U5. 5' U3-R-U5 should be identical to 3' U3-R-U5. HIV sequences with different sequences denote possible sequencing artifacts or synthetic constructs. E.g., NL4-3 made from two HIV-1 cDNAs spliced together during molecular cloning, as opposed to during an infective cycle. Note, the current historical HIV reference genome GenBank:K03455.1 has sequencing artifacts in its U3, R, and U5:

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
