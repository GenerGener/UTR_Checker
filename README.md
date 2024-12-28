# UTR_Checker
Assess for U3-R-U5 presence in user-provided HIV sequences.

# UTR_Checker development notes
Problem: HIV sequences are often plagued with sequence defects and/or sequencing/assembly artifacts.

Aim: To develop a tool to screen for HIV-1 LTR sequences, including U3 region, R region, U5 region at the 5' and 3' ends. 

Background: HIV DNA has LTR-HIV-LTR structure, more specifically U3-R-U5-HIV-U3-R-U5. 5' U3-R-U5 should be identical to 3' U3-R-U5. HIV sequences with different sequences denote possible sequencing artifacts or synthetic constructs. E.g., NL4-3 made from two HIV-1 cDNAs spliced together during molecular cloning, as opposed to during an infective cycle. Note, the current historical HIV reference genome GenBank:K03455.1 has sequencing artifacts in its U3 and U5:

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
These sequencing artifacts do not align with the reported identity of HXB2, which was believed to be a captured host-integrated HIV-1 provirus from early in the prepandemic HIV-1 epidemic era. Based on the current model of retroviral replication, LTRs homogenize after RT and prior to integration. More recent resequencing of a legacy HXB2 molecular clone pHXB2 D (pHXB2_D as Genbank:[MW079479.1](https://www.ncbi.nlm.nih.gov/nuccore/MW079479)) enabled correction of this incongruity. Because it was made with high-depth long-read sequencing, it is the most accurate representation of HXB2. As such, landmarks were used for the development of UTR_Checker "utr-checker.py". This tool aims to determine how much of U3, R, U5 are present at 3 and 5' ends of user-provided HIV-1 sequence. Based on this, the tool may then classify input as HIV DNA, RNA, or partial/incomplete.

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

MW079479.1:1-634          cctgcatataagcagctgctttttgcctgtactgg[gtctctctggttagaccagatctga
MW079479.1:9086-9719      cctgcatataagcagctgctttttgcctgtactgg[gtctctctggttagaccagatctga
                          ***********************************[*************************

MW079479.1:1-634          gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
MW079479.1:9086-9719      gcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgcct
                          ************************************************************

MW079479.1:1-634          tgagtgcttca]agtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
MW079479.1:9086-9719      tgagtgcttca]agtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
                          ***********]*************************************************

MW079479.1:1-634          agacccttttagtcagtgtggaaaatctctagca
MW079479.1:9086-9719      agacccttttagtcagtgtggaaaatctctagca
                          **********************************

[R]
```
HXB2 LTR heterogeneity discovered after pHXB2-D resequencing was first described in version 1 of Gener, 2019[^1].
[^1]: Alejandro R. Gener. "Full-coverage sequencing of HIV-1 provirus from a reference plasmid" bioRxiv 611848; doi: https://doi.org/10.1101/611848.

HXB2 as pHXB2_D plasmid assembly and data release as [MW079479.1](https://www.ncbi.nlm.nih.gov/nuccore/MW079479) was described in the most current work: Gener et al. 2021[^2].

[^2]: Alejandro R. Gener, Wei Zou, Brian T. Foley, Deborah P. Hyink, Paul E. Klotman. "Reference plasmid pHXB2_D is an HIV-1 molecular clone that exhibits identical LTRs and a single integration site indicative of an HIV provirus" bioRxiv 611848; doi: https://doi.org/10.1101/611848.

(Pairwise alignment done with MAFFT online server. Method FFT-NS-i (Standard). Command: mafft --reorder --auto input.)

References:
   Katoh et al. (2002) describes FFT-NS-1, FFT-NS-2 and FFT-NS-i.
   
   Kuraku et al. (2013) outlines this web service.
   
MAFFT home:
   https://mafft.cbrc.jp/alignment/software/

# UTR-Checker Tutorial

## Overview
UTR-Checker is a Python script designed to analyze HIV sequences for the presence and arrangement of U3, R, and U5 regions. It uses minimap2 for initial candidate identification followed by detailed alignment analysis, making it efficient for both short and long sequences. Version 13 "utr-checker-13.py" works as expected on HXB2 DNA (K03455) as ACGT and NL4-3 mRNA as ACGT MZ242719.

## Dependencies

### Required Python Packages
```bash
pip install biopython mappy
```

The script requires:
- Python 3.6 or higher
- Biopython (for sequence handling and alignment)
- mappy (Python bindings for minimap2)
- Built-in libraries: tempfile, os, logging, typing, argparse

## Installation

1. Save the script as `utr-checker.py`
2. Ensure your input sequences are in FASTA format
3. Make the script executable (Unix/Linux):
   ```bash
   chmod +x utr-checker.py
   ```

## Usage

### Basic Usage
```bash
python utr-checker.py input.fasta
```

The script accepts any file containing FASTA-formatted sequences, regardless of extension (.fasta, .fsa, .fa, .txt).

### Advanced Options
```bash
python utr-checker.py --minimap-threshold 0.55 --final-threshold 0.65 --gap-open -2 --gap-extend -0.5 --debug input.fasta
```

Parameters:
- `--minimap-threshold`: Initial screening threshold (default: 0.60)
- `--final-threshold`: Final similarity threshold (default: 0.70)
- `--gap-open`: Gap opening penalty (default: -2)
- `--gap-extend`: Gap extension penalty (default: -0.5)
- `--debug`: Enable debug output
- `--format`: Input file format (default: fasta)

## Example Analysis

Using sample test files:

### Example 1: HIV-1 mRNA (MZ242719.1)
```bash
user@computer:~$ python utr-checker.py MZ242719.fasta --minimap-threshold 0.55 --final-threshold 0.65

Analyzing sequence: MZ242719.1
Best match found on forward strand
Classification: Likely viral RNA
Overall confidence: 99.45%

Details:
- Multiple R regions detected (5' and 3' ends)
- U5 region present near 5' end
- U3 region present near 3' end
- U3 occurrences:
-   1. 98.35% similarity at position 8622-9077
- R occurrences:
-   1. 100.00% similarity at position 9077-9173
-   2. 97.92% similarity at position 2-95
- U5 occurrences:
-   1. 100.00% similarity at position 98-181
```

### Example 2: HIV-1 DNA (HXB2)
```bash
user@computer:~$ python utr-checker.py HIV-1_HXB2.fasta --minimap-threshold 0.55 --final-threshold 0.65

Analyzing sequence: K03455.1
Best match found on forward strand
Classification: Incomplete/Unclear
Overall confidence: 99.78%

Details:
- Partial or unclear LTR pattern
- U3 occurrences:
-   1. 99.34% similarity at position 9085-9540
-   2. 97.03% similarity at position 0-455
- R occurrences:
-   1. 100.00% similarity at position 455-551
-   2. 100.00% similarity at position 9540-9636
- U5 occurrences:
-   1. 100.00% similarity at position 551-634
-   2. 100.00% similarity at position 9636-9719
```

## Sequence Classification

The script classifies sequences into:

1. "Likely viral RNA"
   - R regions at both ends
   - U5 near 5' end
   - U3 near 3' end
   - Expected pattern: R-U5-genome-U3-R

2. "Likely genomic DNA"
   - Complete U3-R-U5 pattern
   - Found in correct order
   - May be present at both ends

3. "Incomplete/Unclear"
   - Regions present but in unexpected arrangement
   - Missing expected regions
   - Ambiguous pattern

## Performance Considerations

1. Two-Step Analysis
   - Initial fast screening using minimap2
   - Detailed alignment for candidate regions
   - Adjustable thresholds for both steps

2. Speed vs Accuracy
   - Lower thresholds increase sensitivity but may add false positives
   - Higher thresholds increase specificity but might miss divergent sequences
   - Default values optimized for HIV-1 group M

3. Memory Usage
   - Efficient with long sequences due to minimap2
   - Memory scales with sequence length
   - Temporary files used for minimap2 analysis

## Best Practices

1. Threshold Selection
   - Start with default thresholds
   - Lower minimap-threshold for divergent sequences
   - Adjust final-threshold based on expected similarity

2. Result Interpretation
   - Check both similarity scores and positions
   - Verify region order matches expected pattern
   - Consider biological context (RNA vs DNA)

3. Troubleshooting
   - Use --debug flag for detailed output
   - Check for sequence quality issues
   - Verify FASTA format is correct

## Limitations

1. Reference Sequences
   - Based on HIV-1 HXB2 references
   - May have reduced sensitivity for highly divergent strains
   - Best suited for HIV-1 group M subtype B analysis

2. Structure Detection
   - Analyzes full sequence for all elements (U3, R, U5)
   - Detects both terminal and internal matches
   - Uses position information for pattern classification (e.g., R-U5-genome-U3-R for RNA)

3. Format Requirements
   - Input can be any file containing FASTA-formatted sequences
   - Common extensions (.fasta, .fsa, .fa, .txt) all supported
   - Can process multiple sequences in a single file
   - Handles both DNA and RNA sequences


