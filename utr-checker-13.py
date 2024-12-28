#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
import mappy as mp
import logging
from typing import Tuple, Dict, List
import tempfile
import os

class UTRChecker:
    def __init__(self, minimap_threshold: float = 0.60, final_threshold: float = 0.70):
        """
        Initialize UTR Checker with configurable thresholds
        
        Args:
            minimap_threshold: Initial screening threshold for minimap2 (more permissive)
            final_threshold: Final threshold for detailed alignment
        """
        self.utr_sequences = {
            'U3': "tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccacacacaaggctacttccctgattagcagaactacacaccagggccaggggtcagatatccactgacctttggatggtgctacaagctagtaccagttgagccagataaggtagaagaggccaataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccggagagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacgtggcccgagagctgcatccggagtacttcaagaactgctgatatcgagcttgctacaagggactttccgctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagatcctgcatataagcagctgctttttgcctgtactgg",
            'R': "gtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttca",
            'U5': "agtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagca"
        }
        
        self.minimap_threshold = minimap_threshold
        self.final_threshold = final_threshold
        self.sequence_lengths = {region: len(seq) for region, seq in self.utr_sequences.items()}
        
        # Parameters for detailed alignment
        self.match_score = 2
        self.mismatch_score = -1
        self.gap_open = -2
        self.gap_extend = -0.5
        
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

    def clean_sequence(self, sequence: str) -> str:
        """Clean input sequence by removing whitespace and converting to uppercase"""
        return ''.join(sequence.split()).upper()

    def find_candidate_regions(self, sequence: str, region_type: str) -> List[Tuple[int, int]]:
        """
        Use minimap2 to find all potential candidate regions
        Returns list of (start, end) positions
        """
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as ref_file, \
             tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
            
            # Write reference sequence
            ref_file.write(f">{region_type}\n{self.utr_sequences[region_type]}\n")
            ref_file.flush()
            
            # Write query sequence
            query_file.write(f">query\n{sequence}\n")
            query_file.flush()
            
            # Initialize minimap2 aligner with parameters for finding multiple matches
            # k=10: smaller k-mer size for higher sensitivity
            # w=1: minimize minimizer window size
            # f=0.0001: very low minimum minimizer frequency threshold
            # max-chain-skip=100: allow more skipping in chaining
            aligner = mp.Aligner(ref_file.name, preset="map-ont", k=10, w=1, 
                               min_chain_score=20, min_dp_score=20, best_n=10)
            if not aligner:
                raise Exception("Failed to initialize minimap2 aligner")
            
            candidates = []
            sequence_length = len(sequence)
            
            # Process all hits from minimap2
            for hit in aligner.map(sequence):
                is_terminal = hit.q_st < 1000 or hit.q_en > sequence_length - 1000
                map_quality = hit.mapq / 60.0  # Normalize mapping quality
                
                if region_type in ['R', 'U5']:
                    # For R and U5, focus on terminal regions but keep high-quality internal matches
                    if is_terminal or map_quality > self.minimap_threshold:
                        candidates.append((hit.q_st, hit.q_en, map_quality))
                else:
                    # For U3, consider all regions meeting quality threshold
                    if map_quality > self.minimap_threshold:
                        candidates.append((hit.q_st, hit.q_en, map_quality))
            
            # Clean up temp files
            os.unlink(ref_file.name)
            os.unlink(query_file.name)
            
            # Sort by position and quality, remove duplicates
            candidates.sort(key=lambda x: (x[0], -x[2]))
            unique_candidates = []
            for start, end, quality in candidates:
                if not any(abs(start - prev_start) < 50 for prev_start, _, _ in unique_candidates):
                    unique_candidates.append((start, end, quality))
            
            return [(start, end) for start, end, _ in unique_candidates]

    def calculate_alignment_score(self, seq1: str, seq2: str) -> Tuple[float, str, str]:
        """Calculate detailed alignment score using pairwise2"""
        alignments = pairwise2.align.localms(
            seq1, seq2,
            self.match_score,
            self.mismatch_score,
            self.gap_open,
            self.gap_extend,
            one_alignment_only=True
        )
        
        if not alignments:
            return 0.0, "", ""
        
        alignment = alignments[0]
        align_seq1, align_seq2 = alignment[0], alignment[1]
        max_possible_score = len(seq1) * self.match_score
        normalized_score = alignment[2] / max_possible_score
        
        return normalized_score, align_seq1, align_seq2

    def analyze_candidate_region(self, reference: str, sequence: str, start: int, end: int) -> Tuple[float, int, int]:
        """Analyze a candidate region with detailed alignment"""
        # Add padding to capture potential alignment extensions
        padding = 30
        window_start = max(0, start - padding)
        window_end = min(len(sequence), end + padding)
        window = sequence[window_start:window_end]
        
        similarity, align_ref, align_window = self.calculate_alignment_score(reference, window)
        
        # Adjust positions to account for padding
        if similarity >= self.final_threshold:
            match_start = start  # Use minimap2's positions as they're often more accurate
            match_end = end
            return similarity, match_start, match_end
        
        return similarity, -1, -1

    def find_region_matches(self, reference: str, sequence: str) -> List[Tuple[float, int, int]]:
        """Find all matches above threshold using minimap2 pre-filtering"""
        sequence = self.clean_sequence(sequence)
        reference = self.clean_sequence(reference)
        matches = []
        
        # Determine region type
        region_type = 'U3' if len(reference) > 200 else ('R' if "AAGCTT" in reference else 'U5')
        
        # Find candidate regions
        candidates = self.find_candidate_regions(sequence, region_type)
        logging.debug(f"Found {len(candidates)} candidate regions for {region_type}")
        
        # Analyze each candidate region
        for start, end in candidates:
            similarity, match_start, match_end = self.analyze_candidate_region(
                reference, sequence, start, end)
            
            if similarity >= self.final_threshold:
                matches.append((similarity, match_start, match_end))
                logging.debug(f"Confirmed match for {region_type}: {match_start}-{match_end} ({similarity:.2%})")
        
        return sorted(matches, key=lambda x: x[0], reverse=True)

    def analyze_sequence(self, sequence: str) -> Dict:
        """Analyze sequence for UTR regions"""
        sequence = self.clean_sequence(sequence)
        results = {
            'regions': {},
            'classification': '',
            'details': []
        }
        
        for region, ref_seq in self.utr_sequences.items():
            matches = self.find_region_matches(ref_seq, sequence)
            
            if matches:
                all_matches = [{'similarity': m[0], 'start': m[1], 'end': m[2]} 
                             for m in matches]
                results['regions'][region] = {
                    'present': True,
                    'matches': all_matches,
                    'expected_length': self.sequence_lengths[region],
                    'similarity': matches[0][0],
                    'start': matches[0][1],
                    'end': matches[0][2],
                    'length': matches[0][2] - matches[0][1]
                }
            else:
                results['regions'][region] = {
                    'present': False,
                    'matches': [],
                    'similarity': 0.0,
                    'start': -1,
                    'end': -1,
                    'length': 0,
                    'expected_length': self.sequence_lengths[region]
                }
        
        results.update(self._classify_sequence(results['regions']))
        return results

    def _find_terminal_r_regions(self, r_matches: List[Dict], u5_matches: List[Dict]) -> Tuple[Dict, Dict]:
        """Find 5' and 3' R regions"""
        five_prime_r = None
        three_prime_r = None
        
        # Find 5' R region (near start)
        for match in r_matches:
            if match['start'] < 200:
                five_prime_r = match
                break
        
        # Find 3' R region (far from start)
        threshold = 5000
        if u5_matches:
            threshold = u5_matches[0]['start'] + 1000
            
        for match in r_matches:
            if match['start'] > threshold:
                three_prime_r = match
                break
                
        return five_prime_r, three_prime_r

    def _classify_sequence(self, regions: Dict) -> Dict:
        """Classify sequence based on detected UTR regions"""
        classification = ''
        details = []
        
        if not any(region['present'] for region in regions.values()):
            return {
                'classification': 'Non-LTR sequence',
                'details': ["No LTR regions detected with sufficient similarity"]
            }
        
        r_matches = regions['R']['matches'] if regions['R']['present'] else []
        u5_matches = regions['U5']['matches'] if regions['U5']['present'] else []
        u3_matches = regions['U3']['matches'] if regions['U3']['present'] else []
        
        # Check for RNA pattern
        if len(r_matches) >= 2:
            five_prime_r, three_prime_r = self._find_terminal_r_regions(r_matches, u5_matches)
            
            if five_prime_r and three_prime_r:
                classification = 'Likely viral RNA'
                details.append("Multiple R regions detected (5' and 3' ends)")
                if u5_matches:
                    details.append("U5 region present near 5' end")
                if u3_matches:
                    details.append("U3 region present near 3' end")
                details.extend(self._generate_region_details(regions))
                return {'classification': classification, 'details': details}
        
        # Check for DNA pattern
        if u3_matches and r_matches and u5_matches:
            u3_pos = u3_matches[0]['start']
            r_pos = r_matches[0]['start']
            u5_pos = u5_matches[0]['start']
            
            if u3_pos < r_pos < u5_pos:
                return {
                    'classification': 'Likely genomic DNA',
                    'details': ["Complete U3-R-U5 regions found in correct order"] + 
                              self._generate_region_details(regions)
                }
        
        return {
            'classification': 'Incomplete/Unclear',
            'details': ["Partial or unclear LTR pattern"] + 
                      self._generate_region_details(regions)
        }

    def _generate_region_details(self, regions: Dict) -> List[str]:
        """Generate detailed analysis of each detected region"""
        details = []
        for region, data in regions.items():
            if data['present']:
                details.append(f"{region} occurrences:")
                for idx, match in enumerate(data['matches'], 1):
                    details.append(f"  {idx}. {match['similarity']:.2%} similarity at position "
                                 f"{match['start']}-{match['end']}")
        return details

    def check_both_strands(self, sequence: str) -> Dict:
        """Check both forward and reverse complement strands"""
        forward_results = self.analyze_sequence(sequence)
        rev_comp = str(Seq(sequence).reverse_complement())
        reverse_results = self.analyze_sequence(rev_comp)
        
        forward_confidence = sum(region['similarity'] for region in forward_results['regions'].values()) / len(self.utr_sequences)
        reverse_confidence = sum(region['similarity'] for region in reverse_results['regions'].values()) / len(self.utr_sequences)
        
        if reverse_confidence > forward_confidence:
            return {
                'strand': 'reverse',
                'results': reverse_results,
                'forward_results': forward_results,
                'overall_confidence': reverse_confidence
            }
            
        return {
            'strand': 'forward',
            'results': forward_results,
            'reverse_results': reverse_results,
            'overall_confidence': forward_confidence
        }

def main():
    parser = argparse.ArgumentParser(description='Check HIV sequences for UTR regions')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('--format', default='fasta', help='Input file format (default: fasta)')
    parser.add_argument('--minimap-threshold', type=float, default=0.60,
                       help='Similarity threshold for minimap2 screening (0.0-1.0)')
    parser.add_argument('--final-threshold', type=float, default=0.70,
                       help='Similarity threshold for final sequence matching (0.0-1.0)')
    parser.add_argument('--gap-open', type=float, default=-2,
                       help='Gap opening penalty (default: -2)')
    parser.add_argument('--gap-extend', type=float, default=-0.5,
                       help='Gap extension penalty (default: -0.5)')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    checker = UTRChecker(
        minimap_threshold=args.minimap_threshold,
        final_threshold=args.final_threshold
    )
    checker.gap_open = args.gap_open
    checker.gap_extend = args.gap_extend
    
    try:
        for record in SeqIO.parse(args.input_file, args.format):
            print(f"\nAnalyzing sequence: {record.id}")
            results = checker.check_both_strands(str(record.seq))
            
            print(f"Best match found on {results['strand']} strand")
            print(f"Classification: {results['results']['classification']}")
            print(f"Overall confidence: {results['overall_confidence']:.2%}")
            print("\nDetails:")
            for detail in results['results']['details']:
                print(f"- {detail}")
            
            print("\nRegion detection details:")
            for region, data in results['results']['regions'].items():
                print(f"\n{region}:")
                if data['present']:
                    print("  Present: True")
                    for idx, match in enumerate(data['matches'], 1):
                        print(f"  Match {idx}:")
                        print(f"    Similarity: {match['similarity']:.2%}")
                        print(f"    Position: {match['start']}-{match['end']}")
                else:
                    print("  Present: False")
                print(f"  Expected length: {data['expected_length']} bp")
            
            print("-" * 50)
            
    except Exception as e:
        logging.error(f"Error processing file: {str(e)}")
        raise

if __name__ == "__main__":
    main()