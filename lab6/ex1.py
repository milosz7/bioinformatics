from Bio import Blast, Entrez, SeqIO
from Bio.Blast import NCBIXML
import argparse


def load_sequence(path):
    with open(path, "r") as file:
        return file.read().strip()


def fetch_entrez_sequence(seq_id, db):
    with Entrez.efetch(db=db, id=seq_id, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
    return record.seq


def query_blast(algorithm, db, sequence, n_results=1):
    result_handle = Blast.qblast(algorithm, db, sequence, hitlist_size=n_results)
    return result_handle


def parse_blast_results(result_handle):
    parsed_results = NCBIXML.parse(result_handle)
    return parsed_results


def print_alignments(blast_records):
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            print('Alignment----------------------------------')
            print('title:', alignment.title)
            print('length:', alignment.length)
            for hsp in alignment.hsps:
                print('HSP : ')
                print('e value:', hsp.expect)
                print(hsp.query[0:75] + '...')
                print(hsp.match[0:75] + '...')
                print(hsp.sbjct[0:75] + '...')


def run_demo():
    Entrez.email = input("Email (for Entrez): ")
    sequence_1 = """
    GTACCTTGATTTCGTATTCTGAGAGGCTGCTGCTTAGCGGTAGCCCCTTGGTTTCCGTGGCAACGGAAAA
    GCGCGGGAATTACAGATAAATTAAAACTGCGACTGCGCGGCGTGAGCTCGCTGAGACTTCCTGGACGGGG
    GACAGGCTGTGGGGTTTCTCAGATAACTGGGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTA
    AAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAAT
    GCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTG
    ACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTT
    ATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTA
    TTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAA
    AGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCG
    TGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTGCAGGAAACCAGTCTCAGTGTCCAACTC
    TCTAACCTTGGAACTGTGAGAACTCTGAGGACAAAGCAGCGGATACAACCTCAAAAGACGTCTGTCTACA
    TTGAATTGGGATCTGATTCTTCTGAAGATACCGTTAATAAGGCAACTTATTGCAGTGTGGGAGATCAAGA
    ATTGTTACAAATCACCCCTCAAGGAACCAGGGATGAAATCAGTTTGGATTCTGCAAAAAAGGCTGCTTGT
    GAATTTTCTGAGACGGATGTAACAAATACTGAACATCATCAACCCAGTAATAATGATTTGAACACCACTG
    AGAAGCGTGCAGCTGAGAGGCATCCAGAAAAGTATCAGGGTAGTTCTGTTTCAAACTTGCATGTGGAGCC
    ATGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACTAAAGACAGA
    ATGAATGTAGAAAAGGCTGAATTCTGTAATAAAAGCAAACAGCCTGGCTTAGCAAGGAGCCAACATAACA
    GATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACTCCCAGCACAGAAAAAAAGGTAGATCTGAA
    TGCTGATCCCCTGTGTGAGAGAAAAGAATGGAATAAGCAGAAACTGCCATGCTCAGAGAATCCTAGAGAT
    ACTGAAGATGTTCCTTGGATAACACTAAATAGCAGCATTCAGAAAGTTAATGAGTGGTTTTCCAGAAGTG
    ATGAACTGTTAGGTTCTGATGACTCACATGATGGGGAGTCTGAATCAAATGCCAAAGTAGCTGATGTATT
    GGACGTTCTAAATGAGGTAGATGAATATTCTGGTTCTTCAGAGAAAATAGACTTACTGGCCAGTGATCCT
    CATGAGGCTTTAATATGTAAAAGTGAAAGAGTTCACTCCAAATCAGTAGAGAGTAATATTGAAGACAAAA
    TATTTGGGAAAACCTATCGGAAGAAGGCAAGCCTCCCCAACTTAAGCCATGTAACTGAAAATCTAATTAT
    AGGAGCATTTGTTACTGAGCCACAGATAATACAAGAGCGTCCCCTC
    """
    result_handle = query_blast("blastn", "nt", sequence_1)
    results_1 = parse_blast_results(result_handle)
    print_alignments(results_1)

    sequence_2 = """
    MKSILDGLADTTFRTITTDLLGSPFQEKMTAGDNPQLVPADQVNITEFYNKSLSSFKENEENIQCGENFM
    DIECFMVLNPSQQLAIAVLSLTLGTFTVLENLLVLCVILHSRSLRCRPSYHFIGSLAVADLLGSVIFVYS
    FIDFHVFHRKDSRNVFLFKLGGVTASFTASVGSLFLTAIDRYISIHRPLAYKRIVTRPKAVVAFCLMWTI
    AIVIAVLPLLGWNCEKLQSVCSDIFPHIDETYLMFWIGVTSVLLLFIVYAYMYILWKAHSHAVRMIQRGT
    QKSIIIHTSEDGKVQVTRPDQARMDIRLAKTLVLILVVLIICWGPLLAIMVYDVFGKMNKLIKTVFAFCS
    MLCLLNSTVNPIIYALRSKDLRHAFRSMFPSCEGTAQPLDNSMGDSDCLHKHANNAASVHRAAESCIKST
    VKIAKVTMSVSTDTSAEAL 
    """
    result_handle = query_blast("blastp", "nr", sequence_2)
    results_2 = parse_blast_results(result_handle)
    print_alignments(results_2)

    sequence_3_id = "NM_000539.3"
    sequence_3 = fetch_entrez_sequence(sequence_3_id, "nucleotide")
    result_handle = query_blast("tblastx", "refseq_select_rna", sequence_3, n_results=5)
    results_3 = parse_blast_results(result_handle)
    print_alignments(results_3)


def main(args):
    if not args.demo:
        if args.file:
            sequence = load_sequence(args.file)
        elif args.seq_id:
            Entrez.email = input("Email (for Entrez): ")
            sequence = fetch_entrez_sequence(args.seq_id, args.seq_db)
        else:
            raise ValueError(
                "Either --file or --seq_id must be provided unless --demo is used. Use --help for more information.")

        result_handle = query_blast(args.algorithm, args.db, sequence, args.n_results)
        blast_record = parse_blast_results(result_handle)
        print_alignments(blast_record)
    else:
        run_demo()

if __name__ == "__main__":
    # Use -h to see usage information
    parser = argparse.ArgumentParser(description='BLAST query alignment parser')
    parser.add_argument('--algorithm', type=str, help='BLAST algorithm to use (default - blastn)', default='blastn')
    parser.add_argument('--db', type=str, help='BLAST database to use (default - nt)', default='nt')
    parser.add_argument("--demo", action='store_true', help='Run demo with predefined sequences')
    parser.add_argument("--n_results", type=int, help='Number of BLAST results to retrieve (default - 1)', default=1)
    parser.add_argument('--file', type=str, help='Path to the file containing the sequence to query')
    parser.add_argument('--seq_id', type=str, help='ID of the sequence to query')
    parser.add_argument('--seq_db', type=str, help='Database for fetching sequence by ID (default - nucleotide)',
                        default='nucleotide')
    args = parser.parse_args()
    main(args)
