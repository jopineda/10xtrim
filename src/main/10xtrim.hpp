#include <stdio.h>
#include <stdlib.h>
using namespace std;
void parse_args(int argc, char *argv[]);
int calculateScore(string seq, string seq_rc, int match_score, int penalty_score);
string cigar_ops_to_string(const vector<uint32_t>& ops);	
void write_new_alignment(bam_hdr_t *header, bam1_t *read, htsFile *outfile, vector<uint32_t> new_cigar, bool is_unmapped, int num_ref_based_ops_trimmed);
enum { OPT_VERSION };
