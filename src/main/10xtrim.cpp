#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "10xtrim.hpp"
#include "../overlapper.hpp"
#include "../common.hpp"
#include "../cigar.hpp"


#define VERSION "beta"
#define PROGRAM "10xtrim"

using namespace std;

namespace opt
{
    static unsigned int verbose;
    static string bamfile = "";
    static string output_prefix = "default";
    static int min_score = 20;
    static int padding = 0;
}

void parse_args ( int argc, char *argv[])
{
    // getopt
    extern char *optarg;
    extern int optind, optopt;
    const char* const short_opts = "hvb:p:m:o:";
    const option long_opts[] = {
        {"verbose",             no_argument,        NULL,   'v'},
        {"version",             no_argument,        NULL,   OPT_VERSION},
        {"bam",                 required_argument,  NULL,   'b'},
        {"out",                 required_argument,  NULL,   'o'},
        {"min-score",           required_argument,  NULL,   'm'},
        {"padding",             required_argument,  NULL,   'p'},
        {"help",                no_argument,        NULL,   'h'},
        { NULL, 0, NULL, 0 }
    };


    static const char* VERSION_MESSAGE =
    PROGRAM " version " VERSION "\n"
    "Written by Joanna Pineda.\n"
    "\n";

    static const char* USAGE_MESSAGE =
    "usage: 10xtrim [OPTIONS] --bam test.bam \n"
    "Trim 10x-artifacts and create new bam file\n"
    "\n"
    "    -v, --verbose              Display verbose output\n"
    "        --version              Display version\n"
    "    -b, --bam                  BAM file containing alignment information\n"
    "    -o, --out                  Output prefix for new BAM file\n"
    "        --min-score            Minimum overlap score\n"
    "    -p, --padding              Number of bases added to overlap start when trimming\n";

    int pflag=0; int mflag=0;
    int c;
    while ( (c = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1 ) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch(c) {
        case 'v':
            opt::verbose = 1; // set verbose flag
            break;
        case OPT_VERSION:
            std::cout << VERSION_MESSAGE << endl;
            exit(0);
        case 'b':
            if ( opt::bamfile != "" ) {
                fprintf(stderr, VERSION_MESSAGE);
                fprintf(stderr, "10xtrim: multiple instances of option -b,--bam. \n\n");
                fprintf(stderr, USAGE_MESSAGE, argv[0]);
                exit(EXIT_FAILURE);
            }
            arg >> opt::bamfile;
            //std::cout << "BAM file = " + opt::bamfile << endl;
            break;
       case 'p':
            if ( pflag == 1 ) {
                fprintf(stderr, PROGRAM ": multiple instances of option -p,--padding. \n\n");
                fprintf(stderr, USAGE_MESSAGE, argv[0]);
                exit(EXIT_FAILURE);
            }
            pflag = 1;
            arg >> opt::padding;
       case 'm':
            if ( mflag == 1 ) {
                fprintf(stderr, PROGRAM ": multiple instances of option -m,--min-score. \n\n");
                fprintf(stderr, USAGE_MESSAGE, argv[0]);
                exit(EXIT_FAILURE);
            }
            mflag = 1;
            arg >> opt::min_score;
       case 'o':
            arg >> opt::output_prefix;
            break;
        }
   }

    // check mandatory variables and assign defaults
    if ( opt::bamfile == "" ) {
        fprintf(stderr, "10xtrim: missing -b,--bam option\n\n");
        fprintf(stderr, USAGE_MESSAGE, argv[0]);
        exit(EXIT_FAILURE);
    }

    if (optind < argc) {
        printf("WARNING: invalid option thus ignored: ");
        while (optind < argc)
            printf("%s ", argv[optind++]);
        printf("\n");
    }

}

void trim() {
    // open input bamfile
    const char* in_bamfilename = opt::bamfile.c_str();
    htsFile *infile = hts_open(in_bamfilename,"rb"); //open bam file
    if(infile==NULL) {
        fprintf(stderr, "10xtrim: could not open input bamfile\n\n");
        exit(EXIT_FAILURE);       
    }
    bam1_t *read = bam_init1();
    bam_hdr_t *header = sam_hdr_read(infile);
    // initialize output bamfile
    string temp =  opt::output_prefix + ".bam";
    const char* out_filename = temp.c_str();
    htsFile *outfile= hts_open(out_filename, "wb");
    int ret_val = sam_hdr_write(outfile, header); // copy header
    if ( ret_val < 0 ) {
        fprintf(stderr, "10xtrim: error copying header from input bamfile to new bamfile\n\n");
        exit(EXIT_FAILURE);
    }

    cout << "read_name\tread_length\ttotal_bases_removed\told_cigar\tnew_cigar\toverlap_score\toverlap.match[0].start\toverlap.match[1].start\thairpin_beginning\tseq\n";
    while(sam_read1(infile, header, read) >= 0) {
        if( (read->core.flag & BAM_FUNMAP) == 0 ) {
            // get basic information
            string rname = bam_get_qname(read);
            const bam1_core_t *c = &read->core;

            // get the sequence and its complement
            int qlen = read->core.l_qseq;
            string seq = "";
            string seq_rc = "";
            char nuc;
            //cout << bam_get_seq(read) << "\n";
            //exit(0);
            size_t n = (size_t) read->core.l_qseq;
            for (size_t i = 0; i < n; ++i) {
                nuc =  get_nuc(bam_seqi(bam_get_seq(read), i));
                seq += nuc;
                seq_rc += get_complement_base(nuc); 
            }

            // convert complement to reverse complement
            reverse(seq_rc.begin(), seq_rc.end());
          
            // calculate the overlaps between rev comp and original seq
            SequenceOverlap overlap = Overlapper::computeOverlap(seq, seq_rc);
            //if ( rname == "A00469:17:HFK7YDSXX:1:1422:11071:34162" ) { cout << overlap.score << "\n"; overlap.printAlignment(seq, seq_rc); exit(0); }
            //cout << rname << "\t" << seq <<"\t" << overlap.score << "\t" <<  overlap.match[0].start << "\t" << overlap.match[1].start <<"\n";

            // hairpin detected if
            // - overlap score is greater than thresh
            // - overlap length is greater than 10
            // start modifying cigar
            bool hairpin_detected = ( c->n_cigar && ((int)overlap.score > opt::min_score) && (overlap.total_columns > 10) );
            // old cigar information
            uint32_t *cigar = bam_get_cigar(read);
            vector<char> expanded_cigar = bam_expand_cigar(cigar, (unsigned int)c->n_cigar);
            string old_cigar_str = "";
            uint32_t old_cigar_idx = 0; // position in cigar
            while ( old_cigar_idx < c->n_cigar ) {
                //int type = bam_cigar_op(cigar[old_cigar_idx]);
                int len = bam_cigar_oplen(cigar[old_cigar_idx]);
                char op = BAM_CIGAR_STR[bam_cigar_op(cigar[old_cigar_idx])];
                old_cigar_str += to_string(len) + op;
                old_cigar_idx++;
            }

            bool is_unmapped;
            if ( hairpin_detected && overlap.match[0].start <  5 &&  (overlap.match[1].end > qlen - 5) ) {
                // CASE 1 : HAIRPIN AT THE BEGINNING
                // original:       ====-----------
                // rev comp: ------====
                int clip_pos = overlap.match[0].end + opt::padding;
                if ( clip_pos >= qlen ) clip_pos = overlap.match[0].end;
                int total_bases_to_trim = clip_pos + 1;
                int cigar_pos = 0;
                int read_pos = 0;
                int ref_bases_consumed = 0;
                while ( read_pos < clip_pos ) {
                    char op = expanded_cigar[cigar_pos];
                    if ( op == 'M' || op == 'S' || op == 'I' || op == '=' || op == 'X' ) {
                        expanded_cigar[cigar_pos] = 'S';
                        read_pos++;
                    }
                    if ( op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X' ) ref_bases_consumed++;
                    cigar_pos++;
                }
                vector<uint32_t> new_cigar = compress_expanded_cigar(expanded_cigar);
                // check if unmapped
                if ( new_cigar.size() == 1 ) {
                    is_unmapped = true;
                }
                //cout << clip_pos << "\t" << old_cigar_str  << "\t"<< cigar_ops_to_string(new_cigar) << "\t" << ref_bases_consumed << "\n";   
                cout << rname << "\t" << to_string(qlen) << "\t" << total_bases_to_trim << "\t" << old_cigar_str << "\t" << cigar_ops_to_string(new_cigar) << "\t" << overlap.score << "\t" << overlap.match[0].start << "\t" << overlap.match[1].start << "\t" << "1" << "\t" <<  seq  <<"\n";
                write_new_alignment(header, read, outfile, new_cigar, is_unmapped, ref_bases_consumed);
                // CASE 2 : HAIRPIN AT THE END
                // original: ------====
                // rev comp:       ====-----------
            } else if ( hairpin_detected && overlap.match[0].end < 5 && (overlap.match[1].end > qlen - 5) ) {
                int clip_pos = overlap.match[0].start - opt::padding;
                if ( clip_pos >= 0 ) { 
                    clip_pos = 0;
                    is_unmapped = true;
                }
                int total_bases_to_trim = qlen - clip_pos + 1;
                int cigar_pos = expanded_cigar.size();
                int read_pos = qlen - 1;
                while ( read_pos >= clip_pos ) {
                    char op = expanded_cigar[cigar_pos];
                    if ( op == 'M' || op == 'S' || op == 'I' || op == '=' || op == 'X' ) {
                        expanded_cigar[cigar_pos] = 'S';
                        read_pos--;
                    }
                    cigar_pos--;
                }
                //print_expanded_cigar(expanded_cigar);
                vector<uint32_t> new_cigar = compress_expanded_cigar(expanded_cigar);
                // check if unmapped
                if ( new_cigar.size() == 1 ) {
                    is_unmapped = true;
                }
                cout << rname << "\t" << to_string(qlen) << "\t" <<  total_bases_to_trim << "\t"<< old_cigar_str << "\t" << cigar_ops_to_string(new_cigar) << "\t" << overlap.score << "\t" << overlap.match[0].start << "\t" << overlap.match[1].start << "\t" << "0" << "\t" <<  seq  <<"\n";
                write_new_alignment(header, read, outfile, new_cigar, is_unmapped, 0);
                    //cout << clip_pos << "\t" << old_cigar_str  << "\t"<< cigar_ops_to_string(new_cigar) << "\n";
            } else {
                int ret = sam_write1(outfile, header, read);
                if(ret < 0) {
                    fprintf(stderr, "10xtrim: error writing sam record\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    bam_hdr_destroy(header);
    bam_destroy1(read);
    hts_close(infile);
    hts_close(outfile);
}


void write_new_alignment(bam_hdr_t *header, bam1_t *read, htsFile *outfile, vector<uint32_t> new_cigar, bool is_unmapped, int num_ref_based_ops_trimmed) {
    // make copy of alignment record
    bam1_t *new_record = bam_init1();

    // Variable-length data
    string qname = bam_get_qname(read);

    // basic stats
    new_record->core.tid = read->core.tid;
    new_record->core.pos = read->core.pos + num_ref_based_ops_trimmed;
    new_record->core.qual = read->core.qual;
    new_record->core.l_qname = qname.length() + 1; // must be null-terminated

    new_record->core.flag = read->core.flag;
    if (is_unmapped) {
       new_record->core.flag = (read->core.flag |= 0x0004);
       new_record->core.pos = 0;
    }
    new_record->core.l_qseq = read->core.l_qseq;
    // this is in bases
    new_record->core.l_qseq = read->core.l_qseq;
    uint32_t stored_qseq_len = (new_record->core.l_qseq + 1) / 2;

    new_record->core.mtid = read->core.mtid;
    new_record->core.mpos = read->core.mpos;
    new_record->core.isize = read->core.isize;

    vector<uint32_t> cigar1 = new_cigar;
    new_record->core.n_cigar = cigar1.size();

    // calculate length of incoming data
    new_record->m_data = new_record->core.l_qname + // query name
                          new_record->core.n_cigar * 4 + // 4 bytes per cigar op
                          new_record->core.l_qseq + // query seq
                          new_record->core.l_qseq + // query quality
                          bam_get_l_aux(read); // aux data

    // nothing copied yet
    new_record->l_data = 0;

    // allocate data
    new_record->data = (uint8_t*)malloc(new_record->m_data);

    // copy q name
    assert(new_record->core.l_qname <= new_record->m_data);
    memcpy(bam_get_qname(new_record), bam_get_qname(read), new_record->core.l_qname);
    new_record->l_data += new_record->core.l_qname;

    // cigar
    assert(new_record->l_data + new_record->core.n_cigar * 4 <= new_record->m_data);
    memcpy(bam_get_cigar(new_record), &new_cigar[0], new_record->core.n_cigar * 4);
    new_record->l_data += new_record->core.n_cigar * 4;

    // set bin
    bam1_core_t *c = &new_record->core;
    if (c->n_cigar > 0) { // recompute "bin" and check CIGAR-qlen consistency
        hts_pos_t rlen =  bam_cigar2rlen(c->n_cigar, bam_get_cigar(new_record));
        hts_pos_t qlen = bam_cigar2qlen(c->n_cigar, bam_get_cigar(new_record));
        if ((new_record->core.flag & BAM_FUNMAP)) rlen=1;
        new_record->core.bin = hts_reg2bin(new_record->core.pos, new_record->core.pos + rlen, 14, 5);
        // Sanity check for broken CIGAR alignments
        if (c->l_qseq > 0 && !(c->flag & BAM_FUNMAP) && qlen != c->l_qseq) {
            hts_log_error("CIGAR and query sequence lengths differ for %s",
                    bam_get_qname(new_record));
            exit(1);
        }
    }

    // sequence
    memcpy(bam_get_seq(new_record), bam_get_seq(read), stored_qseq_len );
    new_record->l_data += stored_qseq_len;

    // qualities
    memcpy(bam_get_qual(new_record), bam_get_qual(read), new_record->core.l_qseq );
    new_record->l_data += new_record->core.l_qseq;

    //cout << bam_get_aux(read) << "\n";

    // aux
    memcpy(bam_get_aux(new_record), bam_get_aux(read), bam_get_l_aux(read));
    new_record->l_data += bam_get_l_aux(read);

    int ret = sam_write1(outfile, header, new_record);

    if(ret < 0) {
        fprintf(stderr, "10xtrim: error writing sam record\n");
        exit(EXIT_FAILURE);
    }

    bam_destroy1(new_record); // automatically frees malloc'd segment
}


// main() is where program execution begins.
int main(int argc, char *argv[]) {
   parse_args(argc, argv);
   trim();
   return 0;
}


