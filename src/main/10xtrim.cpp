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

    int bflag=0; int pflag=0; int mflag=0;
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
            if ( bflag == 1 ) {
                fprintf(stderr, VERSION_MESSAGE);
                fprintf(stderr, "10xtrim: multiple instances of option -b,--bam. \n\n");
                fprintf(stderr, USAGE_MESSAGE, argv[0]);
                exit(EXIT_FAILURE);
            }
            bflag = 1;
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
    if ( bflag == 0 ) {
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

const char* cigar_op_to_string(int type) {
    switch(type) {
        case BAM_CMATCH:  return "M";
        case BAM_CINS: return "I";
        case BAM_CDEL: return "D";
        case BAM_CREF_SKIP: return "N";
        case BAM_CSOFT_CLIP: return "S";
        case BAM_CHARD_CLIP: return "H";
        case BAM_CPAD: return "P";
        case BAM_CEQUAL: return "=";
        case BAM_CDIFF: return "X";
    }
}

uint32_t string_op_to_cigar(int len, string op) {
    uint32_t incoming;
    incoming = len << BAM_CIGAR_SHIFT;
    if (op == "M") incoming |= BAM_CMATCH;
    else if (op == "I") incoming |= BAM_CINS;
    else if (op == "=") incoming |= BAM_CEQUAL;
    else if (op == "X") incoming |= BAM_CDIFF;
    else if (op == "D") incoming |= BAM_CDEL;
    else if (op == "H") incoming |= BAM_CHARD_CLIP;
    else if (op == "N") incoming |= BAM_CREF_SKIP;
    else if (op == "S") incoming |= BAM_CSOFT_CLIP;
    return incoming;
}

void trim() {
    // open input bamfile
    const char* in_bamfilename = opt::bamfile.c_str();
    htsFile *infile = hts_open(in_bamfilename,"rb"); //open bam file
    bam1_t *read = bam_init1();
    bam_hdr_t *header = sam_hdr_read(infile);
    // initialize output bamfile
    string temp =  opt::output_prefix + ".bam";
    const char* out_filename = temp.c_str();
    htsFile *outfile= hts_open(out_filename, "wb");
    int ret_val = sam_hdr_write(outfile, header); // copy header

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
            if ( c->n_cigar && ((int)overlap.score > opt::min_score) && (overlap.total_columns > 10) ) {
                //if ( rname == "A00469:17:HFK7YDSXX:1:1422:11071:34162" ) { cout <<seq <<"\n"; cout << seq_rc <<"\n"; cout << overlap.score << "\n"; overlap.printAlignment(seq, seq_rc); exit(0); }

                // old cigar information
                uint32_t *cigar = bam_get_cigar(read);
                vector<char> expanded_cigar = bam_expand_cigar(cigar, (unsigned int)c->n_cigar);
                string old_cigar_str = "";
                uint32_t old_cigar_idx = 0; // position in cigar
                while ( old_cigar_idx < c->n_cigar ) {
                    int type = bam_cigar_op(cigar[old_cigar_idx]);
                    int len = bam_cigar_oplen(cigar[old_cigar_idx]);
                    string op = cigar_op_to_string(type);
                    old_cigar_str += to_string(len) + op;
                    old_cigar_idx++;
                }

                bool is_unmapped;
                // CASE 1 : HAIRPIN AT THE BEGINNING
                // original:       ====-----------
                // rev comp: ------====
                if ( overlap.match[0].start <  5 &&  (overlap.match[1].end > qlen - 5) ) {
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
                    //cout << clip_pos << "\t" << old_cigar_str  << "\t"<< cigar_ops_to_string(new_cigar) << "\t" << ref_bases_consumed << "\n";   
                    cout << rname << "\t" << to_string(qlen) << "\t" << total_bases_to_trim << "\t" << old_cigar_str << "\t" << cigar_ops_to_string(new_cigar) << "\t" << overlap.score << "\t" << overlap.match[0].start << "\t" << overlap.match[1].start << "\t" << "1" << "\t" <<  seq  <<"\n";
                    write_new_alignment(header, read, outfile, new_cigar, is_unmapped, ref_bases_consumed);
                // CASE 2 : HAIRPIN AT THE END
                // original: ------====
                // rev comp:       ====-----------
                } else if ( overlap.match[0].end < 5 && (overlap.match[1].end > qlen - 5) ) {
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
                    cout << rname << "\t" << to_string(qlen) << "\t" <<  total_bases_to_trim << "\t"<< old_cigar_str << "\t" << cigar_ops_to_string(new_cigar) << "\t" << overlap.score << "\t" << overlap.match[0].start << "\t" << overlap.match[1].start << "\t" << "0" << "\t" <<  seq  <<"\n";
                    write_new_alignment(header, read, outfile, new_cigar, is_unmapped, 0);
                    //cout << clip_pos << "\t" << old_cigar_str  << "\t"<< cigar_ops_to_string(new_cigar) << "\n";
                }
                



                // get old cigar for debugging
                /*string old_cigar_str = "";
                uint32_t old_cigar_idx = 0; // position in cigar
                while ( old_cigar_idx < c->n_cigar ) {
                    int type = bam_cigar_op(cigar[old_cigar_idx]);
                    int len = bam_cigar_oplen(cigar[old_cigar_idx]);
                    string op = cigar_op_to_string(type);
                    old_cigar_str += to_string(len) + op;
                    old_cigar_idx++;
                }

                // initialize new cigar
                string new_cigar_str = "";
                int type = 0; int len = 0; string op = "";
                vector<uint32_t> new_cigar;
                uint32_t incoming;

                // initialize variables necessary for trimming
                int clip_pos = 0;               // position to clip read from
                bool is_unmapped = false;       // marks if entire read gets softclipped
                int total_bases_trimmed = 0;    // total bases trimmed
                int num_bases_to_remove = 0;    // number of bases remaining to remove
                int num_ref_based_ops_trimmed = 0; // may need to update left-most map pos
                uint32_t cigar_idx = 0;         // cigar operation index
 
                // CASE 1 : HAIRPIN AT THE BEGINNING
                // original:       ====-----------
                // rev comp: ------====
                if ( overlap.match[0].start < 5 && overlap.match[1].end > qlen - 5) {
                    // find where to clip and add padding if requested
                    clip_pos = overlap.match[0].end + opt::padding;
                    if ( clip_pos >= qlen ) clip_pos = overlap.match[0].end;
                    assert(clip_pos < qlen);
                    num_bases_to_remove = clip_pos + 1;
                    total_bases_trimmed = num_bases_to_remove;

                    // first cigar operation
                    type = bam_cigar_op(cigar[0]);
                    len = bam_cigar_oplen(cigar[0]);
                    op = cigar_op_to_string(type); 

                    // add softclip to the beginning
                    int init = 0;
                    int cigar_len = 0;
                    // if already softclipped, and current softclip len does not include
                    // the bases we want to trim, we need to add to softclip
                    if ( op == "S" && len < clip_pos + 1 ) {
                        num_bases_to_remove -= len;
                        incoming = string_op_to_cigar(clip_pos + 1, "S");
                        cigar_len += clip_pos + 1;
                        init++;
                    // if already softclipped, and current softclip len already includes
                    // the bases we want to trim, we don't do anything
                    } else if ( op == "S" && len >= clip_pos + 1) {
                        num_bases_to_remove = 0;
                        cigar_len += len;
                        incoming = string_op_to_cigar(len, "S");
                    // if we don't start with a softclip, then add softclip to beginning
                    } else {
                        cigar_len += clip_pos + 1;
                        incoming = string_op_to_cigar(clip_pos + 1, "S"); 
                    }
                    new_cigar.push_back(incoming);
 
                    // remove bases from front until removed num_bases_to_remove
                    while ( num_bases_to_remove > 0 && init < c->n_cigar) {
                        // collect cigar operation info
                        type = bam_cigar_op(cigar[init]);
                        len = bam_cigar_oplen(cigar[init]);
                        op = cigar_op_to_string(type);

                        // if we still have to remove bases and the next operation is softclipped
                        // the entire read is being softclipped, we need to mark as unmapped
                        if (op == "S") {
                            init = c->n_cigar; // init = next position in cigar to start adding 
                            new_cigar.clear();
                            incoming = string_op_to_cigar(qlen, "S");
                            new_cigar.push_back(incoming); 
                            is_unmapped = true;
                        }
                        // otherwise move on to next read-operation
                        if (op == "M" || op == "I" || op == "=" || op == "X" ) {
                            int new_len = len - num_bases_to_remove;
                            // stop softclipping
                            if ( len > num_bases_to_remove ) { // remove part of the cigar operation only
                                incoming = string_op_to_cigar(new_len, op);
                                new_cigar.push_back(incoming);
                                cigar_len += new_len;
                                if ( op == "M" || op == "=" || op == "X" ) {
                                    num_ref_based_ops_trimmed += num_bases_to_remove;
                                }
                                break;
                            } else if ( len == num_bases_to_remove ) {
                                num_bases_to_remove -= len;
                                break;
                            } else { // remove cigar operation completely
                                num_bases_to_remove -= len;
                            }
                        } else {
                            // reference-based operation, just add to new cigar
                            incoming = string_op_to_cigar(len, op);
                            new_cigar.push_back(incoming);
                            num_ref_based_ops_trimmed += len;
                        }
                        init++;
                    }
                    
                    if (init <= c->n_cigar && init != 0) {
                        type = bam_cigar_op(cigar[init]);
                        len = bam_cigar_oplen(cigar[init]);
                        op = cigar_op_to_string(type);

                        // if we still have to remove bases and the next operation is softclipped
                        // the entire read is being softclipped, we need to mark as unmapped
                        if (op == "S") { 
                            is_unmapped = true;
                            init = c->n_cigar; 
                        }
                    }
                    // add the remaining cigar string with modified beginning
                    for (cigar_idx = init + 1; cigar_idx < c->n_cigar; cigar_idx++) {
                        type = bam_cigar_op(cigar[cigar_idx]);
                        len = bam_cigar_oplen(cigar[cigar_idx]);
                        op = cigar_op_to_string(type);
                        incoming = string_op_to_cigar(len, op);
                        new_cigar.push_back(incoming);
                        //cout << op << to_string(len)<< "\n";
                        if (op == "M" || op == "I" || op == "=" || op == "X" || op =="S" ) {
                            cigar_len += len;
                        }
                    }
                    if ( is_unmapped ) {
                        new_cigar_str = to_string(qlen) + "S";
                        cigar_len = qlen;
                        new_cigar.clear();
                        incoming = string_op_to_cigar(qlen, "S");
                        new_cigar.push_back(incoming);
                    }
                    cout << rname << "\t" << to_string(qlen) << "\t" << to_string(total_bases_trimmed) << "\t" << old_cigar_str << "\t" << cigar_ops_to_string(new_cigar) << "\t" << overlap.score << "\t" << overlap.match[0].start << "\t" << overlap.match[1].start << "\t" << "1" << "\t" <<  seq  <<"\n";
                    assert(cigar_len == qlen);
                    write_new_alignment(header, read, outfile, new_cigar, is_unmapped, num_ref_based_ops_trimmed);
                // CASE 2 : HAIRPIN AT THE END
                // original: ------====
                // rev comp:       ====-----------
                } else {
                    // determine position to start trimming from
                    clip_pos = overlap.match[0].start - opt::padding;
                    if ( clip_pos < 0 ) {
                        clip_pos = 0; // will clip everything
                        is_unmapped = true;
                    }
                    assert( clip_pos >= 0 );
                    total_bases_trimmed = qlen - clip_pos - 1;

                    int pos = 0; int cigar_len = 0; int len = 0;
                    cigar_idx = type = 0; op = "";
                    // add everything up until clip position
                    while ( cigar_idx < c->n_cigar ) {
                        type = bam_cigar_op(cigar[cigar_idx]);
                        len = bam_cigar_oplen(cigar[cigar_idx]);
                        op = cigar_op_to_string(type);
                        if (op == "M" || op == "I" || op == "=" || op == "X" || op == "S") {
                            pos += len;
                        }
                        if ( pos > clip_pos + 1 ) { // reached position to stop clipping
                            break; 
                        } else if ((pos == clip_pos + 1) && (op == "M" || op == "I" || op == "=" || op == "X" || op == "S")){
                            incoming = string_op_to_cigar(len, op);
                            new_cigar.push_back(incoming);
                            cigar_len += len;
                            //cigar_idx++;
                            break;
                        }
                        // keep adding to the cigar strings
                        incoming = string_op_to_cigar(len, op);
                        new_cigar.push_back(incoming);
                        if (op == "M" || op == "I" || op == "=" || op == "X" || op == "S" ) {
                            cigar_len += len;
                        }
                        cigar_idx++;
                    }


                    // Start trimming depending on case
                    if ( op != "S" ) {
                        // CASE 1: Last operation needs to be trimmed
                        // Ex. CIGAR = 10S15M, clip pos = 15, NEW CIGAR = 10S5M10S
                        if ( pos  > clip_pos + 1 ) {
                            int remove_len = pos - clip_pos;
                            int new_len = len - remove_len + 1;
                            cigar_len += new_len;
                            incoming = string_op_to_cigar(new_len, op);
                            new_cigar.push_back(incoming);
                        }

                        // CASE 2: Last operation needs to be completely trimmed
                        // Ex. CIGAR = 10S5I15M, clip pos = 15, NEW CIGAR = 10S5I15S 
                        int softclip = qlen - clip_pos - 1;
                        incoming = string_op_to_cigar(softclip, "S");
                        new_cigar.push_back(incoming);
                        cigar_len += softclip;
                    } else {
                        // CASE 3: Last operation is softclipped and includes trimmed bases
                        // Ex. CIGAR = 10S15M10S, clip pos = 30, NEW CIGAR = 10S15M10S
                        incoming = string_op_to_cigar(len, op);
                        new_cigar.push_back(incoming);
                        cigar_len += len;
                        // CASE 4: Trimming from already softclipped position at the front
                        // Ex. CIGAR = 85S15M, clip pos = 10, NEW CIGAR = 100S
                        if ( cigar_idx == 0 ) {
                            is_unmapped = true;
                        } 
                    }
                    if ( is_unmapped ) {
                        new_cigar.clear();
                        incoming = string_op_to_cigar(qlen, "S");
                        new_cigar.push_back(incoming);
                        cigar_len = qlen;
                    }
                                
                    cout << rname << "\t" << to_string(qlen) << "\t" << to_string(total_bases_trimmed) << "\t" << old_cigar_str << "\t" << cigar_ops_to_string(new_cigar) << "\t" << overlap.score << "\t" << overlap.match[0].start << "\t" << overlap.match[1].start << "\t" << "0" << "\t" <<  seq  <<"\n";
                    assert(cigar_len == qlen);
                    write_new_alignment(header, read, outfile, new_cigar, is_unmapped, 0);
                }*/
            } else {
                int ret = sam_write1(outfile, header, read);
                if(ret < 0) {
                    fprintf(stderr, "error writing sam record\n");
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
        fprintf(stderr, "error writing sam record\n");
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


