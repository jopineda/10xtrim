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

#define VERSION "beta"
#define SUBPROGRAM "10xtrim"

using namespace std;

namespace opt
{
    static unsigned int verbose;
    static string bamfile = "";
    static string output_prefix = "default";
}

void parse_args ( int argc, char *argv[])
{
    // getopt
    extern char *optarg;
    extern int optind, optopt;
    const char* const short_opts = "hvb:";
    const option long_opts[] = {
        {"verbose",             no_argument,        NULL,   'v'},
        {"version",             no_argument,        NULL,   OPT_VERSION},
        {"bam",                 required_argument,  NULL,   'b'},
        {"out",                 required_argument,  NULL,   'o'}, 
        {"help",                no_argument,        NULL,   'h'},
        { NULL, 0, NULL, 0 }
    };


    static const char* VERSION_MESSAGE =
    "10xtrim " SUBPROGRAM " version " VERSION "\n"
    "Written by Joanna Pineda.\n"
    "\n";

    static const char* USAGE_MESSAGE =
    "usage: 10xtrim [OPTIONS] --bam test.bam \n"
    "Create new bam file\n"
    "\n"
    "    -v, --verbose              Display verbose output\n"
    "        --version              Display version\n"
    "    -b, --bam                  BAM file containing alignment information\n"
    "    -o, --out                  Output prefix\n";

    int bflag=0;
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
                fprintf(stderr, "10xtrim: multiple instances of option -b,--bam. \n\n");
                fprintf(stderr, USAGE_MESSAGE, argv[0]);
                exit(EXIT_FAILURE);
            }
            bflag = 1;
            arg >> opt::bamfile;
            //std::cout << "BAM file = " + opt::bamfile << endl;
            break;
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

char get_nuc(int x) {
    switch (x) {
        case 1: return 'A';
        case 2: return 'C';
        case 4: return 'G';
        case 8: return 'T';
        case 15: return 'N';
        default: return '.';
    }
}

char get_complement_base(char x) {
    switch (x) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default: return '.';
    }
}

const char* get_string_op(int type) {
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
//code from: https://github.com/ernfrid/diagnose_dups/blob/17f63ed3d07c63f9b55dc7431f6528707d30709f/src/lib/common/Utility.cpp
/*vector<uint32_t> parse_string_to_cigar_vector(char const* cigar_string) {
        char const* beg = cigar_string;
        std::size_t len = strlen(beg);
        char const* end = beg + len;

        std::vector<uint32_t> cigar;
        // 2x performance increase
        cigar.reserve(len);

        for(; *beg != '\0'; ++beg) {
            uint32_t oplen = 0;
            // qi numeric parsing is much faster than strtoul
            // it takes beg by reference and moves it to the first non-numeric
            // character
            if (UNLIKELY(!auto_parse(beg, end, oplen) || !valid_cigar_len(oplen))) {
                // do you want to get mad if this isn't a number?
                throw std::runtime_error(str(format(
                    "Error parsing cigar string %1%: expected number at position %2%"
                    ) % cigar_string % (beg - cigar_string)));
            }

            int op = opcode_for_char(*beg);

            if (UNLIKELY(op < 0)) {
                throw std::runtime_error(str(format(
                    "Error parsing cigar string %1%: invalid cigar op char at position %2%"
                    ) % cigar_string % (beg - cigar_string)));
            }

            cigar.push_back(bam_cigar_gen(oplen, op));
        }
        return cigar;
}*/

uint32_t get_cigar_op(int len, string op) {
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

void trim(string file) {
    const char* in_bamfilename = file.c_str();
    htsFile *infile = hts_open(in_bamfilename,"rb"); //open bam file
    bam1_t *read = bam_init1();
    bam_hdr_t *header = sam_hdr_read(infile);
    string temp = "output.sam";
    const char* out_filename = temp.c_str();
    htsFile *outfile= hts_open(out_filename, "w");
    int ret_val = sam_hdr_write(outfile, header); // copy header
    int i, r = 0;

    while(sam_read1(infile, header, read) >= 0) {
        if( (read->core.flag & BAM_FUNMAP) == 0 ) {
            // get the sequence
            // translated seq A=>0,C=>1,G=>2,T=>3,other=>4
            string seq = "";
            string seq_rc = "";
            char nuc;
            size_t n = (size_t) read->core.l_qseq;
            for (size_t i = 0; i < n; ++i) {
                nuc =  get_nuc(bam_seqi(bam_get_seq(read), i));
                seq += nuc;
                seq_rc += get_complement_base(nuc); 
            }
            // get the reverese complement of the sequence
            reverse(seq_rc.begin(), seq_rc.end());
          
            // calculate the overlaps between rev comp and original seq
            // auto score = calculateScore(seq, seq_rc, 2, -1);           
            SequenceOverlap overlap = Overlapper::computeOverlap(seq, seq_rc);
            int score = overlap.score;
            int qlen = read->core.l_qseq;

            // TESTING:
            string target_name = "A00469:17:HFK7YDSXX:2:1138:22245:2754";
            string rname = bam_get_qname(read);
            if (rname == target_name) {
                overlap.printAlignment(seq, seq_rc);
                cout << rname << "\t" << seq <<"\t" << overlap.score << "\t" <<  overlap.match[0].start << "\t" << overlap.match[1].start <<"\n";
            }

            // update cigar if overlap score is greater than thresh
            const bam1_core_t *c = &read->core;
            if ( c->n_cigar && overlap.score > 20 ) { 
                // old cigar information
                uint32_t *cigar = bam_get_cigar(read);
                int clip_pos = 0; // position to clip read from

                // initialize new cigar
                string new_cigar_str = "";
                string op = "";
                vector<uint32_t> new_cigar;
                uint32_t incoming;

                // get old cigar for debugging
                string old_cigar_str = "";
                uint32_t cigar_idx = 0; // position in cigar
                while ( cigar_idx < c->n_cigar ) {
                    int type = bam_cigar_op(cigar[cigar_idx]);
                    int len = bam_cigar_oplen(cigar[cigar_idx]);
                    string op = get_string_op(type);
                    old_cigar_str += to_string(len) + op;
                    cigar_idx++;
                }

                // alignment at the beginning of original seq = hairpin at the beginning
                cigar_idx = 0;
                // original:       ====-----------
                // rev comp: ------====
                if ( overlap.match[0].start < 5 ) {
                    // find where to clip
                    clip_pos = overlap.match[0].end;
                    int num_bases_to_remove = clip_pos - 1;

                    // first cigar operation
                    int type = bam_cigar_op(cigar[0]);
                    int len = bam_cigar_oplen(cigar[0]);
                    string op = get_string_op(type); 

                    // if first operation is already softclip just add to it
                    int init = 0;
                    // if softclipped len < clip_pos
                    if ( op == "S" && len < clip_pos ) {
                        num_bases_to_remove -= len;
                        // add softclip to beginning
                        new_cigar_str = to_string(clip_pos) + "S";
                        incoming = get_cigar_op(clip_pos, "S");
                        init++;
                    } else if ( op == "S" && len >= clip_pos ) {
                        // if softclip len >= clip_pos
                        // add the length of the softclip
                        new_cigar_str = to_string(len) + "S"; 
                        num_bases_to_remove = 0;
                        incoming = get_cigar_op(len, "S");
                    } else {
                        // add softclip to beginning
                        new_cigar_str = to_string(clip_pos) + "S";
                        incoming = get_cigar_op(len, "S"); 
                    }
                    new_cigar.push_back(incoming); 
                    // remove bases from front until removed num_bases_to_remove
                    while ( num_bases_to_remove > 0 && init < c->n_cigar) {
                        //cout << " num_bases_to_remove= " << num_bases_to_remove << " op= " << op << " len= " << len << " cigar="<<new_cigar_str << "\n"; 
                        type = bam_cigar_op(cigar[init]);
                        len = bam_cigar_oplen(cigar[init]);
                        op = get_string_op(type);
                        // if we still have to remove bases and the next operation is softclipped
                        // the entire read is being softclipped, we need to mark as unmapped
                        if (op == "S") { new_cigar_str = to_string(qlen) + "S"; init = c->n_cigar; break; }
                        // otherwise move on to next read-operation
                        if (op == "M" || op == "I" || op == "=" || op == "X" ) {
                            int new_len = len - num_bases_to_remove - 1;
                            // stop softclipping
                            if ( len > num_bases_to_remove ) { // remove part of the cigar operation only
                                incoming = get_cigar_op(new_len, op);
                                new_cigar.push_back(incoming);
                                new_cigar_str += to_string(new_len) + op;
                                break;
                            } else { // remove cigar operation completely
                                num_bases_to_remove -= len;
                            }
                        } else {
                            // reference-based operation, just add to new cigar
                            new_cigar_str += to_string(len) + op;
                            incoming = get_cigar_op(len, op);
                            new_cigar.push_back(incoming);       
                        }
                        init++;
                    }
                    // add the remaining cigar string with modified beginning
                    int pos = clip_pos;
                    for (cigar_idx = init + 1; cigar_idx < c->n_cigar;cigar_idx++) {
                        type = bam_cigar_op(cigar[cigar_idx]);
                        len = bam_cigar_oplen(cigar[cigar_idx]);
                        op = get_string_op(type);
                        incoming = get_cigar_op(len, op);
                        new_cigar.push_back(incoming);
                        new_cigar_str += to_string(len) + op;
                        //cigar_ops_.push_back(CigarOp(bam_cigar_opchr(cigars[i]), bam_cigar_oplen(cigars[i])));
                    }
                    if (rname == target_name) {
                        cout << rname << " with len of " << to_string(qlen) <<"\n"; 
                        cout << "tried to trim: " << clip_pos  << "\n";
                        cout << "old: " << old_cigar_str  << "\n";
                        cout << "new: " << new_cigar_str << "\n";
                        cout << "new: " << cigar_ops_to_string(new_cigar) << "\n";
                        //exit(0);
                    }
                    write_new_alignment(infile, header, read, outfile, new_cigar);
                    //cout << new_cigar_str << "\n";
                //} else if ( overlap.match[1].start < 5 ) { // hairpin at the end
                } else {
                    // original: ------====
                    // rev comp:       ====-----------
                    int clip_pos = overlap.match[0].start; // clip everything after this position
                    int pos = 0; 
                    int cigar_idx = 0;
                    string new_cigar_str = "";
                    int type = 0; int len = 0; string op = "";
                    //cout << "before while\n";
                    while ( cigar_idx < c->n_cigar ) {
                        //cout << "in while\n";
                        type = bam_cigar_op(cigar[cigar_idx]);
                        len = bam_cigar_oplen(cigar[cigar_idx]);
                        op = get_string_op(type);
                        //cout << "operation " << op << "\n";
                        if (op == "M" || op == "I" || op == "=" || op == "X" || op == "S") {
                            pos += len;
                        }
                        //cout << "pos= " << pos << " clip_pos= " << clip_pos << " op= " << op << " len= " << len << " cigar="<<new_cigar_str << "\n"; 
                        if ( pos >= clip_pos ) { break; }
                        // keep adding to the cigar strings
                        new_cigar_str += to_string(len) + op;
                        incoming = get_cigar_op(len, op);
                        new_cigar.push_back(incoming);
                        cigar_idx++;
                    }
                    //cout <<"done while\n";
                    //cout << "old*: " << old_cigar_str  << "\n";
                    if ( pos > clip_pos ) {
                        int remove_len = pos - clip_pos;
                        int new_len = len - remove_len;
                        incoming = get_cigar_op(new_len, op);
                        new_cigar.push_back(incoming);
                        new_cigar_str += to_string(new_len) + op;
                    }
                    int qlen = bam_cigar2qlen(c->n_cigar, cigar);
                    int softclip = qlen - clip_pos;
                    incoming = get_cigar_op(softclip, "S");
                    new_cigar.push_back(incoming);
                    new_cigar_str += to_string(softclip) + "S";

                                
                    if (rname == "A00469:17:HFK7YDSXX:2:2213:30779:4726") {
                        cout << rname << " with len of " << to_string(qlen) <<"\n"; 
                        cout << "tried to trim: " << softclip  << "\n";
                        cout << "old: " << old_cigar_str  << "\n";
                        cout << "new: " << new_cigar_str << "\n";
                        cout << "new: " << cigar_ops_to_string(new_cigar) << "\n";
                    }
                    //cout << "done\n";
                    write_new_alignment(infile, header, read, outfile, new_cigar);
                }
            } else {
                int ret = sam_write1(outfile, header, read);

                if(ret < 0) {
                    fprintf(stderr, "error writing sam record\n");
                    exit(EXIT_FAILURE);
                }
            }
            r++;
            //cout <<"read number " << r <<"\n";
        }
    }

    bam_hdr_destroy(header);
    bam_destroy1(read);
    hts_close(infile);
    hts_close(outfile);
}

string cigar_ops_to_string(const vector<uint32_t>& ops)
{
    stringstream ss;
    for(size_t i = 0; i < ops.size(); ++i) {
        ss << bam_cigar_oplen(ops[i]);
        ss << BAM_CIGAR_STR[bam_cigar_op(ops[i])];
    }
    return ss.str();
}

void write_new_alignment(htsFile *infile, bam_hdr_t *header, bam1_t *read, htsFile *outfile, vector<uint32_t> new_cigar) {
    //cout << "writing\n";
    // make copy of alignment record
    bam1_t *new_record = bam_init1();

    // Variable-length data
    string qname = bam_get_qname(read);

    // basic stats
    new_record->core.tid = read->core.tid;
    new_record->core.pos = read->core.pos;
    new_record->core.qual = read->core.qual;
    new_record->core.l_qname = qname.length() + 1; // must be null-terminated

    new_record->core.flag = read->core.flag;

    new_record->core.l_qseq = read->core.l_qseq;

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
    strncpy(bam_get_qname(new_record),
           qname.c_str(),
           new_record->core.l_qname);
    new_record->l_data += new_record->core.l_qname;

    // cigar
    assert(new_record->l_data + new_record->core.n_cigar * 4 <= new_record->m_data);
    memcpy(bam_get_cigar(new_record), &new_cigar[0], new_record->core.n_cigar * 4);
    new_record->l_data += new_record->core.n_cigar * 4;

    // sequence
    memcpy(bam_get_seq(new_record), bam_get_seq(read), new_record->core.l_qseq );
    new_record->l_data += new_record->core.l_qseq;

    // qualities
    memcpy(bam_get_qual(new_record), bam_get_qual(read), new_record->core.l_qseq );
    new_record->l_data += new_record->core.l_qseq;

    //cout << bam_get_aux(read) << "\n";

    // aux
    memcpy(bam_get_aux(new_record), bam_get_aux(read), bam_get_l_aux(read));
    //new_record->l_data += bam_get_l_aux(new_record);

    // no copy for seq and qual
    //assert(new_record->l_data <= new_record->m_data);

    int ret = sam_write1(outfile, header, new_record);

    if(ret < 0) {
        fprintf(stderr, "error writing sam record\n");
        exit(EXIT_FAILURE);
    }

    bam_destroy1(new_record); // automatically frees malloc'd segment
}

int calculateScore(string seq, string seq_rc, int match_score, int penalty_score) 
{ 
    int m = seq.length(); // length of gene1 
    int n = seq_rc.length(); // length of gene2 
      
    // table for storing optimal substructure answers 
    int score[n] = {0}; 
  
    // get scores
    int i = 0;
    if ( seq[i] == seq_rc[i] ) { 
        score[i] = match_score;
    } else {
        score[i] = penalty_score;
    }
    
    // get max value in score
    int current_max = score[0];
    int current_max_index = 0;
    for (i = 1; i < n; i++) 
    {
        //cout << seq[i] << " and " << seq_rc[i] << "\n";
        if ( seq[i] == seq_rc[i] ) { 
            // sequences match
            score[i] = score[i - 1] + match_score;
        } else {
            // sequence doesn't match
            score[i] = score[i - 1] + penalty_score;
        }
        if ( score[i] >= current_max ) {
            current_max = score[i];
            current_max_index = i;
        }
        //cout << score[i] << "\n";
    }

    //cout << current_max << "\n";


  
    return current_max; 
} 


// main() is where program execution begins.
int main(int argc, char *argv[]) {
   parse_args(argc, argv);
   //calculateScore("ACCCT", "ACCCT", 1, -2);
   trim(opt::bamfile);

   return 0;
}


