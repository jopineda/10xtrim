#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <htslib/sam.h>
#include <htslib/hts.h>

using namespace std;

vector<char> bam_expand_cigar(uint32_t *cigar, uint32_t n_cigar) {
    vector<char> expanded_cigar;
    uint32_t i = 0; // position in cigar
    while ( i < n_cigar ) {
        char op = BAM_CIGAR_STR[bam_cigar_op(cigar[i])];
        int len = bam_cigar_oplen(cigar[i]);
        //string op = bam_cigar_table[]
        for ( int c = 0; c < len; c++ ) { 
            expanded_cigar.push_back(op);
        }
        i++;
    }
	return expanded_cigar;
}

void print_expanded_cigar(vector<char> expanded_cigar) {
	for (auto const& i: expanded_cigar) {
		cout << i;
	}
	cout <<"\n";
}

uint32_t string_op_to_cigar(int len, char op) {
    uint32_t incoming;
    incoming = len << BAM_CIGAR_SHIFT;
    if (op == 'M') incoming |= BAM_CMATCH;
    else if (op == 'I') incoming |= BAM_CINS;
    else if (op == '=') incoming |= BAM_CEQUAL;
    else if (op == 'X') incoming |= BAM_CDIFF;
    else if (op == 'D') incoming |= BAM_CDEL;
    else if (op == 'H') incoming |= BAM_CHARD_CLIP;
    else if (op == 'N') incoming |= BAM_CREF_SKIP;
    else if (op == 'S') incoming |= BAM_CSOFT_CLIP;
    return incoming;
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

vector<uint32_t> compress_expanded_cigar(vector<char> expanded_cigar) {
	int len = 0; 
	char op;
    vector<uint32_t> new_cigar;
    uint32_t incoming;
    char prev_op = expanded_cigar[0];
    for ( const char op : expanded_cigar ) {
    	if ( op != prev_op ) {
    	    incoming = string_op_to_cigar(len, prev_op);
    	    len = 0;
    	    new_cigar.push_back(incoming);
    	    prev_op = op;
    	}
    	len += 1;
    }
    incoming = string_op_to_cigar(len, prev_op);
    new_cigar.push_back(incoming);
    return new_cigar;
}
