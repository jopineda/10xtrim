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
    for (auto const& i: expanded_cigar) {
		cout << i << " ";
	}
}