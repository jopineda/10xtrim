#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

using namespace std;

vector<char> bam_expand_cigar(uint32_t *cigar, uint32_t n_cigar);