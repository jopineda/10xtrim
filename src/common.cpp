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
