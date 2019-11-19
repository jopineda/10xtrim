# 10xtrim

Program for trimming 10x-specific artifacts.

## Dependencies

* A compiler that supports C++11
* htslib

## Installation Instructions
You can download and compile the latest code from github as follows:

```
git clone --recursive https://github.com/jopineda/10xtrim.git
cd 10xtrim
make
```

## Workflow Examples

### Data preprocessing

```
# mark duplicates
# index bam
```

### Generate new BAM with additional softclip
```
./10xtrim -b input.bam -o output > 10xtrim.stats.tsv
```
