# STO-MS-plus
Quantitative profiling of stoichiometry with 15nt mass tags and LFQ

## STO-MS workflow
A chemoproteomic strategy called “STO-MS” to systematically quantify the PTM stoichiometry in complex biological samples
[github](https://github.com/wangchulab/STO-MS) [paper](https://pubs.rsc.org/en/content/articlelanding/2022/cb/d2cb00179a)

## Upgrade instructions

1. A novel nucleic acid mass tag with better resolution.
2. Label-free quantificaton method allows for simpler sample preparation and higher sample throughput.
3. Gaussian fitting of band intensity allows for better identification of modification species.

## Example:
```python
cd example
python ../scripts/maxq_compare_proGrp_LFQ.py \ 
uniprot_sprot_HUMAN_160504.fasta \
exp.list \
markers.txt \
. \
DMSO_ \
EXP_
```
Options explanation:
```
    fasta with all protein sequnces
    id of each fraction
    shift of each marker
    working directory
    prefix of the first sample
    prefix of the second sample
    ...
```
