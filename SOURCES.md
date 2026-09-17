# Sources

## NCBI GEO dataset

GSE58331: Gene expression in human orbit  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58331

GEO reports:
- 175 total samples in the full series
- GPL570 Affymetrix Human Genome U133 Plus 2.0 platform
- processed data in each sample table
- `VALUE` as log2 RMA signal
- multiple disease classes, including TED and normal controls
- repeat / technical samples for some subjects

Example public GEO sample pages used in `public_data_example.csv`:
- GSM1407197, TED Anterior Orbit 15
- GSM1407223, Normal Anterior Orbit 6

## Amgen portfolio context

TEPEZZA (teprotumumab) is an Amgen rare-disease medicine for thyroid eye disease.

Amgen 2026 Q1 results:
https://investors.amgen.com/news-releases/news-release-details/amgen-reports-first-quarter-2026-financial-results

## Published biomarker / transcriptomic analyses

Machine learning-based prediction of diagnostic markers for Graves' orbitopathy:
https://pubmed.ncbi.nlm.nih.gov/37059863/

TEDML:
https://pubmed.ncbi.nlm.nih.gov/40028837/

A more recent transcriptomic re-analysis of orbital tissue also used GSE58331 and curated TED / normal orbital samples:
https://pmc.ncbi.nlm.nih.gov/articles/PMC12783674/
