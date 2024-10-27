> The Python scripts in this repository were run using Python v.3.9.

> To run the Python script variant_search.py, the following libraries need to be installed:

Built-in Python libraries:

os - https://docs.python.org/3/library/os.html
sys - https://docs.python.org/3/library/sys.html
random - https://docs.python.org/3/library/random.html
datetime - https://docs.python.org/3/library/datetime.html
from collections import defaultdict, Counter - https://docs.python.org/3/library/collections.html

External libraries:

pysam v.0.16.0.1 - https://pysam.readthedocs.io/en/latest/
References:
The Sequence Alignment/Map format and SAMtools. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data       Processing Subgroup. Bioinformatics. 2009 Aug 15;25(16):2078-9. Epub 2009 Jun 8 btp352. PMID: 19505943.
HTSlib: C library for reading/writing high-throughput sequencing data. Bonfield JK, Marshall J, Danecek P, Li H, Ohan V, Whitwham A, Keane T, Davies RM. GigaScience (2021) 10(2) giab007. PMID: 33594436.
Twelve years of SAMtools and BCFtools. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. GigaScience (2021) 10(2) giab008. PMID: 33590861.

numpy v.1.23.4 - https://numpy.org/
    Reference: Oliphant, T. (n.d.). Guide to NumPy: 2nd Edition. CreateSpace. 2006

matplotlib.pyplot and matplotlib.colors v.3.4.2 - https://matplotlib.org/stable/
    Reference: J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.


> To run the Python script ibdseq.py, the following built-in library is required:

bisect - https://docs.python.org/3/library/bisect.html


> The bash script generate_ibds.sh was run using the bash version 5.0.17(1). The following tools are required to use this script:

VCFtools v.0.1.17 - Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R., 1000 Genomes Project Analysis Group, 2011. The variant call format and VCFtools. Bioinformatics 27, 2156–2158. https://doi.org/10.1093/bioinformatics/btr330
Java v.17.0.10 - Arnold, K., Gosling, J. & Holmes, D., 2005. The Java programming language, Addison Wesley Professional.
Beagle v.5.4 - Browning, B.L., Tian, X., Zhou, Y., Browning, S.R., 2021. Fast two-stage phasing of large-scale sequence data. The American Journal of Human Genetics 108, 1880–1890. https://doi.org/10.1016/j.ajhg.2021.08.005
RaPID v.1.7 - Naseri, A., Liu, X., Tang, K., Zhang, S., Zhi, D., 2019. RaPID: ultra-fast, powerful, and accurate detection of segments identical by descent (IBD) in biobank-scale cohorts. Genome Biology 20, 143. https://doi.org/10.1186/s13059-019-1754-8
