## SNV pipeline for SMaHT consortium

### Currently in development

The current scripts in this directory were all originally authored by Jiny Ha and
are for postprocessing of SNV calls for deep sequenced Illumina and PacBio
data.

Specifically common SNPs from dbSNP and gnomad are removed as mSNV candidates. Then
PacBio long reads are used to remove and mSNVs that are common between many samples (PoE
approach). Next long reads are phased to determine if expected haplotypes are observed
at mSNV locations.

More details to follow...
