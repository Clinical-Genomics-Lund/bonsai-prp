## [Unreleased]

### Added

 - Pytest for Mycobacterium tuberculosis

### Fixed

 - Reordered variant models in ElementTypeResult
 - Fixed assignment of TbProfiler variant type

### Changed

 - Renamed mutations to variants
 - Removed prc_dup and dup from postAlnQc output
 - Added reference, note fields to PhenotypeInfo
 - Mykrobe output parser handles csv format instead of json
 - Split generic ResistanceVariant and ResistanceGene to one specified for individual tools
 - Include non resistance yeilding variants from tbprofiler

## [0.2.0]

### Added

 - CLI command for generating QC report for CDM
 - Parsing of Virulencefinder STX typing

### Fixed

 - Bug that prevented parsing of AMRfinder results

### Changed

 - Resfinder variant describe changed NT instead of codon.
 - Antibiotic class is stored with phenotypes.
 - Accession nr, gene symbol and AA change is stored for Resfinder variants

## [0.1.0]

### Added

 - Added PRP application from https://github.com/genomic-medicine-sweden/JASEN
 - CHANGELOG reminder github action for PRs
 - Test and test files for S.aureus and E.coli
 - Pytest github action

### Fixed

### Changed
