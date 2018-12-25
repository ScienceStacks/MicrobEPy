# Sequence Data
## Mutation Data
## \*allsamples_attributes_date.xlsx
Below is listed the columns and their descriptions.
- variant_id - genome
- experiment - when/how the sample was obtained
  - EPD_seq - endpoint dilution sequencing
  - after_300g - after 300 generations, not 1000 generations
  - clonal-isolate - details in the sample column
  - 1000-gen - 1K generation
  - Ancestors - Ancestor sequences
    - AN_00_S7_Stahl - DvH ancestor from Stahl lab sequencing
    - AN_Coculture-Ancestor - Co-culture ancestor from Hillesland lab sequencing
    - AN_Dv-Ancestor-1 - DvH ancestor from Hillesland lab sequencing (DvH (sensitive to naladixic acid))
    - AN_Dv-Ancestor-2 - DvH ancestor from Hillesland lab sequencing (DvH (abx-resistant))
    - AN_Mm-Ancestor-1 - Mmp ancestor from Hillesland lab sequencing ( Mmp (sensitive to neomycin))
    - AN_Mm-Ancestor-2 - Mmp ancestor from Hillesland lab sequencing (Mmp (abx-resistant))
- sample - provides line, sample
  - EG: Early generations, EP: EPDs, CI: Clonal Isolates, AN: Ancestor, TG: 1000-gen
  - for EPD_seq and after_300g: line (e.g., UA3) and transfer number (e.g., 10)
  - for clonal-isolates: _NN_ (e.g., _28_) is the pairing number that
    links to an isolate using the files names_clonal_isolates_*.csv
- source - chromosome or plasmid
- position - position on the genome source
- initial_nt_1st - 1st base in the original nucleotides sequence of the reference genome
- initial_nt - original nucleotide(s) of the reference genome
- changed_nt - nucleotides after the mutation
- mutation - characterization of the mutation (e.g., stop codon added)
- effect (see:  http://snpeff.sourceforge.net/SnpEff_manual.html#input for reference)
  - MODIFIER - intergenic
  - LOW - e.g., synonymous mutation
  - MODERATE - e.g., same kind of amino acid (e.g., hydrophobic)
  - HIGH - e.g., possible knockout (e.g., stop gained)
- type - characterization of the mutation
- codon_change - old codon/new codon
- aa_change - change in the Amino acid
- gene_id - gene idenfier
- region - CODING or INTERGENIC
- freq - IGNORE
- predictor - list of tools used in the mutation call
- freq_predictor - mutation occurrence frequency for the tool
  in the same colon-separated position in the predictor column
- read_number - number of reads captured from the variant call that support the variant
## Updates 03022018
### \*allsamples_attributes_03022018
- Created an .xlsx file
- Deleted AN_Stahl from the Dvh file
- Deleted CI_36, CI_37 from Mmp file (since bogus as per email from Serdar)
### names_clonal_isolates_*.csv
- Changed Dvh/Mmp WT in the "evolution line" column to "AN".
### Issues
- The mutation file (\*allsamples_attributes\*) has lines for AN_\*1 and AN_\*2. It's unclear if these related to the clonal isolates with AN.
## Updates 03282018
- created XXX_mutations_allsmaples_03282018.xlsx, where XXX is either Dvh or Mmp
  - Removed CI_00_S7
  - Added the column name "other_predictor" for unnamed column T
  - Mmp
    - deleted CI_36, CI_37
  - Dvh
    - Deleted AN_00_S7_Stahl

## Updates 09042018
Contains refinements and sequencing of additional generations.

- The original data (in save) has samples with 118.early and
  118.new. The xlsx file only uses 118.early.
- FREQ1, FREQ2 are present in the raw data but deleted since they are not part of the data model
- Deleted references to other_predictor
- Changed DVU\_ to DVU
- Changed names_clonal_isolates_\*. "UA3" -> "UE3"
