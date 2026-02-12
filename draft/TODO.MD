# TODO: Addressing Reviewer Comments

### 1. Enhanced Synthetic Dataset Analysis
**Reviewers:** R1 (points 10, 14, 15, 16), R2 (point 9)
**Tasks:**
- [x] **Create more realistic synthetic datasets** with:
  - [x] Heterogeneous contig lengths (not uniform)
  - [x] Compositional bias (non-uniform GC content)
- [ ] **Implement false positive analysis** in synthetic data:
  - [x] Treat spurious alignments as false positives (not true positives)
  - [x] use spurious alignments on simulated data as a proxy for an FP rate of equivalent real data
  - [x] Calculate precision metrics alongside recall
  - [ ] Analyze sensitivity vs specificity trade-offs (not sure this is really needed)

### 2. Real Dataset Limitations and Analysis
**Reviewers:** R1 (points 12, 36, 37), R2 (point 7)
**Tasks:**
- [ ] **Address positive set definition limitations:**
  - [ ] Discuss that union of all tools may miss true positives
  - [ ] Estimate frequency of missed cases
  - [ ] Add discussion of this limitation
- [ ] **Clarify false negative identification** in real dataset
- [ ] **Test DUST filtering effects:**
  - [ ] Apply DUST masking to sequences
  - [ ] Compare performance with/without masking
  - [ ] Analyze if high-abundance sequences are low-complexity

### 3. Tool Performance Analysis and Parameter Testing
**Reviewers:** R1 (points 18, 32), R2 (points 5, 6, 7)
**Tasks:**
- [x] **Add new tools:**
  - [x] sassy (perfect recall, supports arbitrary edit distances)
  - [x] x-mapper
  - [x] indelfree.sh (both bruteforce and indexed versions)
- [x] **Sassy-specific analysis:**
  - [x] Test edit distances up to 5
  - [x] Document computational requirements (1M CPU seconds, 4.7TB output)
  - [x] Show perfect recall but prohibitive for large datasets
  - [x] Recommend for small datasets and certain experimental setups (in discussion)
- [ ] **Fix minimap2 performance:**
  - [x] Test with lower chaining score threshold (-m parameter)
  - [x] Document other parameter changes 
- [ ] **explain BLASTn-short"est" analysis:**
  - [ ] Document default parameters used (and compare with CRISPRTarget parameters).
- [ ] **Better Tool categorization and heuristics analysis:**
  - [ ] Add to introduction brief overview of algorithms (hamming vs edit-based)
  - [ ] Add tool design features to Table 1 (exhaustive vs heuristic)
  - [ ] Explain why tools perform differently (Discussion) (Connect heuristics / algo to performance patterns )
  - [ ] Restructure as Rick suggested: "Comparing tool algorithms" then "Benchmarking hamming/heuristic based tools"

### 4. Precision vs Recall Analysis
**Reviewers:** R1 (points 14, 15, 16), R2 (point 8)
**Tasks:**
- [x] **Implement precision calculations:**
  - [x] Calculate precision for each tool using simulated data
  - [x] Analyze trade-offs for different use cases (hamming vs edit distance)
  - [x] Demonstrate that edit distance >3 leads to >10% false positives
  - [x] Show that hamming distance >3 leads to >1% false positives
- [x] **Hamming vs Edit Distance Analysis (using sassy):**
  - [x] Test edit distances up to 5
  - [x] Compare FP rates between hamming and edit distance
  - [x] Demonstrate computational costs (sassy: 4.7TB output for 5% subsample)
- [x] **Tool recommendations for different scenarios:**
  - [x] High precision scenarios: hamming distance ≤3
  - [x] High recall scenarios: sassy for small datasets with edit distance
  - [x] Large-scale scenarios: bowtie1/indelfree_indexed with hamming ≤3
  - [x] Document sassy's perfect recall but computational limitations

### 5. Methodological Clarifications
**Reviewers:** R1 (points 24, 26, 28, 38, 39), R2 (point 2)
**Tasks:**
- [x] **Define mismatch threshold explicitly**
- [x] **Fix distance metric terminology:**
  - [x] Explain choice of alignment method (hamming vs edit)
  - [x] Clarify that our parasail re-alignments as we do them are explicitly set to calculate hamming distance, penalizing gaps.
  - [x] Document that edit distance allows indels, hamming only substitutions
  - [x] Show empirically that edit distance increases FP rate dramatically
- [x] **Biological justification for hamming distance:**
  - [x] Literature review: indels rarer than substitutions (~4x in bacteria, add citations)
  - [x] Experimental evidence: most escape mutations are 1 mismatch
  - [x] Phage biology: coding-dense genomes make indels likely lethal (frameshift)
  - [x] Note: some indels reported but not quantified in literature 
  - [x] Large genomic rearrangements exist but won't be caught by edit distance
  - [ ] Add citations for escape mutation studies
- [x] **Improve synthetic dataset description:**
  - [x] Move detailed description to Methods section
  - [x] Document GC% matched to real data (~49% for spacers, ~46% for contigs)
  - [x] Document length and GC distributions matched to IMG/VR4
- [ ] **Sequencing errors discussion:**
  - [ ] Argue that with sufficient depth, sequencing-induced indels should be rare
  - [ ] Distinguish biological mutations from sequencing artifacts

### 6. Historical Comparison
**Reviewers:** R1 (point 22), R2 (point 4)
**Tasks:**
- [ ] **Perform IMG/VR4 historical comparison:**
  - [ ] Compare with original IMG/VR4 spacer-protospacer results
  - [ ] Document differences in methodology
  - [ ] Add results to manuscript
- [ ] **Update BLAST discussion:**
  - [ ] Remove outdated 2018 issues
  - [ ] Focus on current parameter considerations

### 7. Figure and Presentation Improvements
**Reviewers:** R1 (points 40, 42, 43, 44, 48, 50, 52, 54), R2 (point 6)
**Tasks:**
- [ ] **Fix figure inconsistencies:**
  - [x] Standardize colors/symbols across figures
  - [ ] Fix Figure 2 caption (detection fraction definition)
  - [ ] Improve Figure 3 readability (log scaling, abbreviations)
- [ ] **Fix section headers:**
  - [ ] Reorganize Methods section structure
  - [ ] Fix section 4.2 header terminology
- [ ] **Clarify figure references:**
  - [ ] Fix Figure 2 reference for perfect matches
  - [ ] Clarify contig dependent vs independent recall

### 8. Paper structure Improvements
**Reviewers:** R2 (points 1, 8)
**Tasks:**
- [x] **Revise title:**
- [ ] **Clarify guidelines in discussion:**
  - [x] **Hamming vs Edit Distance:** Demonstrate empirically that edit distance >3 leads to >10% FPs
  - [x] **Biological justification:** Indels rarer than substitutions, most escape mutations are 1 mismatch, phage genomes coding-dense (indels likely lethal)
  - [x] **Sassy recommendation:** Only tool supporting arbitrary edit distances; perfect recall; recommend for:
    - Small datasets where computational cost is acceptable
    - Experimental setups where mutation type is of interest
    - Comparative studies of escape mutation types (substitution vs indel vs rearrangement)
    - Testing with ≥3 substitutions (limitation of some other tools)
    - if using low-accuracy raw long reads such as Oxford Nanopore, where indels are a larger concern, some minimal edit distance might be considered.
  - [x] **Sassy limitations:** Note computational costs (1M CPU seconds, 4.7TB outputs) make it prohibitive for large-scale analysis
  - [x] **Large-scale recommendation:** indelfree_indexed and bowtie1 for high speed and high recall with hamming ≤3
  - [x] **False positive rates:** >4 substitutions leads to >1% FPs; >3 edits leads to much higher FPs
  - [ ] **Literature support:** Add citations for:
    - Escape mutation studies showing predominantly 1 mismatch
    - Indel rarity in bacteria (4x less common than substitutions)
    - Experimental phage-host systems
  - [ ] **Sequencing vs biology:** Distinguish biological indels from sequencing errors; argue sequencing-induced indels are rare with sufficient depth for most illumina based NGS datasets, especially when consensus (assembled contigs) sequences are uses. Note that for other sequencing technologies (e.g., nanopore) or low depth datasets, or when raw reads are used, some indels may need to be accounted for. Note that it is tempting to use long-read technologies as CRISPR arrays are repetitive and thus prone to misassembly in short-read data - but that eveutally the alignemnts can not be better than the underlying sequencing accuracy, and than is another subject altogether.  
  - [ ] **PAM-proximal mutations:** Note that most studies focus on substitutions, and that these seem to identify that substitions closer to PAM ("seed" region) are most observed in escape mutations. Note that these historical studies might suffer from "assumption bias" as they revolve around substitutions - suggesting that future experimental work could explore whether indels are at all more common, compared to substitutions and large genomic rearrangements.

### 9. Supplementary Data Organization
**Reviewers:** R1 (points 40, 41), R2 (point 6)
**Tasks:**
- [ ] **Present synthetic results in main text:**
  - [ ] Move key synthetic results to section 4.2
  - [ ] Explain multiplicity differences between datasets
- [ ] **Fix supplementary figures:**
  - [x] Fix "percision" typo in S3
  - [ ] Improve S4 legend visibility
  - [ ] Better organize supplementary data

### 10. Additional Analysis
**Reviewers:** R1 (point 30), R2 (point 9)
**Tasks:**
- [ ] **Reference genome analysis:** - explain these spacers (now) are from iphop, combining reference genomes and metagenomes, and that the virus contigs are also from both.
  - [x] Test with spacers from reference genomes
  - [x] Compare with metagenomic results
- [ ] **Address high multiplicity:**
  - [ ] Discuss biological relevance of >1000 target spacers ("frozen" prophages)
  - [ ] Focus on 1-1000 range for most applications

### 11. Methods Restructuring (Rick's suggestion)
**Tasks:**
- [ ] **Section 1: "Comparing tool algorithms" (or better heading)**
  - [ ] Set biological expectation: hamming distance is expected
  - [ ] Explain why some approches will allways report more matches (affine/edit-based algorithms vs heuristic/hamming)
  - [ ] Brief description of what different algorithms do and implications
  - [ ] Use the subsamples of real data, the synthetic and semi-synthetic data to showcase the changes in False Positives when using different distance metrics and values.
  - [ ] show that hamming distance is more aligned with biological expectation at a substantial reduction in FPs compared to edit distance.
  - [ ] Run the tools that can finish in reasonable time on the complete real dataset.
  - [ ] Show tools differ not because they're "bad" but solve different problems (in discussion).
  - [ ] Frame as "algorithm not aligned with certain biological question" (in discussion).
- [ ] **Section 2: "Benchmarking heuristic/hamming/edit based tools"**
  - [ ] After subsample establishes we want hamming distance
  - [ ] Identify which tools (X, Y, Z) use hamming distance
  - [ ] Run on full dataset
  - [ ] Compare speed and results
 Then have an overview in the discussion about which tools to use for what purpose and in what scenarios.
