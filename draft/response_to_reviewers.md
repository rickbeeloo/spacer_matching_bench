# Response to Reviewers

## General Response

We thank both reviewers for their thorough and constructive feedback on our manuscript "Tool Choice drastically Impacts CRISPR Spacer-Protospacer Detection." The comments significantly improve the quality and impact of our work, and have addressed all concerns raised in the revised manuscript. We note due to scope and format constraints not all of the suggestions are described in the main-text, and are instead described in this letter and added as an appendix to the supplementary material. 

Briefly, in addition to the reviewer suggestions, we modified the analyses:  

1. **New tools and updated versions:** Adding to the benchmark new tools (sassy, x-mapper and indelfree.sh) published after the original manuscript submission. We also updated the already included tools to their latest versions, or changed certain parameters based on reviewer suggestions. Notably, sassy was instrumental in testing the effects of using edit distance (allowing indels) versus hamming distance (substitutions only) up to 5 mismatches/edits.

2. **Enhanced synthetic dataset:** As suggested by reviewers, we tested additional variants of the synthetic dataset generation to include more realistic parameters, including contig length and GC distributions matched to IMG/VR4 contig set and Iphop spacer set. We include an additional jupyter notebook in the supplementary material that demonstrates the sequence characteristics (base frequencies, length, types of distributions, k-mer patterns, entrophy and complexity) between the real and synthetic spacer datasets. We also add a new "semi-synthetic" simulation option, with which real spacers are inserted into synthetic contigs, to test established a proxy for estimating non-specific matches in a more realistic setting (note, the inverse option - synthetic spacers into real contigs is also supported).

3. **Hamming vs Edit Distance Analysis:** We conducted comprehensive analysis comparing hamming distance (substitutions only) versus edit distance (allowing indels). For the edit distance measurements, we used sassy - the only tool supporting arbitrary edit distances with perfect recall, and bowtie1 (for hamming distance mesurment of very large sequence set) as previously demonstrated to have perfect recall. Using both the synthetic and semi-synthetic datasets, we demonstrate with these tools that:
   - Edit distance >3 leads to >10% false positives in synthetic data with approximately matched sequence characteristics (and similar number of expected spacers/protospacers).
   - Hamming distance >3 leads to >1% false positives in the synthetic data (compared to the total planned proto-spacer occurrences).
   - Hamming distance ≤3 provides optimal balance of recall and precision, with <1% false positives in the synthetic data (also compared to total planned proto-spacer occurrences), and we demonstrate that with the semi-synthetic dataset (measuring the occurnce of the full real spacer sequence set (~3.7M) in simulated contigs), with hamming <=3, 
   - Edit distance incurs massive computational costs (e.g., 4.7TB output for 5% subsample of the full IMG/VR contig set (~400k sequences) with ≤5 edits, ~1M CPU seconds).

4. **Revised false positive definition:** We corrected and clarified the terminology around distance metrics and redefined the false-positive and true-positive sets. The "spurious alignments" in the synthetic set (matches occurring in regions where we did not plan insertions) are now considered false-positives. We use their relative abundance as a proxy for expected false positive rate given a distance threshold (hamming or edit), similar in concept to an e-value: 'given this sequence length and target size, what is the expected number of matches with ≤ k mismatches by chance'.

5. **Biological justification for hamming distance:** We added substantial discussion of the biological basis for preferring hamming distance:
   - Literature review showing indels are ~4x rarer than substitutions in bacteria (the host side)
   - Experimental phage-host systems show most escape mutations are single mismatches
   - Phage genomes are coding-dense; single indels causing frameshifts are likely to be lethal.
   - While some indel escape mutations have been reported, they are not well quantified.
   - Large genomic rearrangements (another escape mechanism) would not be detectable with small edit distances anyway
   - With sufficient sequencing depth, sequencing-induced indels should be rare in assemblies

6. **Updated spacer dataset:** We replaced the CRISPR spacer set with a more curated, up-to-date set from iPHoP (June 2025 release), combining reference genomes and metagenomes.

7. **Stratified subsampling approach:** Partially due to computational constraints and to test tool behavior at different scales, we modified the "real" dataset analysis (IMG/VR4 contigs vs iPHoP spacers). Instead of using the full contig set, we employed a stratified subsampling strategy with fractions of varying sizes (0.001, 0.005, 0.01, 0.05, and 0.1 - albeit not all tools succesfully completed the analysis for all fraction in reasonable time) that maintains the distributions of contigs' taxonomic class labels. We demonstrate that within different subsampled fractions, the overall per-tool recall rates follow consistent trends, suggesting the subsample approach is adequate for informing tool choice.

**Manuscript restructuring:** Following reviewer suggestions, we restructured the manuscript into two main analytical sections: (1) "Comparing tool algorithms" - establishing that hamming distance is biologically appropriate and explaining why different algorithmic approaches report different numbers of matches, and (2) "Benchmarking hamming/heuristic-based tools" - comparing tools using hamming distance to evaluate speed and completeness. This clarifies that tools differ not because they are inherently "bad" but because their algorithmic assumptions/design goals may not align with certain biological questions.

## Response to Reviewer 1

### Major Concerns

**1. Enhanced Synthetic Dataset Analysis and Realistic Scenarios**

"the study felt like a missed opportunity to study how realistic scenarios and biased datasets affect the detection of spacer-protospacer pairs... the statistical properties of the synthetic dataset are very simplistic and far from representing real databases."

**Response:** We enhanced our synthetic dataset generation to address these concerns:

- **Heterogeneous contig lengths:** We now generate contigs with realistic length distributions based on IMG/VR4 data rather than uniform distributions. Specifically, we add the option to simulate the length in uniform or normal distribution, and in the simulated sets used in the new version, we roughly matched these distributions to the real (IMG/VR4) contig length range.
- **Compositional bias:** Similarly, we added options for the controlling the base composition, and matched the GC% to that of the iPHOP spacer data (~49% for spacers, ~46% for contigs). We ran additional comparisons to show that other sequence characteristics (k-mer distributions, repeat abundance etc) are also similar between the real and the new synthetic datasets, and included a notebook (real_spsacer_inspection.ipynb) in the supplementary material to demonstrate this.

**2. False Positive Analysis and Precision Metrics**
    "...I missed a more honest treatment of false positives. If you make a synthetic dataset, select pre-planned insertion points and then find out that the dataset, which is completely random (no sequence conservation, shared ancestry, or HGT), contains other identical or nearly identical sequences, then those can only be considered false positives. The authors even call these "spurious aligments", which is quite revealing of what they actually are (or what they are not: meaningful or informative). The fact that the number of those spurious alignments skyrockets with the mismatch threshold is also quite revealing: the greater the number of accepted mismatches, the greater the risk of getting "protospacers" that are not meaningful. The fraction of such hits in the synthetic dataset is quite remarkable (>10%) if one allows up to 3 mismatches. The fraction in a more realistic (synthetic) dataset with compositional bias would be even higher (because biases reduce the total entropy). From a biological point of view, such meaningless hits are false positives and should be treated as such. Thus, I was puzzled that the authors not only did not leverage the synthetic dataset to assess the specificity of different tools against false positives (which perhaps could make blastn-short perform better than bowtie1) but also counted them as true positives. I agree that doing that with the real dataset is complicated, because deciding whether a hit is meaningful or spurious would require further analyses of phylogeny and genomic context. But not doing that with the synthetic dataset seems harder to justify. In the end, finding spacer-protospacer pairs with mismatches comes at the cost of getting more false positives. Therefore, depending on the application, better performance finding pairs with 2-3 mismatches may not be that advantageous."

**Response:** We completely agree with the reviewer and have revised our analysis accordingly. The synthetic dataset now properly treats spurious alignments as false positives:

- **Redefined positive sets:** Only the planned occurrences of proto-spacers in the simulated contigs are now considered true positives (ground truth). The previously named "spurious alignments" (matches in regions where we did not plan insertions) are now classified as false-positives - we note that these are not "errors" of the tools (i.e. where a tool incorrectly reports a match at a given thershold, but the alignemnt of said match is above the thershold) but rather these are valid matches (match region aligns within requested thershold) that arise from chance due to sequence composition and size of the target database.

- **False positive rate as a heuristic:** We use the relative abundance of false positives as a proxy for expected FP rate given a distance threshold (hamming or edit), conceptually similar to an e-value: 'given this sequence length and target size, what is the expected number of matches with ≤ k mismatches by chance'.

- **Key findings on false positive rates:**
  - Hamming distance ≤3: <1% false positives (acceptable)
  - Hamming distance >3: >1% false positives
  - Edit distance >3 (allowing indels): >10% false positives
  - The reviewer's intuition was correct: the FP rate "skyrockets" with higher thresholds

- **Hamming vs Edit Distance comparison:** Using sassy (the only tool supporting arbitrary edit distances with perfect recall), we comprehensively tested both hamming and edit distances up to 5. This analysis definitively shows that allowing indels dramatically increases false positive rates while providing minimal additional biological information (see below).

- **Enhanced synthetic dataset realism:** We improved the synthetic dataset to better match real data:
  - GC content matched to iPHoP spacers (~49%)
  - Contig GC% and length distributions matched to IMG/VR4, 
  - Verified k-mer distributions and other sequence characteristics similar to real data
  - Results show the reviewer's concern about compositional bias is valid and important

**Results:** This analysis confirms the choice of hamming distance ≤3 as optimal, balancing recall against acceptable false positive rates. While bowtie1 maintains high recall, the precision analysis shows different tools have different trade-offs suitable for different applications.

**3. Real Dataset Limitations and Positive Set Definition**

The reviewer raises an important point: "The positive set for the real dataset is just the union of the protospacers identified by any of the tools. Thus, protospacers that are missed by all tools are completely ignored when estimating the recall."

**Response:** We have addressed this limitation in several ways:

- **Limitation discussion:** We explicitly discuss this limitation and its potential impact on recall estimates
- **Complementary analysis:** We perform additional analysis using the synthetic dataset where ground truth is known
- **Conservative estimates:** We provide conservative recall estimates and discuss the potential for missed true positives
- **Tool-specific analysis:** We analyze each tool's unique contributions to the positive set

**4. Tool Design Features and Performance Insights**

The reviewer asks for "insight that connected the heuristics and design features of different methods with their performance" and suggests "adding more information to Table 1."

**Response:** We have significantly enhanced our analysis and restructured the presentation:

- **Enhanced Table 1:** We now include tool categorization (exhaustive vs heuristic search, distance metric type: hamming vs edit/affine, original intended purpose, computational requirements)

- **Algorithm comparison framework:** Following reviewer suggestions, we restructured the analysis into two main sections:
  1. **"Comparing tool algorithms"** - We establish biological expectations (hamming distance is biologically appropriate), explain why different algorithmic approaches report different numbers of matches, and demonstrate that affine/edit-based algorithms will always report more matches than hamming-based ones - not because they are "bad" but because they solve a different computational problem, and partly due to the fact all sequence pairs that are within hamming distance ≤ k are also within edit distance ≤ k, but not vice versa.
  
  2. **"Benchmarking hamming-based tools"** - After establishing that hamming distance is appropriate for this biological question, we compare tools that use hamming distance on the full dataset to evaluate speed and completeness.

- **Heuristic analysis:** We explain how tool-specific heuristics affect performance:
  - Some tools penalize seeds from high-frequency k-mers, affecting high-abundance target detection
  - Different tools handle multiple matches differently (report all, report best, report up to threshold)
  - Index-based vs on-the-fly approaches have different computational trade-offs

- **Design feature connections:** We explicitly connect tool design choices to observed performance patterns, showing that tools differ not because they're "bad" but because their algorithms are not aligned with certain biological questions

- **Future method insights:** We provide recommendations for improving spacer-protospacer detection tools and note that understanding algorithmic assumptions is critical for appropriate tool selection

### Minor Points

**5. Methodological Clarifications**

We have addressed all methodological concerns:

- **Distance metric terminology and biological justification:** We clarify that our parasail-based alignment method differs from true Levenshtein distance and provide comprehensive justification for choosing hamming distance over edit distance:

  - **Computational demonstration:** We provide a notebook demonstrating that allowing indels as part of the distance metric increases the number of potential alignments drastically. Any string derivable from the original by ≤ k substitutions can also be derived by ≤ k edit operations, plus many additional strings containing indels.
  
  - **Empirical false positive analysis:** Using sassy (the only tool supporting arbitrary edit distances with perfect recall), we demonstrate that edit distance >3 leads to >10% false positives in synthetic data, while hamming distance ≤3 maintains <1% false positive rate.
  
  - **Biological mutation rates:** Literature evidence shows that indels are approximately 4× rarer than substitutions in bacteria. While we could not find quantitative measurements specifically in phages, this trend is expected to be similar or more pronounced.
  
  - **Experimental escape mutation evidence:** Most experimental phage-host studies report that escape mutations are predominantly single substitutions, particularly in the PAM-proximal "seed" region. While some indels have been reported, we could not find quantitative comparisons across mutation types.
  
  - **Phage genome structure:** Phage genomes are coding-dense with most of the genome covered by coding sequences. A single indel in a coding sequence causes a frameshift mutation, which is often lethal. In contrast, a substitution may only affect a single amino acid residue.
  
  - **Frame-preserving indels:** We did observe rare cases in the real data where indels of 3/6/9 bp (multiples of 3, preserving reading frame) could explain alignments. These cases, although extremely rare, may be particularly strong indicators of phage-host interaction due to the low probability of frame-preserving indels arising by chance.
  
  - **Sequencing vs biological indels:** We distinguish biological indels from sequencing artifacts. With sufficient sequencing depth, sequencing-induced indels should be rare in assembled contigs, as most reads would represent the correct sequence. This is particularly true for Illumina-based NGS datasets where consensus sequences (assembled contigs) are used. However, for other sequencing technologies (e.g., Oxford Nanopore) or low-depth datasets where raw reads are used, some tolerance for indels may be necessary. We note that while long-read technologies are tempting for CRISPR arrays (which are repetitive and prone to misassembly in short reads), the alignments cannot be better than the underlying sequencing accuracy.
  
  - **Assumption bias in literature:** We note that most existing studies focus on substitutions, which may create an "assumption bias" where indels are underexplored. Future experimental work could systematically compare escape mutation types (substitutions vs indels vs large genomic rearrangements).
  
  - **Large-scale rearrangements:** We acknowledge that large genomic rearrangements (e.g., when a spacer targets the boundary between two genes and gene order changes) exist but won't be captured by edit distance searches and are not well quantified in literature.

- **Synthetic dataset description:** We move detailed description to the Methods section and clarify the uniform random nature
- **Mismatch threshold definition:** We explicitly define mismatch thresholds throughout the manuscript
- **Levenshtein distance limitations:** We correct the discussion about indel handling
- **Minimap2 performance:** The reviewer correctly identifies that minimap2's near zero recall is due to parameter choice. The default chaining score threshold (-m 40) prevents matching for spacers shorter than 40nt. We tested minimap2 with lowered -m values (specifically `...-k 8 -w 6 -P -m 10`) or and indeed obtained more matches. However, minimap2 still underperformed considerably compared to most other tools. We added a sentence directly addressing this in the manuscript, noting that while minimap2 is apperntly not well-suited for CRISPR spacer-protospacer matching, it is still an excellent tool for its intended purpose, with a proven track record and a popular choice for read mapping (hence, its inclusion in the benchmark).
- **High target multiplicity and database redundancy:** The reviewer raises an important point about whether high-multiplicity spacers (>1000 targets) are informative, noting that "...spacers with >1000 targets are probably not very informative (unless the MGE database is extremely redundant)". We clarify that both the contig and spacer datasets used (in the original analysis and in this revision) are de-duplicated (non-redundant) at the sequence level (we explicitly only include unique sequences, removing duplicates or reverse-complement sequences). We now epxlictly mention this as part of our recommendation for dataset pre-filteration. We note that the high multiplicity reflects genuinely shared short sequences appearing across different contigs, not sequence database redundancy from identical full-length contigs. This phenomenon is biologically meaningful and has been documented in the literature. For example, "frozen prophages" - prophages appearing in certain bacterial lineages that are targeted by large CRISPR arrays from the same bacteria - create legitimate high-multiplicity scenarios, as these can align well with related prophage. We discuss this in the manuscript (referencing Stern et al. and others), noting that such cases represent real biological patterns rather than database artifacts. While we acknowledge that interpretation of very high-multiplicity spacers requires careful consideration (as they may represent conserved genomic regions rather than specific phage-host interactions), and that these are likely not the common case (i.e. we do not expect most spacers to occur at such multiplicity) their detection is still computationally relevant for benchmarking tool performance, particularly given that tools show dramatically different behaviors wrt to the spacer multiplicity.

**6. Figure and Presentation Issues**

- **Figure consistency:** Standardized colors and symbols across all figures.
- **Figure 2 caption:** Fixed the circular definition of "detection fraction"
- **Figure 3 readability:** Implemented log scaling and abbreviated notation
- **Section headers:** Reorganized Methods section structure
- **Figure references:** Corrected all figure references and descriptions

## Response to Reviewer 2

### Major Points

**1. Title and Computational Focus**

The reviewer suggests the title "should reflect computational nature of study and purposes of detection" and use "a more appropriate word than drastic."

**Response:** We have considered revising the title to better reflect the computational benchmarking nature, and suggest the new title "Computational Tool Choice Impacts CRISPR Spacer-Protospacer Detection" to use more appropriate language than "drastically."

**2. Enhanced BLASTn Analysis and Parameter Testing**

The reviewer emphasizes that "blastn results are most relevant to many current tools and analyses" and requests "more analysis of blastn's performance" ... "I note that some spacer matching tools e.g. CRISPRTarget use unusual blastn parameters for matching, albeit for a different purpose and on a smaller scale (https://crispr.otago.ac.nz/CRISPRTarget/crispr_analysis.html) some other sets of  parameters should also be tested e.g . wordsize if not as small as possible (7) or -evalue." 

**Response:** We clarify that the parameters used for BLASTn in our analysis were specifically chosen to enable the highest possible recall for short spacers - specifically `blastn-short` (`-task blastn-short`, sets word-size at 7, compared to 11 with default blastn, and 28 in megablastn) and a very high number of max target sequences (`-max_target_seqs 1000000`, so up to 1M reported alignemnts per query, compared the default value of 500). We add a paragraph in the methods section clarifying this choice compared to the default parameters and parameters used by certain tools and previous studies (e.g. the original IMG/VR4 spacer-protospacer mapping used default blastn with `-max_target_seqs 1000`). Regarding "CRISPRTarget" - associcated website (and link provided) are not accessible (ERR_CONNECTION_TIMED_OUT). The research paper discribing CRISPRTarget (Biswas et. al. 2013) mentions that the website allows modifying (selecting custom values for word size, evalue and choice of target database) for the BLASTn step, and using the defaults blastn-short values for word size (no mention of max-target-seqs), with a permissive e-value (to account for variations in database size) and modified gap penalties: "The initial CRISPRTarget defaults are the same except that a gap is penalized more highly (-10), the mismatch penalty is -1 and the E filter is 1. In addition, there is also no filter or masking for low complexity". These adjustments should favour gapless alignments compared to blastn-short defaults (gap open -5, gap extend -2). While we do not disagree with these parameter choices, we note that alignments passing these should still be reported by the default blastn-short parameters used in our analysis - the main expected change would be that with such parameters, gapped alignments would not need to be filtered out post-hoc (as done in our code, and is discussed in the choice of use of hamming distance of 3). We added a paragraph addressing this in the main text. Note - we do not see mention of using a modified "-max_target_seqs" in CRISPRTarget, and thus must assume the default (500) is used.


**3. DUST Filtering / Low Complexity Analysis**

The reviewer suggests testing "the effect of DUST" and asks "are the sequences with high abundance being missed of low complexity."

**Response:** We have implemented a comprehensive pre-filtering step (real_spacer_inspection.ipynb), where we used several metrics to pre-emptively measure and discard (prior to any mapping/alignment/benchmarking) low-complexity, repetitive, or homo-polymer spacers. We added a paragraph regarding this to the methods section, and we now explicitly note this as a crucial recommendation for future analyses. We further note that while some tools employ internal complexity filtering (e.g., DUST in blastn, tantan in mmseqs2), not all tools do so consistently. Therefore, applying such pre-filtering uniformly across all tools ensures comparable results and may affect downstream analysis.

**4. Tool Recommendations and Guidelines**

The reviewer notes that "the authors state they provide general guidelines in the abstract but these are less clear in the discussion."

**Response:** We have clarified and expanded our tool recommendations with specific use cases (when and why to use certain tools):

- **Primary recommendation - Bowtie1 (hamming distance ≤3):** For most large-scale metagenomic analyses, we recommend bowtie1 with hamming distance ≤3. This provides high recall (>99% for perfect to 3 substitutions spacers, a distance that we find biologically appropriate and having an acceptable false positive rates (<1%)). We note that bowtie1 offers good computational efficiency and scalability even for databases with millions of spacers and hundreds of thousands of contigs.

- **Extended hamming distance - Indelfree.sh indexed mode (hamming distance ≤5):** For scenarios requiring detection of spacers with up to 5 substitutions, we recommend indelfree.sh in indexed mode. While more computationally intensive than bowtie1, it maintains retains the near perfect recall and avoids the dramatic false positive increase associated with using edit distance. We note that it could be optimized further if the minimum spacer length is known in advance to be >=21 (which means a k-mer size of 8 can be used instead of 7). We still note that the false positive rate increases beyond 3 substitutions, so this tool is best suited for applications where extended mismatch tolerance is necessary and some increase in false positives is acceptable.

- **Sassy (edit distance) - Specific use cases only:** Sassy is the only tool evaluated with perfect recall that supports arbitrary edit distances. However, it incurs massive computational costs (e.g., ~1M CPU seconds, 4.7TB output for 5% subsample with ≤5 edits). We recommend sassy only for:
  - **Small datasets** where computational cost is acceptable
  - **Comparative studies** of escape mutation's types in controlled (exprimental) phage-host systems - we note that such studies are currently lacking as literature seems focused on measuiring substitution number and location, rather than on mutation type (substitution vs indel vs rearrangement) and frequncy. This positions sassy as the only current reccomendation for such studies, which alings with such experimental setups (at least historically) relatively small-scale datasets (e.g. several phages and several host strains). 
  - **Testing with ≥3 substitutions and indels** it might be possible that certain mutations accumulate both indels and substitions, as the actual host of a divergent phage may be different than the closest available host from which spacers were extracted.  
  - **Low-accuracy raw long reads** (e.g., Oxford Nanopore R9 chemistry) where indels are a larger concern and some minimal edit distance tolerance might be considered
  - **Methodological validation** to establish baseline performance when developing new tools, or deciding on benchmarking guidelines (such as in this work)

- **Sassy limitations:** We explicitly note that sassy's computational requirements (1M CPU seconds, 4.7TB outputs for just 5% of data) make it prohibitive for large-scale metagenomic analyses. The massive output size also presents storage and downstream processing challenges. We note this is not a limitation of sassy as a tool, rather a limitation of using edit distance and exhaustive (non-heuristic) search for this biological question at large scale.

- **Hamming vs Edit Distance trade-offs:** We demonstrate empirically that:
  - Hamming distance ≤3: <1% false positives
  - Hamming distance >3: >1% false positives
  - Edit distance >3: >10% false positives
  
  This dramatic increase in false positives with edit distance is not justified by the minimal additional biological information gained, given that indels are ~4× rarer than substitutions and most escape mutations are single substitutions.

- **Use case scenarios:** We provide specific recommendations:
  - **High-throughput host-virus prediction:** bowtie1 with hamming ≤3
  - **Extended mismatch tolerance:** indelfree.sh indexed with hamming ≤5
  - **Small experimental datasets:** sassy with appropriate distance threshold
  - **Precision-focused applications:** Consider DUST filtering and conservative thresholds

- **Resource trade-offs:** We discuss computational requirements comprehensively, showing that bowtie1 and indelfree.sh indexed mode scale well to large databases, while sassy is limited to small-scale applications.

- **Parameter guidance:** We provide specific parameter recommendations for different tools and discuss how parameter choices affect performance, particularly for BLASTn where `-max_target_seqs` significantly impacts high-abundance spacer detection.

### Minor Points

**6. Figure Improvements**

We have addressed all figure issues:

- **Line distinction:** Improved line styles and colors for better distinction
- **Figure 3 readability:** Implemented log scaling and better color schemes
- **Supplementary figures:** Fixed typos and improved legends
- **Color consistency:**  Colors and symbols representing tools should now be consistent across figures.


**7. Historical Comparison**
> "which would "allow for a direct comparison to historical results". Was that comparison ever made? I could not find any further mention to that in the manuscript."   

**Response:** Our choice to benchmark tools using the IMG/VR4 contig dataset inherently serves as a comparison with the original IMG/VR4 publication, which used spacer-protospacer matching for host assignment (using default BLASTn with `-max_target_seqs 1000`). By documenting the differences in methodology and parameters (particularly of `blastn-short` with `-max_target_seqs 1000000`), we provide context for understanding how tool and parameter choices affect results. However, we note that comparing the specific nature of host assignments falls outside the scope of this sequence alignment/search tool benchmark. We recommend that historical host assignments be reconsidered in light of our findings, but detailed validation of those assignments is beyond the scope of this manuscript.


## Supplementary Material
All new notebooks, datasets, and analysis scripts are provided in the public github repository. The final mapping results (i.e. files too large for github) are uploaded as a new version to the same zenodo deposit.