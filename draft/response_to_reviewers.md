# Response to Reviewers
We thank both reviewers for their thorough and constructive feedback on
our manuscript "Tool Choice drastically Impacts CRISPR
Spacer-Protospacer Detection." The comments have significantly improved the
quality and impact of our work. In this revised version, we present a large overhaul of the
project that addresses all concerns raised. We
note that due to scope and format constraints, not all suggestions are
described in the main text, and are instead described in this letter and
added as an appendix to the supplementary material.

Briefly, following the reviewer suggestions we modified the analyses:

1.  **New tools and updated versions:** We added to the benchmark new
    tools (sassy, x-mapper and indelfree.sh) published after the
    original manuscript submission (and the preprint version). We also
    updated the already included tools to their latest versions, or
    changed certain parameters based on reviewer suggestions. Notably,
    sassy was instrumental in testing the effects of using edit distance
    (allowing indels) versus hamming distance (substitutions only) up to
    larger mismatches/edits.
2.  **More detailed overview of the tested tools:** We now provide a brief
    overview in the introduction section of the general logic and types
    of sequence alignment/search tools. We also updated table 1 to
    include some of these features. Briefly, we explain the difference
    between exhaustive and heuristic tools, the different types of
    alignments (hamming, gap-affine and edit) they focus on, and the
    computational resources (namely memory and CPU time) they require.
3.  **Explicitly addressing resource usage:** Because we have added two
    exhaustive tools (sassy and indelfree\_bruteforce) to the analysis
    (when possible, see note below about subsampling), we now include a
    section with the runtime/memory usage (vs distance threshold and
    approximate search space size), and we explicitly discuss the
    tradeoffs of heuristic vs exhaustive tools.
4.  **Explicitly suggesting tool choice:** To clarify the new results and
    different scenarios where spacer-protospacer might be investigated,
    we have added a flowchart (diagram, figure nnn) that factors in the
    dataset size, the availability of validation / experiment design, and
    the expected performance.
5.  **Updated spacer dataset:** We replaced the "real-data" CRISPR
    spacer set with a more curated, up-to-date set from iPHoP (June 2025
    release), combining spacers from arrays detected in reference
    genomes and metagenomes. We now explicitly recommend that users and
    potential tool developers employ minimal low-complexity filtering.
    
    From a benchmark perspective, this helps isolate tool-specific
    complexity handling (not all tools allow manually disabling such
    settings). From a practical perspective, we identified that certain
    exhaustive tools (or tools with immutable cost matrices) create
    massive outputs when aligning sequences containing ambiguous "N" bases.
    
    In the current version of the "real-data" analysis, we employ very
    light complexity filtering: no ambiguous bases, no spacers where a
    single base frequency exceeds 90% of the sequence (e.g., homopolymers
    such as polyA/T/G/C), and no spacers with highly repetitive
    subsequences.
6.  **Stratified subsampling approach:** Partially due to computational
    constraints and to test tool behavior at different scales, we
    modified the "real" dataset analysis (IMG/VR4 contigs vs iPHoP
    spacers). Instead of using the full contig set, we employed a
    stratified subsampling strategy with fractions of varying sizes
    (0.001, 0.005, 0.01, 0.05, 0.1 and 1 - albeit not all tools
    successfully completed the analysis for all fractions in reasonable
    time) that maintains the distributions of contigs' taxonomic class
    labels. We demonstrate that within different subsampled fractions,
    the overall per-tool recall rates follow consistent trends for
    recall and resource usage. A notable exception were strobealign,
    mmseqs, and x-mapper. However, even in subsamples where these
    performed relatively well (compared to themselves), they are not
    currently compelling tools for this task. We therefore avoid
    investigating their inconsistencies as it is outside the scope of
    this project. The same consistency for the majority of the tools
    was observed in the different fully simulated datasets (see Supplementary
    Figures).
7.  **Enhanced synthetic dataset:** As suggested by reviewers, we
    modified the synthetic dataset generation to include more realistic
    parameters, including contig length and GC distributions matched to
    IMG/VR4 viral contig set and the Iphop spacer set. We include an
    additional jupyter notebook in the supplementary material that
    demonstrates the sequence characteristics (nucleobase frequencies,
    lengths, types of distributions, k-mer patterns, entropy and
    complexity) between the real and synthetic spacer datasets. Briefly,
    overall the real and simulated sequences appear to have similar
    feature distributions, although the real spacer set tends to include
    some low complexity sequences (albeit these are relatively rare, and
    we now attempt to minimize their effect on the analyses).
    
    We also added a new "semi-synthetic" simulation run, where the
    presence of real spacers in fully synthetic contigs is tested. We
    use this to estimate the rate of non-specific matches (previously
    termed "spurious" matches) in a more realistic setting.

8.  **Hamming vs Edit Distance Analysis:** We corrected and clarified
    the terminology around distance metrics, true-positive sets, and
    invalid sets. Specifically, the previously termed "spurious
    alignments" in the synthetic set (matches occurring in regions
    where the simulation did not explicitly insert a spacer) are now
    used as a proxy for estimating the rate of matches arising by chance,
    given a distance threshold (hamming or edit). This is conceptually
    similar to an e-value: 'given a query sequence and a target sequence
    set, what is the expected number of matches with ≤ k distance by chance?'
    
    We conducted a comprehensive analysis comparing hamming distance
    (substitutions only) versus edit distance (allowing indels), and
    reviewed existing literature on phage "escape" mutations (added to
    the discussion and introduction). For edit distance measurements, we
    used sassy—the only tool supporting arbitrary edit distances with
    perfect recall (exhaustive)—and bowtie1 for hamming distance
    measurement of very large sequence sets (previously demonstrated to
    have good recall).
    
    Using both synthetic and semi-synthetic datasets, we estimated the
    frequency of valid-but-unplanned matches (i.e., alignments verified
    to be within a given distance threshold, regardless of which tool
    predicted them, but not explicitly inserted in the simulated contigs).
    Based on the simulated runs, we observe that such unplanned matches
    are typically an order of magnitude more frequent when searching up
    to n edit distance compared to n hamming distance, and this trend
    increases considerably above distance 3 (see Supplementary Figure nnn).
    
    While this is expected (all hamming alignments are a subset of all
    edit-distance alignments), the exact rate and difference in the
    number of permitted alignments under these metrics is difficult to
    estimate. For practical purposes, a value is often approximated for
    the frequency of matches arising by chance (e.g., e-value). Our use
    of empirical data (the simulation runs) is partially justified by the
    new simulation's ability to mimic real biological sequences more
    closely (with regard to base distribution, lengths, etc.). For
    example, using the semi-synthetic dataset (measuring the occurrence
    of ~3.7M real spacers in 400k simulated contigs, matching the IMG/VR4
    set size), we identified specific numbers of real spacers forming
    alignments within hamming distance ≤3 [specific values to be added].

9.  **Biological justification for hamming distance:** We added a brief
    literature review and argument for preferring hamming distance in most
    scenarios:

    -   **Experimental escape mutation evidence:** Most experimental
        phage-host studies report that escape mutations are
        predominantly single substitutions, particularly in the
        PAM-proximal "seed" region. While some indels have been
        reported, we could not find quantitative comparisons across
        mutation types. We only identified one study mentioning indel
        escape mutations in protospacers, however that study does not
        perform a comparative quantification of them, and to the best of
        our understanding, the provided information in the paper [citation nnn]
        mentions only a single indel mutation.
    -   **Phage genome structure:** Phage genomes are coding-dense with
        most of the genome covered by coding sequences. A single indel
        in a coding sequence causes a frameshift mutation, which is
        often lethal. In contrast, a substitution may only affect a
        single amino acid residue.
    -   **Frame-preserving indels:** Literature suggests indels are ~4x
        rarer than substitutions in bacteria (the host), and we
        expect this to be similar in phages.
    -   **Sequencing vs biological indels:** We distinguish biological
        indels from sequencing artifacts. With sufficient sequencing
        depth, sequencing-induced indels should be rare in assembled
        contigs, as most reads would represent the correct sequence.
        This is particularly true for Illumina-based NGS datasets where
        consensus sequences (assembled contigs) are used. However, for
        other sequencing technologies (e.g., Oxford Nanopore) or
        low-depth datasets where raw reads are used, some tolerance for
        indels may be necessary. We note that while long-read
        technologies are tempting for CRISPR arrays (which are
        repetitive and prone to misassembly in short reads), the
        alignments cannot be better than the underlying sequencing
        accuracy.
    -   **Assumption bias in literature:** We note that most existing
        studies focus on substitutions, which may create an "assumption
        bias" where indels are underexplored. Future experimental work
        could systematically compare escape mutation types
        (substitutions vs indels vs large genomic rearrangements).
    -   **Large-scale rearrangements:** Another reported escape
        mechanism is large genomic rearrangements (e.g., when a spacer
        targets the boundary between two genes and gene order changes)
        [citation nnn]. While these are rare and theoretically identifiable
        by gap-affine methods, they would not be detectable with a small
        edit distance threshold (considering that these rearrangements
        involve substantial gene order changes).

We explicitly address cases where edit distance (or gap-aware tools)
should be considered, such as when using low-accuracy long read data
(e.g., Oxford Nanopore sequencing) or when investigating escape mutation
types under verifiable experimental settings.

## Response to Reviewer 1

### Major Concerns

**1. Enhanced Synthetic Dataset Analysis and Realistic Scenarios**

"the study felt like a missed opportunity to study how realistic
scenarios and biased datasets affect the detection of spacer-protospacer
pairs... the statistical properties of the synthetic dataset are very
simplistic and far from representing real databases."

**Response:** We enhanced our synthetic dataset generation to address
these concerns:

-   **Heterogeneous contig lengths:** We now generate contigs with
    realistic length distributions based on IMG/VR4 data rather than
    uniform distributions. Specifically, we added (to the CLI tool) the
    option to simulate the length in uniform or normal distribution, and
    in the simulated sets used in the new version, we roughly matched
    these distributions to the real (IMG/VR4) contig length range.
-   **Compositional bias:** Similarly, we added options for
    controlling the base composition, and matched the GC% to that of the
    iPHOP spacer data (~49% for spacers, ~46% for contigs). We ran
    additional comparisons to show that other sequence characteristics
    (k-mer distributions, repeat abundance, etc.) are also similar between
    the real and the new synthetic datasets, and included a notebook
    (spacer_inspection.ipynb) in the supplementary material to
    demonstrate this (from which Supplementary Figures nnn were generated).

**2. False Positive Analysis and Precision Metrics** "...I missed a more
honest treatment of false positives. If you make a synthetic dataset,
select pre-planned insertion points and then find out that the dataset,
which is completely random (no sequence conservation, shared ancestry,
or HGT), contains other identical or nearly identical sequences, then
those can only be considered false positives. The authors even call
these "spurious alignments", which is quite revealing of what they
actually are (or what they are not: meaningful or informative). The fact
that the number of those spurious alignments skyrockets with the
mismatch threshold is also quite revealing: the greater the number of
accepted mismatches, the greater the risk of getting "protospacers\" that
are not meaningful. The fraction of such hits in the synthetic dataset
is quite remarkable (\>10%) if one allows up to 3 mismatches. The
fraction in a more realistic (synthetic) dataset with compositional bias
would be even higher (because biases reduce the total entropy). From a
biological point of view, such meaningless hits are false positives and
should be treated as such. Thus, I was puzzled that the authors not only
did not leverage the synthetic dataset to assess the specificity of
different tools against false positives (which perhaps could make
blastn-short perform better than bowtie1) but also counted them as true
positives. I agree that doing that with the real dataset is complicated,
because deciding whether a hit is meaningful or spurious would require
further analyses of phylogeny and genomic context. But not doing that
with the synthetic dataset seems harder to justify. In the end, finding
spacer-protospacer pairs with mismatches comes at the cost of getting
more false positives. Therefore, depending on the application, better
performance finding pairs with 2-3 mismatches may not be that
advantageous.\"

**Response:** We completely agree with the reviewer and have revised our
analysis accordingly. As mentioned above, we have clarified the terminology and
metrics used. We now use these "spurious alignments" (after validating
that they are not erroneous tool reports) as a proxy for the rate of
non-informative matches with regard to biological origin (such as would
be expected for true host-phage interactions).

As the reviewer notes, the ">10%" rate is indeed quite remarkable. However,
this percentage estimates the rate of non-planned matches relative to
planned spacer occurrences. Since we control the number of planned
occurrences, this percentage is not the most useful metric. We now
instead report either the absolute number of such cases or their
frequency normalized by the target (contigs) and query (spacers) set sizes.

Specifically, using the semi-synthetic dataset (~3.7 million real spacers
as queries, searched against 400k simulated contigs matching IMG/VR4
sequence characteristics), we identified: 54,388 unique alignments at
hamming distance ≤3, 2,217 unique alignments at ≤2, 47 alignments at ≤1,
and 1 exact match (distance 0). We observe similar rates (same order of
magnitude, when adjusted for sequence length) in the fully synthetic runs.

**3. Real Dataset Limitations and Positive Set Definition**

The reviewer raises an important point: "The positive set for the real
dataset is just the union of the protospacers identified by any of the
tools. Thus, protospacers that are missed by all tools are completely
ignored when estimating the recall."

**Response:** We now note this limitation explicitly in the text and
explain that the majority of the tools tested are not exhaustive (i.e.,
they use heuristics, which means they are not guaranteed to report all matches).

We have specifically added subsample analyses small enough to include
exhaustive tools, and their results are part of the union of all results.
For larger subsets, these tools were unable to finish even after 1M
CPU seconds, so we do not include them in those plots. Overall, the
relative recall rates between tools appear consistent with and without
adding the union of the exhaustive tool results. That said, we now
explicitly mention that without the exhaustive tool results, we are
likely underestimating the actual number of valid alignments.

**4. Tool Design Features and Performance Insights**

The reviewer asks for "insight that connected the heuristics and design
features of different methods with their performance" and suggests
"adding more information to Table 1."

**Response:** We have significantly enhanced our analysis and
restructured the presentation:

-   **Enhanced Table 1:** We now include tool categorization (exhaustive
    vs heuristic search, distance metric type: hamming vs edit/affine,
    original intended purpose, computational requirements etc).
-   **Algorithm comparison framework:** Following the reviewer's suggestions,
    we now explicitly state when we are able to compare results to exhaustive
    tools. We explain how tool-specific heuristics can affect
    performance, noting some common examples (such as how penalizing
    seeds from high-frequency k-mers can reduce detection of high-abundance
    targets). We also discuss how some tools are designed to work with
    pre-generated indices, and how this affects runtime and memory
    requirements (e.g., whether the index can be fetched remotely or
    must be computed on the fly).
-   **Clarifying intended use and biological relevance:** We added a
    brief overview of how biological expectations relate to tool design
    and intended use. Specifically, we argue that hamming distance is
    appropriate in most cases, explain why different algorithmic
    approaches may limit the number of reported matches (as both a
    heuristic strategy and to match user expectations), and note that
    affine/edit-based algorithms will always report more matches than
    hamming-based ones. This addition clarifies that we do not
    consider any tool inherently "bad," but rather as tools designed to solve
    different computational problems that may not align with all analysis goals.

### Minor Points

**5. Methodological Clarifications**

We have addressed all methodological concerns:

-   **Distance metric terminology and biological justification:** We
    clarify that the parasail-based alignment method differed from true
    Levenshtein distance, and in the new version, we have compared and
    used hamming distance only when the underlying alignment did not
    include gaps. We added explicit options for the CLI tool for
    gap-affine (with cost matrix, gap open, and extension values),
    minimal edit distance (using edlib), and the simplified, correct
    hamming.
-   **Minimap2 performance:** The reviewer correctly notes that
    minimap2 displays near-zero recall. We adjusted the parameters as
    suggested (lowered -m values: specifically -k 8 -w 6 -P -m 10) and
    indeed obtained more matches. However, minimap2 still underperformed
    considerably, displaying near-zero recall compared to most other
    tools. We added a sentence directly addressing this in the manuscript,
    noting that while minimap2 is apparently not well-suited for CRISPR
    spacer-protospacer matching, it remains an excellent tool for its
    intended purposes and is very popular for read mapping and similar
    tasks (hence its inclusion in the benchmark). We note that we
    modified parameter choices for several tools, clarifying in the
    manuscript that we selected parameters to maximize recall based on
    our understanding of the tools and their guidelines (documentation,
    GitHub issue threads, etc.).
-   **High target multiplicity and database redundancy:** The reviewer
    raises an important point about whether high-multiplicity spacers
    (\>1000 targets) are informative, noting that "...spacers
    with \>1000 targets are probably not very informative (unless the
    MGE database is extremely redundant)". While these are rare (fig
    nnn), we clarify that both the contig and spacer datasets used (in
    the original analysis and in this revision) are de-duplicated
    (non-redundant) at the sequence level (we explicitly only include
    unique sequences, removing duplicates or reverse-complement
    sequences - we now explicitly mention this as part of our
    recommendation for dataset pre-filtering). We note that the high
    multiplicity reflects genuinely shared short sequences appearing
    across different contigs, not sequence database redundancy from
    identical full-length contigs. We argue this phenomenon is
    biologically meaningful and has been documented in the literature.
    For example, "frozen prophages" - prophages appearing in certain
    bacterial lineages that are targeted by large CRISPR arrays from the
    same bacteria - create legitimate high-multiplicity scenarios, as
    these can align well with related prophage. We discuss this in the
    manuscript (referencing Stern et al. and others), noting that such
    cases represent real biological patterns rather than database
    artifacts. While we acknowledge that interpretation of very
    high-multiplicity spacers requires careful consideration (as they
    may represent conserved genomic regions rather than specific
    phage-host interactions), and that these are not the common case, their
    detection is still computationally relevant for benchmarking tool
    performance, particularly given that tools show dramatically
    different behaviors with respect to the spacer multiplicity. In the revised
    manuscript, we still mention the multiplicity as a potential cause
    for tool differences, but we do not frame it as the sole or most
    important determinant of tool choice or effect on tool recall.

**6. Figure and Presentation Issues**

-   **Figure consistency:** Standardized colors and symbols across all
    figures.
-   **Figure 2 caption:** Fixed the circular definition of "detection
    fraction"
-   **Figure 3 readability:** Implemented log scaling and abbreviated
    notation
-   **Section headers:** Reorganized Methods section structure
-   **Figure references:** Corrected all figure references and
    descriptions

## Response to Reviewer 2

### Major Points

**1. Title and Computational Focus**

The reviewer suggests the title "should reflect computational nature of
study and purposes of detection" and use "a more appropriate word than
drastic."

**Response:** We have revised the title to better reflect the
computational benchmarking nature, and suggest the new title
"Computational Tool Choice Impacts CRISPR Spacer-Protospacer Detection"
to use more appropriate language than "drastically."

**2. Enhanced BLASTn Analysis and Parameter Testing**

The reviewer emphasizes that "blastn results are most relevant to many
current tools and analyses" and requests "more analysis of blastn's
performance" ... "I note that some spacer matching tools
e.g. CRISPRTarget use unusual blastn parameters for matching, albeit for
a different purpose and on a smaller scale
(https://crispr.otago.ac.nz/CRISPRTarget/crispr\_analysis.html) some
other sets of parameters should also be tested e.g . wordsize if not as
small as possible (7) or -evalue."

**Response:** We clarify that the parameters used for BLASTn in our
analysis were specifically chosen to enable the highest possible recall
for short spacers - specifically blastn-short (-task blastn-short, sets
word-size at 7, compared to 11 with default blastn, and 28 in
megablastn) and a very high number of max target sequences
(-max\_target\_seqs 1000000, so up to 1M reported alignments per query,
compared to the default value of 500). After some internal experimentation,
we have updated the parameters when calling blastn, adding the
\`\--ungapped\` flag, and we set perc\_identity=84 and
qcov\_hsp\_perc=80. The ungapped flag is used by CRISPRTarget, aligns
with the current project (see above for edit vs hamming etc) and
considerably reduced total cpu time. The perc\_identity=84 and
qcov\_hsp\_perc=80 arguments are mostly used to control the output size,
as we do not know of a way to restrict blastn to a specific "n"
distance. These values are mostly to control the output size, and should
allow for a range of distances to be reported - this range exceeds our
desired threshold in order to capture all valid matches (for example,
85% identity, over at least 80% of the spacer, will allow for 3 distance
for a spacer of length 20 (20 \* 0.15 = 3), but will allow more for a
spacer of length 30 (30 \* 0.15 = 4.5). We add a paragraph in the
methods section clarifying this choice compared to the default
parameters and parameters used by certain tools and previous studies
(e.g. the original IMG/VR4 spacer-protospacer mapping used default
blastn with -max\_target\_seqs 1000). Regarding "CRISPRTarget" -
the associated website (and link provided) is not accessible
(ERR\_CONNECTION\_TIMED\_OUT). The research paper describing
CRISPRTarget (Biswas et. al. 2013) mentions that the website allows
modifying (selecting custom values for word size, evalue and choice of
target database) for the BLASTn step, and using the defaults
blastn-short values for word size (no mention of max-target-seqs), with
a permissive e-value (to account for variations in database size) and
modified gap penalties: "The initial CRISPRTarget defaults are the same
except that a gap is penalized more highly (-10), the mismatch penalty
is -1 and the E filter is 1. In addition, there is also no filter or
masking for low complexity". These adjustments should favour gapless
alignments compared to blastn-short defaults (gap open -5, gap extend
-2). While we do not disagree with these parameter choices, we prefer to
use the qcov and %id options to control the reported alignments. We
added a paragraph addressing this in the main text. Note - we do not see
mention of using a modified "-max\_target\_seqs" in CRISPRTarget, and
thus must assume the default (500) is used.

**3. DUST Filtering / Low Complexity Analysis**

The reviewer suggests testing "the effect of DUST" and asks "are the
sequences with high abundance being missed of low complexity."

**Response:** We have implemented a light pre-filtering step
(spacer_inspection.ipynb) using several metrics to pre-emptively measure
and discard low-complexity, repetitive, or homopolymer spacers prior to
any mapping/alignment/benchmarking. We added a paragraph regarding this
to the methods section and now explicitly note this as a key
recommendation for future analyses.

We further note that while some tools employ internal complexity filtering
(e.g., DUST in blastn, tantan in mmseqs2), not all tools do so consistently
or allow us to disable these filters. Therefore, applying pre-filtering
uniformly across all tools ensures comparable input handling.

Regarding whether the high-abundance spacers being missed are low-complexity:
we note that high-abundance spacers (>1000 occurrences) are still observed
in the reanalysis of this version, indicating they pass the complexity filters.

**4. Tool Recommendations and Guidelines**

The reviewer notes that "the authors state they provide general
guidelines in the abstract but these are less clear in the discussion."

**Response:** We have clarified and expanded our tool recommendations
with specific use cases (when and why to use certain tools) and have
added a flowchart for tool selection under different scenarios (see figure nnn).

-   **Primary recommendation - Bowtie1 (hamming distance ≤3):** For most
    large-scale metagenomic analyses, we recommend bowtie1 with hamming
    distance ≤3. This provides high recall (>99% for spacers with 0-3
    substitutions), a distance threshold that we find biologically
    appropriate. Bowtie1 offers good computational efficiency and
    scalability even for databases with millions of spacers and hundreds
    of thousands of contigs.
-   **Extended hamming distance - Indelfree.sh indexed mode (hamming
    distance ≤5):** For scenarios requiring detection of spacers with up
    to 5 substitutions, we recommend indelfree.sh in indexed mode. While
    more computationally intensive than bowtie1, it maintains near-perfect
    recall and avoids the dramatic false positive increase associated
    with using edit distance. It could be further optimized if the
    minimum spacer length is known in advance to be ≥21 (allowing a
    k-mer size of 8 instead of 7). Note that the false positive rate
    increases beyond 3 substitutions, so this tool is best suited for
    applications where extended mismatch tolerance is necessary and some
    increase in false positives is acceptable.
-   **Sassy (edit distance) - Specific use cases:** Sassy is the only
    tool evaluated with perfect recall that supports arbitrary edit
    distances. However, it incurs massive computational costs (e.g.,
    \~1M CPU seconds for 5% subsample with ≤5 edits). We recommend sassy
    only for:

    -   **Small datasets** where computational cost is acceptable
    -   **Comparative studies** of escape mutation types in controlled
        (experimental) phage-host systems. We note that such studies are
        currently lacking, as the literature seems focused on measuring
        substitution number and location rather than mutation type
        (substitution vs. indel vs. rearrangement) and frequency. This
        positions sassy as the only current recommendation for such
        studies, which aligns with the historically small-scale nature of
        such experimental setups (e.g., several phages and several host
        strains).
    -   **Testing with ≥3 substitutions and indels:** It is possible that
        certain phage mutations accumulate both indels and substitutions,
        particularly when the actual host of a divergent phage differs
        from the closest available host from which spacers were extracted.

    -   **Low-accuracy raw long reads** (e.g., Oxford Nanopore R9
        chemistry) where indels are a larger concern and some minimal
        edit distance tolerance might be considered
    -   **Methodological validation** to establish baseline performance
        when developing new tools, or deciding on benchmarking
        guidelines (such as in this work)

-   **Sassy limitations:** We explicitly note that sassy's computational
    requirements (1M CPU seconds for just 5% of data) make it prohibitive
    for large-scale metagenomic analyses. This is not a limitation of
    sassy as a tool, but rather a fundamental limitation of using edit
    distance with exhaustive (non-heuristic) search for this biological
    question at large scale. Notably, sassy actually has lower runtime
    than some of the evaluated exhaustive tools that only report hamming
    distances, which is quite impressive.
-   **Use case scenarios:** We provide specific recommendations:

    -   **High-throughput host-virus prediction:** bowtie1 with hamming
        ≤3
    -   **Extended mismatch tolerance:** indelfree.sh indexed with
        hamming ≤5
    -   **Small experimental datasets:** sassy with appropriate distance
        threshold
    -   **Precision-focused applications:** Consider DUST filtering and
        conservative thresholds

-   **Resource trade-offs:** We discuss computational requirements
    comprehensively, showing that bowtie1 and indelfree.sh indexed mode
    scale well to large databases, while sassy is limited to small-scale
    applications.
-   **Parameter guidance:** We provide specific parameter
    recommendations for different tools and discuss how parameter
    choices affect performance, particularly for BLASTn where
    -max\_target\_seqs significantly impacts high-abundance spacer
    detection.

### Minor Points

**6. Figure Improvements**

We have addressed all figure issues:

-   **Line distinction:** Improved line styles and colors for better
    distinction
-   **Figure 3 readability:** Implemented log scaling and better color
    schemes
-   **Supplementary figures:** Fixed typos and improved legends
-   **Color consistency:** Colors and symbols representing tools should
    now be consistent across figures.

**7. Historical Comparison**

The reviewer notes: "which would 'allow for a direct comparison to
historical results'. Was that comparison ever made? I could not find any
further mention to that in the manuscript."

**Response:** Our original choice to benchmark tools using the IMG/VR4
contig dataset is inherently a comparison with the original IMG/VR4
publication, which used spacer-protospacer matching for host assignment
(using BLASTn with -max\_target\_seqs 1000). In the current
version, we removed the explicit mention of comparison to historical
results, as we are now using a different spacer sequence set (iPhoP June
2025).

## Supplementary Material
All new notebooks, datasets, and analysis scripts are provided in the
public github repository. The final mapping results (i.e. files too
large for github) are uploaded as a new version to the same zenodo
deposit.
