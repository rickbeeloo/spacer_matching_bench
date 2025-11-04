use indicatif::{ProgressBar, ProgressStyle};
use pyo3::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use bio::io::fasta;
use std::cmp::max;
use std::sync::Arc;
use rand_distr::{Normal, Distribution, Beta};

#[derive(Debug, Clone)]
#[pyclass]
pub enum DistributionType {
    Uniform,
    Normal,
    Bell, // Beta distribution for bell curve
}

#[pymethods]
impl DistributionType {
    #[new]
    fn new(dist_type: &str) -> PyResult<Self> {
        match dist_type.to_lowercase().as_str() {
            "uniform" => Ok(DistributionType::Uniform),
            "normal" => Ok(DistributionType::Normal),
            "bell" => Ok(DistributionType::Bell),
            _ => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Distribution type must be 'uniform', 'normal', or 'bell'"
            ))
        }
    }
}

#[derive(Debug, Clone)]
#[pyclass]
pub struct BaseComposition {
    pub a_frac: f64,
    pub t_frac: f64,
    pub c_frac: f64,
    pub g_frac: f64,
}

impl Default for BaseComposition {
    fn default() -> Self {
        BaseComposition {
            a_frac: 0.25,
            t_frac: 0.25,
            c_frac: 0.25,
            g_frac: 0.25,
        }
    }
}

#[pymethods]
impl BaseComposition {
    #[new]
    fn new(a_frac: f64, t_frac: f64, c_frac: f64, g_frac: f64) -> Self {
        let total = a_frac + t_frac + c_frac + g_frac;
        BaseComposition {
            a_frac: a_frac / total,
            t_frac: t_frac / total,
            c_frac: c_frac / total,
            g_frac: g_frac / total,
        }
    }
    
    #[staticmethod]
    fn from_gc_content(gc_content: f64) -> Self {
        let gc_frac = gc_content / 100.0;
        let at_frac = (1.0 - gc_frac) / 2.0;
        BaseComposition {
            a_frac: at_frac,
            t_frac: at_frac,
            c_frac: gc_frac / 2.0,
            g_frac: gc_frac / 2.0,
        }
    }
    
    #[staticmethod]
    fn from_fractions(a: f64, t: f64, c: f64, g: f64) -> Self {
        let total = a + t + c + g;
        BaseComposition {
            a_frac: a / total,
            t_frac: t / total,
            c_frac: c / total,
            g_frac: g / total,
        }
    }
    
    fn get_weights(&self) -> Vec<f64> {
        vec![self.a_frac, self.t_frac, self.c_frac, self.g_frac]
    }
}


#[pyclass]
struct Simulator {
    nucleotides: Vec<char>,
    base_composition: BaseComposition,
}

#[pymethods]
impl Simulator {
    #[new]
    fn new() -> Self {
        Simulator {
            nucleotides: vec!['A', 'T', 'C', 'G'],
            base_composition: BaseComposition::default(),
        }
    }
    
    // Added method to set base composition (Uri Neri @ 20.10.2025)
    fn set_base_composition(&mut self, composition: BaseComposition) {
        self.base_composition = composition;
    }
    
    // Added method to set base composition from GC content
    fn set_gc_content(&mut self, gc_content: f64) {
        self.base_composition = BaseComposition::from_gc_content(gc_content);
    }
    
    // Added methods to set base composition from fractions and GC content
    fn set_base_fractions(&mut self, a: f64, t: f64, c: f64, g: f64) {
        self.base_composition = BaseComposition::from_fractions(a, t, c, g);
    }

    // modified
    fn generate_random_sequence(&self, length: usize) -> String {
        self.generate_random_sequence_with_composition(length, &self.base_composition)
    }
    
    fn generate_random_sequence_with_composition(&self, length: usize, composition: &BaseComposition) -> String {
        let mut rng = rand::rng();
        let weights = composition.get_weights();
        
        (0..length)
            .map(|_| {
                let r: f64 = rng.random();
                let mut cumulative = 0.0;
                for (i, &weight) in weights.iter().enumerate() {
                    cumulative += weight;
                    if r <= cumulative {
                        return self.nucleotides[i];
                    }
                }
                // Fallback (shouldn't happen with proper weights)
                self.nucleotides[0]
            })
            .collect()
    }

    // Added to allow using non-uniform distributions for contig and / or spacer lengths (default set to uniform for backward consistency)
    fn generate_length_from_distribution(
        &self,
        range: (usize, usize),
        distribution_type: &DistributionType,
    ) -> usize {
        let min_len = range.0;
        let max_len = range.1;
        
        match distribution_type {
            DistributionType::Uniform => {
                let mut rng = rand::rng();
                rng.random_range(min_len..=max_len)
            }
            DistributionType::Normal => {
                let mean = (min_len + max_len) as f64 / 2.0;
                let std_dev = (max_len - min_len) as f64 / 6.0; // 99.7% within range
                let normal = Normal::new(mean, std_dev).unwrap();
                
                loop {
                    let mut rng = rand::rng();
                    let sample = normal.sample(&mut rng) as usize;
                    if sample >= min_len && sample <= max_len {
                        return sample;
                    }
                }
            }
            DistributionType::Bell => {
                // Use Beta distribution for bell curve
                let alpha = 2.0;
                let beta = 2.0;
                let beta_dist = Beta::new(alpha, beta).unwrap();
                let mut rng = rand::rng();
                let sample = beta_dist.sample(&mut rng);
                min_len + ((sample * (max_len - min_len) as f64) as usize)
            }
        }
    }

        // Add these methods to the Simulator impl
    fn set_contig_distribution(&mut self, dist_type: &str) -> PyResult<()> {
        match dist_type.to_lowercase().as_str() {
            "uniform" => Ok(()),
            "normal" => Ok(()),
            "bell" => Ok(()),
            _ => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Distribution type must be 'uniform', 'normal', or 'bell'"
            ))
        }
    }

    fn set_spacer_distribution(&mut self, dist_type: &str) -> PyResult<()> {
        match dist_type.to_lowercase().as_str() {
            "uniform" => Ok(()),
            "normal" => Ok(()),
            "bell" => Ok(()),
            _ => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Distribution type must be 'uniform', 'normal', or 'bell'"
            ))
        }
    }

    /// Read FASTA file and return HashMap of id -> sequence
    fn read_fasta_file(&self, filepath: &str) -> PyResult<HashMap<String, String>> {
        let file = match File::open(filepath) {
            Ok(f) => f,
            Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Could not open file '{}': {}", filepath, e)
            ))
        };
        
        let reader = fasta::Reader::new(file);
        let mut sequences = HashMap::new();
        
        for result in reader.records() {
            let record = match result {
                Ok(r) => r,
                Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Error reading FASTA record: {}", e)
                ))
            };
            
            let id = record.id().to_string();
            let seq = String::from_utf8_lossy(record.seq()).to_string();
            sequences.insert(id, seq);
        }
        
        Ok(sequences)
    }

    /// Apply mutations to a DNA sequence.
    ///
    /// This function modifies a spacer sequence by applying three types of mutations:
    /// 1. Insertion indels: Add random bases at random positions
    /// 2. Deletion indels: Remove bases at random positions  
    /// 3. Substitution mismatches: Replace bases with different bases
    ///
    /// # Parameters
    /// - `sequence`: The original DNA sequence to mutate
    /// - `n_mismatches`: Number of substitution mutations to apply
    /// - `n_insertions`: Number of insertion mutations (bases to ADD)
    /// - `n_deletions`: Number of deletion mutations (bases to REMOVE)
    ///
    /// # Returns
    /// The mutated sequence as a String
    ///
    /// # Important Notes
    /// - Mutations are applied in order: insertions → deletions → mismatches
    /// - This order ensures the final sequence length is: original + insertions - deletions
    /// - The mutated sequence is what gets inserted into contigs
    /// - The original sequence is preserved in simulated_spacers.fa
    fn apply_mutations(&self, sequence: &str, n_mismatches: usize, n_insertions: usize, n_deletions: usize) -> String {
        let mut rng = rand::rng();
        let mut sequence: Vec<char> = sequence.chars().collect();
        
        // STEP 1: Apply insertion indels
        // Each insertion adds one random base at a random position
        for _ in 0..n_insertions {
            let insert_pos = rng.random_range(0..=sequence.len());
            let new_base = *self.nucleotides.choose(&mut rng).unwrap();
            sequence.insert(insert_pos, new_base);
        }
        
        // STEP 2: Apply deletion indels
        // Each deletion removes one base at a random position
        for _ in 0..n_deletions {
            if !sequence.is_empty() {
                let delete_pos = rng.random_range(0..sequence.len());
                sequence.remove(delete_pos);
            }
        }
        
        // STEP 3: Apply substitution mismatches
        // Replace random bases with different bases (ensures no A→A substitutions)
        if n_mismatches > 0 && !sequence.is_empty() {
            let mismatch_positions: HashSet<usize> = (0..sequence.len())
                .collect::<Vec<_>>()
                .choose_multiple(&mut rng, n_mismatches)
                .cloned()
                .collect();

            for pos in mismatch_positions {
                let original = sequence[pos];
                let valid_bases: Vec<_> = self
                    .nucleotides
                    .iter()
                    .filter(|&&b| b != original)
                    .collect();
                let new_base = *valid_bases.choose(&mut rng).unwrap();
                sequence[pos] = *new_base;
            }
        }

        sequence.into_iter().collect()
    }

    fn reverse_complement(&self, sequence: &str) -> String {
        let complement: HashMap<char, char> = [
            ('A', 'T'),
            ('T', 'A'),
            ('C', 'G'),
            ('G', 'C'),
            ('N', 'N'),
            ('R', 'Y'),
            ('Y', 'R'),
            ('W', 'W'),
            ('S', 'S'),
            ('M', 'K'),
            ('K', 'M'),
            ('B', 'V'),
            ('V', 'B'),
            ('D', 'H'),
            ('H', 'D'),
        ]
        .iter()
        .cloned()
        .collect();

        sequence
            .chars()
            .rev()
            .map(|base| *complement.get(&base).unwrap_or(&base))
            .collect()
    }

    // Private method (not exposed to Python)
    fn verify_simulation(
        &self,
        contigs: HashMap<String, String>,
        spacers: HashMap<String, String>,
        ground_truth: Vec<Vec<String>>,
    ) -> bool {
        
        for entry in ground_truth {
            let spacer_id = &entry[0];
            let contig_id = &entry[1];
            let start_pos: usize = entry[2].parse().unwrap();
            let end_pos: usize = entry[3].parse().unwrap();
            let is_rc = entry[4] == "true";
            let expected_mismatches: usize = entry[5].parse().unwrap();

            // Get the original spacer and contig sequences
            let spacer = match spacers.get(spacer_id) {
                Some(s) => s,
                None => {
                    println!("Error: Spacer {} not found", spacer_id);
                    return false;
                }
            };
            let contig = match contigs.get(contig_id) {
                Some(c) => c,
                None => {
                    println!("Error: Contig {} not found", contig_id);
                    return false;
                }
            };

            // Extract the region from the contig
            if start_pos >= contig.len() || end_pos > contig.len() {
                println!("Error: Invalid coordinates for contig {}: start={}, end={}, len={}", 
                    contig_id, start_pos, end_pos, contig.len());
                return false;
            }
            
            let contig_region = &contig[start_pos..end_pos];

            // Compare sequences character by character
            let mut actual_mismatches = 0;
            if is_rc {
                // For reverse complement, we need to compare with the reverse complement of the spacer
                let rc_spacer = self.reverse_complement(spacer);
                for (c1, c2) in rc_spacer.chars().zip(contig_region.chars()) {
                    if c1 != c2 {
                        actual_mismatches += 1;
                    }
                }
            } else {
                // For forward strand, compare directly
                for (c1, c2) in spacer.chars().zip(contig_region.chars()) {
                    if c1 != c2 {
                        actual_mismatches += 1;
                    }
                }
            }

            // Check if the number of mismatches matches
            if actual_mismatches != expected_mismatches {
                println!("Verification failed for spacer {} in contig {}:", spacer_id, contig_id);
                println!("Expected mismatches: {}, Actual mismatches: {}", expected_mismatches, actual_mismatches);
                println!("Spacer: {}", spacer);
                println!("Contig region: {}", contig_region);
                println!("Is RC: {}", is_rc);
                return false;
            }
        }
        true
    }

    fn simulate_data(
        &self,
        contig_length_range: (usize, usize),
        spacer_length_range: (usize, usize),
        n_mismatch_range: (usize, usize),
        mut sample_size_contigs: usize,
        sample_size_spacers: usize,
        insertion_range: (usize, usize),
        n_insertion_range: (usize, usize),
        n_deletion_range: (usize, usize),
        prop_rc: f64,
        threads: usize,
        verify: bool,
        output_dir: String,
        id_prefix: Option<String>,
        // New optional parameters with defaults
        contig_distribution: Option<DistributionType>,
        spacer_distribution: Option<DistributionType>,
        base_composition: Option<BaseComposition>,
        data_subdir: Option<String>,
    ) -> PyResult<(HashMap<String, String>, HashMap<String, String>, Vec<Vec<String>>)> {
        
        // Estimate required contig space
        let avg_spacer_length = (spacer_length_range.0 + spacer_length_range.1) / 2;
        // For a more conservative estimate, use max insertions if the range is wide
        let insertion_ratio = insertion_range.1 as f64 / insertion_range.0 as f64;
        let avg_insertions_per_spacer = if insertion_ratio > 3.0 {
            // Use a weighted average leaning towards the maximum for wide ranges
            (insertion_range.0 as f64 * 0.3 + insertion_range.1 as f64 * 0.7) as usize
        } else {
            (insertion_range.0 + insertion_range.1) / 2
        };
        
        let total_expected_insertion_space = sample_size_spacers * avg_spacer_length * avg_insertions_per_spacer;
        
        let avg_contig_length = (contig_length_range.0 + contig_length_range.1) / 2;
        let expected_total_contig_space = sample_size_contigs * avg_contig_length;
        
        // Calculate the target number of contigs needed
        // Use a safety factor to ensure enough space (50% more than minimum)
        let safety_factor = 1.5;
        let target_contigs = max(
            sample_size_contigs, 
            (total_expected_insertion_space as f64 * safety_factor / avg_contig_length as f64).ceil() as usize
        );
        
        if target_contigs > sample_size_contigs {
            println!("Warning: Estimated required contig space ({}) exceeds available space ({})", 
                total_expected_insertion_space, expected_total_contig_space);
            println!("Increasing number of contigs from {} to {} to accommodate insertions", 
                sample_size_contigs, target_contigs);
            sample_size_contigs = target_contigs;
        }
        
        // Instead of building the global thread pool, build a local one.
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
    
        // All parallel work is now done inside pool.install
        pool.install(|| {
            println!("Generating contigs...");
            let pb = ProgressBar::new(sample_size_contigs as u64);
            pb.set_style(ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap());
    
            // Generate contigs in parallel
            let contigs: HashMap<String, String> = (0..sample_size_contigs)
                .into_par_iter()
                .map(|i| {
                    let mut rng = rand::rng();
                    let length = if let Some(ref dist) = contig_distribution {
                        self.generate_length_from_distribution(contig_length_range, dist)
                    } else {
                        rng.random_range(contig_length_range.0..=contig_length_range.1)
                    };
                    let sequence = if let Some(ref comp) = base_composition {
                        self.generate_random_sequence_with_composition(length, comp)
                    } else {
                        self.generate_random_sequence(length)
                    };
                    pb.inc(1);
                    
                    // Apply prefix if provided
                    let id = match &id_prefix {
                        Some(prefix) => format!("{}_contig_{}", prefix, i),
                        None => format!("contig_{}", i),
                    };
                    
                    (id, sequence)
                })
                .collect();
            pb.finish_with_message("Contigs generated");
    
            println!("Generating spacers...");
            let pb = ProgressBar::new(sample_size_spacers as u64);
            pb.set_style(ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap());
    
            // Generate spacers in parallel
            let spacers: HashMap<String, String> = (0..sample_size_spacers)
                .into_par_iter()
                .map(|i| {
                    let mut rng = rand::rng();
                    let length = if let Some(ref dist) = spacer_distribution {
                        self.generate_length_from_distribution(spacer_length_range, dist)
                    } else {
                        rng.random_range(spacer_length_range.0..=spacer_length_range.1)
                    };
                    let sequence = if let Some(ref comp) = base_composition {
                        self.generate_random_sequence_with_composition(length, comp)
                    } else {
                        self.generate_random_sequence(length)
                    };
                    pb.inc(1);
                    
                    // Apply prefix if provided
                    let id = match &id_prefix {
                        Some(prefix) => format!("{}_spacer_{}", prefix, i),
                        None => format!("spacer_{}", i),
                    };
                    
                    (id, sequence)
                })
                .collect();
            pb.finish_with_message("Spacers generated");
    
            // Predetermine insertion parameters for each spacer
            println!("Predetermining spacer insertion parameters...");
            
            // Structure to hold insertion parameters for each spacer
            #[derive(Debug, Clone)]
            struct SpacerInsertionPlan {
                spacer_id: String,
                n_insertions: usize,
                total_length: usize,
                insertion_plans: Vec<(bool, usize, usize, usize)>, // (is_rc, n_mismatches, n_insertions, n_deletions) tuples
            }
            
            let mut insertion_plans: Vec<SpacerInsertionPlan> = Vec::with_capacity(spacers.len());
            let mut total_insertion_length = 0;
            let mut rng = rand::rng();
            
            for (id, seq) in &spacers {
                // ========================================================================
                // CRITICAL: Understanding the two types of "insertion"
                // ========================================================================
                // 
                // 1. SPACER INSERTION (n_insertions below):
                //    - How many TIMES this spacer gets placed into contigs
                //    - Controlled by `insertion_range` parameter
                //    - Example: n_insertions=2 means this spacer appears in 2 different locations
                //    - This creates the ground truth data we're trying to find
                //
                // 2. INSERTION MUTATION (n_insertions in the loop below):
                //    - How many insertion INDELS (added bases) within the spacer sequence
                //    - Controlled by `n_insertion_range` parameter  
                //    - Example: n_insertions=1 means add 1 extra base to the spacer
                //    - This is a sequence-level mutation that makes alignment harder
                // (these were clear to us but in the paper and in a few other places, we call these "occurrences". But as occurrences is passive, and here we are explicitly creating them, I use "insertions" instead.)
                // Similarly, `n_deletions` refers to deletion MUTATIONS (removed bases)
                // ========================================================================
                
                // Determine number of TIMES to insert this spacer (placement count)
                let n_insertions = rng.random_range(insertion_range.0..=insertion_range.1);
                
                // Create insertion plans - one plan per spacer placement
                let mut plans = Vec::with_capacity(n_insertions);
                for _ in 0..n_insertions {
                    // For each placement, determine:
                    let is_rc = rng.random_bool(prop_rc);  // Reverse complement
                    let n_mismatches = rng.random_range(n_mismatch_range.0..=n_mismatch_range.1);  // Substitutions
                    let n_insertions = rng.random_range(n_insertion_range.0..=n_insertion_range.1);  // Insertion MUTATIONS
                    let n_deletions = rng.random_range(n_deletion_range.0..=n_deletion_range.1);    // Deletion MUTATIONS
                    plans.push((is_rc, n_mismatches, n_insertions, n_deletions));
                }
                
                let total_length = seq.len() * n_insertions;
                total_insertion_length += total_length;
                
                insertion_plans.push(SpacerInsertionPlan {
                    spacer_id: id.clone(),
                    n_insertions,
                    total_length,
                    insertion_plans: plans,
                });
            }
            
            // Sort spacers by total length (descending)
            insertion_plans.sort_by(|a, b| b.total_length.cmp(&a.total_length));
            
            println!("Total planned insertion length: {} bp", total_insertion_length);
            println!("Average insertions per spacer: {:.1}", 
                     insertion_plans.iter().map(|p| p.n_insertions).sum::<usize>() as f64 / insertion_plans.len() as f64);
            
            // Create sorted list of contigs by length (descending)
            let mut sorted_contigs: Vec<(String, String)> = contigs.iter()
                .map(|(id, seq)| (id.clone(), seq.clone()))
                .collect();
            sorted_contigs.sort_by(|a, b| b.1.len().cmp(&a.1.len()));
            
            // Calculate total contig size to verify we have enough space
            let total_contig_size: usize = sorted_contigs.iter().map(|(_, seq)| seq.len()).sum();
            
            // Calculate utilization percentage - how much of the available contig space would be used
            let utilization_percentage = (total_insertion_length as f64 / total_contig_size as f64) * 100.0;
            
            println!("Total contig size: {} bp", total_contig_size);
            println!("Expected contig utilization: {:.1}%", utilization_percentage);
            
            // Different warning levels based on utilization percentage
            if utilization_percentage > 80.0 {
                println!("WARNING: Very high contig utilization expected (>80%)");
                println!("This simulation will likely fail to insert all spacers at the expected rate.");
                println!("Consider increasing contig count, contig size, or reducing spacer insertions.");
                
                // If utilization is extremely high, we could even stop the simulation
                if utilization_percentage > 100.0 {
                    return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                        format!("Expected spacer insertion length ({} bp) exceeds total contig size ({} bp). Increase contig size or count, or reduce spacer insertions.", 
                            total_insertion_length, total_contig_size)
                    ));
                }
            } else if utilization_percentage > 60.0 {
                println!("NOTICE: Moderate to high contig utilization expected (>60%)");
                println!("Some spacers may not be inserted at the expected rate.");
            }
    
            // Limit threads to ensure we don't have more threads than contigs
            let effective_threads = std::cmp::min(threads, sorted_contigs.len());
            if effective_threads < threads {
                println!("WARNING: Reducing number of threads from {} to {} to match available contigs", 
                         threads, effective_threads);
            }

            // Distribute spacers and contigs to threads more intelligently
            let mut thread_spacer_plans: Vec<Vec<SpacerInsertionPlan>> = vec![Vec::new(); effective_threads];
            let mut thread_contigs: Vec<HashMap<String, String>> = vec![HashMap::new(); effective_threads];
            
            // Calculate total insertion length per thread (target)
            let _target_length_per_thread = total_insertion_length / effective_threads;
            
            // First, distribute spacers to balance total insertion length
            let mut thread_lengths = vec![0; effective_threads];
            
            for plan in &insertion_plans {
                // Find the thread with the least total length
                let min_thread = thread_lengths.iter()
                    .enumerate()
                    .min_by_key(|(_, &len)| len)
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);
                
                thread_spacer_plans[min_thread].push(plan.clone());
                thread_lengths[min_thread] += plan.total_length;
            }
            
            // First, ensure all threads with spacers get at least one contig
            let mut thread_indices: Vec<usize> = (0..effective_threads)
                .filter(|&i| !thread_spacer_plans[i].is_empty())
                .collect();

            // Sort threads by workload (highest first)
            thread_indices.sort_by(|&a, &b| thread_lengths[b].cmp(&thread_lengths[a]));

            // Simpler approach: distribute contigs directly to threads with spacers
            let mut sorted_contigs_vec = sorted_contigs;

            // First, make sure every thread with spacers gets at least one contig
            for &thread_idx in &thread_indices {
                if sorted_contigs_vec.is_empty() {
                    break; // No more contigs to distribute
                }
                
                // Give this thread one contig
                let (id, seq) = sorted_contigs_vec.remove(0);
                thread_contigs[thread_idx].insert(id, seq);
            }

            // Keep distributing contigs in the same order until we run out
            let mut i = 0;
            while !sorted_contigs_vec.is_empty() {
                let thread_idx = thread_indices[i % thread_indices.len()];
                let (id, seq) = sorted_contigs_vec.remove(0);
                thread_contigs[thread_idx].insert(id, seq);
                i += 1;
            }
            
            // Print thread workload distribution
            for i in 0..effective_threads {
                let spacer_count = thread_spacer_plans[i].len();
                let contig_count = thread_contigs[i].len();
                let spacer_length: usize = thread_spacer_plans[i].iter().map(|p| p.total_length).sum();
                let contig_length: usize = thread_contigs[i].values().map(|s| s.len()).sum();
                
                println!("Thread {}: {} spacers ({} bp planned) / {} contigs ({} bp)",
                         i, spacer_count, spacer_length, contig_count, contig_length);
            }
    
            // Calculate total insertion attempts 
            let total_insertions: u64 = insertion_plans.iter()
                .map(|p| p.insertion_plans.len())
                .sum::<usize>() as u64;

            // Create a single progress bar for the process
            println!("Processing spacer insertions...");
            let pb = ProgressBar::new(total_insertions);
            pb.set_style(ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap());

            // Wrap it in Arc to share between threads
            let pb = Arc::new(pb);

            // Pass cloned reference to each thread
            let results: Vec<_> = thread_contigs.into_par_iter()
                .enumerate()
                .map(|(thread_idx, mut group_contigs)| {
                    // Clone the progress bar reference for this thread
                    let thread_pb = pb.clone();
                    
                    let mut local_ground_truth: Vec<Vec<String>> = Vec::new();
                    
                    // Instead of tracking used positions, directly track available ranges
                    let mut contig_available_ranges: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
                    
                    // Initialize available ranges for each contig
                    for (contig_id, contig) in &group_contigs {
                        // Each contig starts with one available range from 0 to length
                        contig_available_ranges.insert(contig_id.clone(), vec![(0, contig.len())]);
                    }

                    // Process spacers assigned to this thread
                    let thread_plans = &thread_spacer_plans[thread_idx];
                    
                    // Count total insertion attempts (not just spacers)
                    let total_insertions: u64 = thread_plans.iter()
                        .map(|p| p.insertion_plans.len())
                        .sum::<usize>() as u64;

                    // Create progress bar with correct count
                    let _pb = ProgressBar::new(total_insertions);
                    
                    for plan in thread_plans {
                        let spacer_id = &plan.spacer_id;
                        let spacer = spacers.get(spacer_id).unwrap();
                        let spacer_len = spacer.len();
                        
                        // Cache reverse complement and mismatched versions to avoid repeated calculations
                        let mut cached_variants: HashMap<(bool, usize, usize, usize), String> = HashMap::new();
                        
                        // Process each insertion plan for this spacer
                        for &(is_rc, n_mismatches, n_insertions, n_deletions) in &plan.insertion_plans {
                            // Get or calculate the variant
                            let spacer_variant = cached_variants.entry((is_rc, n_mismatches, n_insertions, n_deletions))
                                .or_insert_with(|| {
                                    let base = if is_rc {
                                        self.reverse_complement(spacer)
                                    } else {
                                        spacer.to_string()
                                    };
                                    self.apply_mutations(&base, n_mismatches, n_insertions, n_deletions)
                                });
                            
                            // Instead of generating all positions:
                            let mut available_ranges: Vec<(String, usize, usize)> = Vec::new();
                            let mut total_range_length = 0;

                            for (contig_id, ranges) in &contig_available_ranges {
                                for &(start, end) in ranges {
                                    if end - start >= spacer_len {
                                        let viable_length = end - start - spacer_len + 1;
                                        total_range_length += viable_length;
                                        available_ranges.push((contig_id.clone(), start, end));
                                    }
                                }
                            }

                            // If no ranges available, skip this insertion
                            if available_ranges.is_empty() {
                                continue;
                            }

                            // Choose a random range weighted by the number of possible insertion points
                            let mut rng = rand::rng();
                            let r = rng.random_range(0..total_range_length);
                            let mut cumulative_length = 0;
                            let mut selected_range = None;

                            for (contig_id, start, end) in &available_ranges {
                                let viable_length = end - start - spacer_len + 1;
                                cumulative_length += viable_length;
                                if r < cumulative_length {
                                    selected_range = Some((contig_id.clone(), *start, *end));
                                    break;
                                }
                            }

                            if let Some((target_contig_id, range_start, range_end)) = selected_range {
                                // Choose a random position within the selected range
                                let max_start = range_end - spacer_len;
                                let start_pos = rng.random_range(range_start..=max_start);
                                let end_pos = start_pos + spacer_len;
                                
                                // Update available ranges
                                let ranges = contig_available_ranges.get_mut(&target_contig_id).unwrap();
                                for i in 0..ranges.len() {
                                    let (r_start, r_end) = ranges[i];
                                    
                                    if start_pos >= r_start && end_pos <= r_end {
                                        // Remove the current range
                                        let _current_range = ranges.remove(i);
                                        
                                        // Add the ranges before and after the insertion, if they exist
                                        if start_pos > r_start {
                                            ranges.push((r_start, start_pos));
                                        }
                                        
                                        if end_pos < r_end {
                                            ranges.push((end_pos, r_end));
                                        }
                                        
                                        // No need to check other ranges
                                        break;
                                    }
                                }
                                
                                // Update contig sequence
                                let target_contig = group_contigs.get_mut(&target_contig_id).unwrap();
                                let mut chars: Vec<char> = target_contig.chars().collect();
                                
                                // Ensure we don't exceed the contig bounds
                                let max_end_pos = chars.len();
                                let splice_end_pos = std::cmp::min(end_pos, max_end_pos);
                                
                                let spacer_variant_chars: Vec<char> = spacer_variant.chars().collect();
                                chars.splice(start_pos..splice_end_pos, spacer_variant_chars);
                                *target_contig = chars.into_iter().collect();

                                // Record in ground truth - use the actual end position based on spacer_variant length
                                let actual_end_pos = start_pos + spacer_variant.len();
                                local_ground_truth.push(vec![
                                    spacer_id.to_string(),
                                    target_contig_id,
                                    start_pos.to_string(),
                                    actual_end_pos.to_string(),
                                    if is_rc { "true".to_string() } else { "false".to_string() },
                                    n_mismatches.to_string()
                                ]);
                            }
                            
                            // Increment the progress bar after attempting an insertion, 
                            // regardless of whether it was successful
                            thread_pb.inc(1);
                        }
                    }
                    
                    (group_contigs, local_ground_truth)
                })
                .collect();
    
            // After all threads are done, finish the progress bar
            pb.finish_with_message("Spacer insertions completed");
    
            // Combine results
            let mut final_contigs = HashMap::new();
            let mut final_ground_truth: Vec<Vec<String>> = Vec::new();
    
            for (group_contigs, group_ground_truth) in results {
                final_contigs.extend(group_contigs);
                final_ground_truth.extend(group_ground_truth);
            }

            println!("Successfully inserted {} spacer instances", final_ground_truth.len());

            // Create data subdirectory if it doesn't exist
            let data_subdir_name = data_subdir.unwrap_or_else(|| "simulated_data".to_string());
            let data_dir = format!("{}/{}", output_dir, data_subdir_name);
            if let Err(e) = std::fs::create_dir_all(&data_dir) {
                return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Could not create data directory '{}': {}", data_subdir_name, e)
                ));
            }
            
            // Write contigs to FASTA file
            println!("Writing contigs to FASTA file...");
            let contig_path = format!("{}/simulated_contigs.fa", data_dir);
            let contig_file = match File::create(&contig_path) {
                Ok(file) => file,
                Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Could not create contigs file: {}", e)
                )),
            };
            
            let mut contig_writer = fasta::Writer::new(contig_file);
            for (id, seq) in &final_contigs {
                if let Err(e) = contig_writer.write(id, None, seq.as_bytes()) {
                    return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                        format!("Error writing contig {}: {}", id, e)
                    ));
                }
            }

            // Write spacers to FASTA file
            println!("Writing spacers to FASTA file...");
            let spacer_path = format!("{}/simulated_spacers.fa", data_dir);
            let spacer_file = match File::create(&spacer_path) {
                Ok(file) => file,
                Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Could not create spacers file: {}", e)
                )),
            };
            
            let mut spacer_writer = fasta::Writer::new(spacer_file);
            for (id, seq) in &spacers {
                if let Err(e) = spacer_writer.write(id, None, seq.as_bytes()) {
                    return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                        format!("Error writing spacer {}: {}", id, e)
                    ));
                }
            }

            // Write planned ground truth to TSV file
            println!("Writing planned ground truth to TSV file...");
            let planned_ground_truth_path = format!("{}/planned_ground_truth.tsv", data_dir);
            let planned_ground_truth_file = match File::create(&planned_ground_truth_path) {
                Ok(file) => file,
                Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Could not create planned ground truth file: {}", e)
                )),
            };
            
            let mut planned_ground_truth_writer = BufWriter::new(planned_ground_truth_file);
            
            // Write header
            if let Err(e) = writeln!(planned_ground_truth_writer, "spacer_id\tcontig_id\tstart\tend\tstrand\tmismatches") {
                return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Error writing planned ground truth header: {}", e)
                ));
            }
            
            // Write data rows
            for row in &final_ground_truth {
                if row.len() != 6 {
                    return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                        format!("Invalid planned ground truth row length: {}", row.len())
                    ));
                }
                
                let line = format!("{}\t{}\t{}\t{}\t{}\t{}\n", 
                    row[0], row[1], row[2], row[3], 
                    row[4].to_lowercase(), row[5]);
                    
                if let Err(e) = write!(planned_ground_truth_writer, "{}", line) {
                    return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                        format!("Error writing planned ground truth line: {}", e)
                    ));
                }
            }

            // Verify the simulation if requested
            if verify {
                println!("Verifying simulation...");
                if !self.verify_simulation(final_contigs.clone(), spacers.clone(), final_ground_truth.clone()) {
                    return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                        "Simulated data verification failed"
                    ));
                }
                println!("Simulation verification passed");
            }

            Ok((final_contigs, spacers, final_ground_truth))
        })
    }

    /// Simulate data using pre-existing contigs and/or spacers from files or provided data
    /// This is a simpler simulation that only applies substitution mismatches (no indels)
    /// Similar to the Python simulate_data_from_input_files function
    fn simulate_from_input_files(
        &self,
        contig_length_range: (usize, usize),
        spacer_length_range: (usize, usize),
        n_mismatch_range: (usize, usize),
        sample_size_contigs: usize,
        sample_size_spacers: usize,
        insertion_range: (usize, usize),
        prop_rc: f64,
        output_dir: String,
        id_prefix: Option<String>,
        contigs_file: Option<String>,
        spacers_file: Option<String>,
        data_subdir: Option<String>,
        debug: bool,
    ) -> PyResult<(HashMap<String, String>, HashMap<String, String>, Vec<Vec<String>>)> {
        
        // Load or generate contigs
        let mut contigs: HashMap<String, String> = if let Some(file_path) = contigs_file {
            if debug {
                println!("Reading contigs from file: {}", file_path);
            }
            self.read_fasta_file(&file_path)?
        } else {
            if debug {
                println!("Generating contigs...");
            }
            let pb = ProgressBar::new(sample_size_contigs as u64);
            pb.set_style(ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap());
            
            let mut generated_contigs = HashMap::new();
            for i in 0..sample_size_contigs {
                let mut rng = rand::rng();
                let length = rng.random_range(contig_length_range.0..=contig_length_range.1);
                let sequence = self.generate_random_sequence(length);
                
                let id = match &id_prefix {
                    Some(prefix) => format!("{}_contig_{}", prefix, i),
                    None => format!("contig_{}", i),
                };
                
                generated_contigs.insert(id, sequence);
                pb.inc(1);
            }
            pb.finish_with_message("Contigs generated");
            generated_contigs
        };
        
        // Load or generate spacers
        let spacers: HashMap<String, String> = if let Some(file_path) = spacers_file {
            if debug {
                println!("Reading spacers from file: {}", file_path);
            }
            self.read_fasta_file(&file_path)?
        } else {
            if debug {
                println!("Generating spacers...");
            }
            let pb = ProgressBar::new(sample_size_spacers as u64);
            pb.set_style(ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap());
            
            let mut generated_spacers = HashMap::new();
            for i in 0..sample_size_spacers {
                let mut rng = rand::rng();
                let length = rng.random_range(spacer_length_range.0..=spacer_length_range.1);
                let sequence = self.generate_random_sequence(length);
                
                let id = match &id_prefix {
                    Some(prefix) => format!("{}_spacer_{}", prefix, i),
                    None => format!("spacer_{}", i),
                };
                
                generated_spacers.insert(id, sequence);
                pb.inc(1);
            }
            pb.finish_with_message("Spacers generated");
            generated_spacers
        };
        
        // Check for empty insertion range
        if insertion_range.0 == 0 && insertion_range.1 == 0 {
            if debug {
                println!("No insertions requested, returning empty ground truth");
            }
            return Ok((contigs, spacers, Vec::new()));
        }
        
        // Track used positions per contig to avoid overlaps
        let mut used_positions: HashMap<String, HashSet<usize>> = HashMap::new();
        for contig_id in contigs.keys() {
            used_positions.insert(contig_id.clone(), HashSet::new());
        }
        
        let mut ground_truth: Vec<Vec<String>> = Vec::new();
        let mut rng = rand::rng();
        
        if debug {
            println!("Starting spacer insertions...");
        }
        
        let pb = ProgressBar::new(spacers.len() as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap());
        
        // Process each spacer
        for (spacer_id, spacer_seq) in &spacers {
            let n_insertions = rng.random_range(insertion_range.0..=insertion_range.1);
            
            // For each insertion of this spacer
            for _ in 0..n_insertions {
                // Select random contig
                let contig_ids: Vec<&String> = contigs.keys().collect();
                let target_contig_id = contig_ids[rng.random_range(0..contig_ids.len())].clone();
                
                // Decide on reverse complement
                let is_rc = rng.random_bool(prop_rc);
                let spacer_to_insert = if is_rc {
                    self.reverse_complement(spacer_seq)
                } else {
                    spacer_seq.to_string()
                };
                
                // Apply mismatches (only substitutions, no indels)
                let n_mismatches = rng.random_range(n_mismatch_range.0..=n_mismatch_range.1);
                let mut spacer_with_errors = self.apply_mutations(&spacer_to_insert, n_mismatches, 0, 0);
                
                // Find non-overlapping position
                let target_contig = contigs.get(&target_contig_id).unwrap();
                let max_attempts = 1000;
                let mut attempt = 0;
                let mut position_found = false;
                let mut start_pos = 0;
                let mut end_pos = 0;
                
                while attempt < max_attempts {
                    let max_pos = target_contig.len().saturating_sub(spacer_with_errors.len());
                    if max_pos == 0 {
                        // Spacer is longer than contig, skip
                        if debug {
                            println!("Warning: spacer {} is longer than contig {}, skipping", 
                                spacer_id, target_contig_id);
                        }
                        break;
                    }
                    
                    start_pos = rng.random_range(0..=max_pos);
                    end_pos = start_pos + spacer_with_errors.len();
                    
                    // Check if this range overlaps with any used positions
                    let position_range: HashSet<usize> = (start_pos..end_pos).collect();
                    let used = used_positions.get(&target_contig_id).unwrap();
                    
                    if position_range.is_disjoint(used) {
                        // Found non-overlapping position
                        position_found = true;
                        break;
                    }
                    
                    attempt += 1;
                }
                
                if !position_found {
                    if debug {
                        println!("Could not find non-overlapping position for spacer {} in contig {} after {} attempts", 
                            spacer_id, target_contig_id, max_attempts);
                    }
                    continue;
                }
                
                // Mark positions as used
                for pos in start_pos..end_pos {
                    used_positions.get_mut(&target_contig_id).unwrap().insert(pos);
                }
                
                // Insert spacer into contig
                let contig = contigs.get_mut(&target_contig_id).unwrap();
                let mut contig_chars: Vec<char> = contig.chars().collect();
                let spacer_chars: Vec<char> = spacer_with_errors.chars().collect();
                
                // Replace the region
                for (i, &c) in spacer_chars.iter().enumerate() {
                    if start_pos + i < contig_chars.len() {
                        contig_chars[start_pos + i] = c;
                    }
                }
                *contig = contig_chars.into_iter().collect();
                
                // Record ground truth
                ground_truth.push(vec![
                    spacer_id.clone(),
                    target_contig_id.clone(),
                    start_pos.to_string(),
                    end_pos.to_string(),
                    if is_rc { "true".to_string() } else { "false".to_string() },
                    n_mismatches.to_string(),
                ]);
            }
            
            pb.inc(1);
        }
        pb.finish_with_message("Spacer insertions completed");
        
        println!("Successfully inserted {} spacer instances", ground_truth.len());
        
        // Create data subdirectory if it doesn't exist
        let data_subdir_name = data_subdir.unwrap_or_else(|| "simulated_data".to_string());
        let data_dir = format!("{}/{}", output_dir, data_subdir_name);
        if let Err(e) = std::fs::create_dir_all(&data_dir) {
            return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Could not create data directory '{}': {}", data_subdir_name, e)
            ));
        }
        
        // Write contigs to FASTA file
        println!("Writing contigs to FASTA file...");
        let contig_path = format!("{}/simulated_contigs.fa", data_dir);
        let contig_file = match File::create(&contig_path) {
            Ok(file) => file,
            Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Could not create contigs file: {}", e)
            )),
        };
        
        let mut contig_writer = fasta::Writer::new(contig_file);
        for (id, seq) in &contigs {
            if let Err(e) = contig_writer.write(id, None, seq.as_bytes()) {
                return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Error writing contig {}: {}", id, e)
                ));
            }
        }
        
        // Write spacers to FASTA file
        println!("Writing spacers to FASTA file...");
        let spacer_path = format!("{}/simulated_spacers.fa", data_dir);
        let spacer_file = match File::create(&spacer_path) {
            Ok(file) => file,
            Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Could not create spacers file: {}", e)
            )),
        };
        
        let mut spacer_writer = fasta::Writer::new(spacer_file);
        for (id, seq) in &spacers {
            if let Err(e) = spacer_writer.write(id, None, seq.as_bytes()) {
                return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Error writing spacer {}: {}", id, e)
                ));
            }
        }
        
        // Write ground truth to TSV file
        println!("Writing planned ground truth to TSV file...");
        let ground_truth_path = format!("{}/planned_ground_truth.tsv", data_dir);
        let ground_truth_file = match File::create(&ground_truth_path) {
            Ok(file) => file,
            Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Could not create ground truth file: {}", e)
            )),
        };
        
        let mut ground_truth_writer = BufWriter::new(ground_truth_file);
        
        // Write header
        if let Err(e) = writeln!(ground_truth_writer, "spacer_id\tcontig_id\tstart\tend\tstrand\tmismatches") {
            return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Error writing ground truth header: {}", e)
            ));
        }
        
        // Write data rows
        for row in &ground_truth {
            if row.len() != 6 {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    format!("Invalid ground truth row length: {}", row.len())
                ));
            }
            
            let line = format!("{}\t{}\t{}\t{}\t{}\t{}\n", 
                row[0], row[1], row[2], row[3], 
                row[4].to_lowercase(), row[5]);
                
            if let Err(e) = write!(ground_truth_writer, "{}", line) {
                return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("Error writing ground truth line: {}", e)
                ));
            }
        }
        
        Ok((contigs, spacers, ground_truth))
    }
}

#[pymodule]
fn rust_simulator(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Simulator>()?;
    m.add_class::<BaseComposition>()?;
    m.add_class::<DistributionType>()?;
    Ok(())
}
