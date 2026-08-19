use anyhow::{Context, Result};
use bio::io::fasta;
use clap::Parser;
use rayon::prelude::*;
use rust_annotale::{dna_to_rvds, is_dna, parse_rvd_sequence};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "TALE Family Builder using UPGMA and Length Mismatch Splitting",
    long_about = None
)]
struct Args {
    /// Input FASTA file (can be DNA sequences or RVD sequences)
    #[arg(short = 'i', long = "input", required = true)]
    input: PathBuf,

    /// Output directory for family assignments and summary
    #[arg(short = 'o', long = "outdir", default_value = ".")]
    outdir: PathBuf,

    /// UPGMA tree cut distance threshold
    #[arg(short = 'c', long = "cut", default_value_t = 6.0)]
    cut: f32,
}

// Custom RVD mismatch cost function
fn rvd_mismatch_cost(rvd1: &str, rvd2: &str) -> f32 {
    if rvd1 == rvd2 {
        return 0.0;
    }
    let chars1: Vec<char> = rvd1.chars().collect();
    let chars2: Vec<char> = rvd2.chars().collect();
    
    let c12_1 = chars1.get(0).copied().unwrap_or('-');
    let c13_1 = chars1.get(1).copied().unwrap_or('-');
    let c12_2 = chars2.get(0).copied().unwrap_or('-');
    let c13_2 = chars2.get(1).copied().unwrap_or('-');
    
    let mut cost = 0.0;
    if c12_1 != c12_2 {
        cost += 0.2;
    }
    if c13_1 != c13_2 {
        cost += 0.8;
    }
    cost
}

// Glocal dynamic programming alignment function
fn glocal_distance(seq1: &[String], seq2: &[String]) -> f32 {
    if seq1.is_empty() && seq2.is_empty() {
        return 0.0;
    }
    if seq1.is_empty() {
        return 5.0 * seq2.len() as f32;
    }
    if seq2.is_empty() {
        return 5.0 * seq1.len() as f32;
    }
    
    let (seq_a, seq_b) = if seq1.len() >= seq2.len() {
        (seq1, seq2)
    } else {
        (seq2, seq1)
    };
    
    let m = seq_a.len();
    let n = seq_b.len();
    
    let mut d = vec![vec![0.0; n + 1]; m + 1];
    for j in 1..=n {
        d[0][j] = 5.0 * j as f32;
    }
    for i in 1..=m {
        d[i][0] = 0.0;
    }
    
    for i in 1..=m {
        for j in 1..=n {
            let cost = rvd_mismatch_cost(&seq_a[i - 1], &seq_b[j - 1]);
            let match_score = d[i - 1][j - 1] + cost;
            let gap_a = d[i - 1][j] + 5.0;
            let gap_b = d[i][j - 1] + 5.0;
            
            let mut val = match_score;
            if gap_a < val { val = gap_a; }
            if gap_b < val { val = gap_b; }
            d[i][j] = val;
        }
    }
    
    let mut min_val = f32::INFINITY;
    let mut i_star = 0;
    for i in 0..=m {
        if d[i][n] < min_val {
            min_val = d[i][n];
            i_star = i;
        }
    }
    
    let mut curr_i = i_star;
    let mut curr_j = n;
    while curr_j > 0 {
        if curr_i > 0 && curr_j > 0 {
            let cost = rvd_mismatch_cost(&seq_a[curr_i - 1], &seq_b[curr_j - 1]);
            if (d[curr_i][curr_j] - (d[curr_i - 1][curr_j - 1] + cost)).abs() < 1e-4 {
                curr_i -= 1;
                curr_j -= 1;
                continue;
            }
        }
        if curr_i > 0 {
            if (d[curr_i][curr_j] - (d[curr_i - 1][curr_j] + 5.0)).abs() < 1e-4 {
                curr_i -= 1;
                continue;
            }
        }
        if curr_j > 0 {
            if (d[curr_i][curr_j] - (d[curr_i][curr_j - 1] + 5.0)).abs() < 1e-4 {
                curr_j -= 1;
                continue;
            }
        }
        if curr_i > 0 && curr_j > 0 {
            curr_i -= 1;
            curr_j -= 1;
        } else if curr_i > 0 {
            curr_i -= 1;
        } else {
            curr_j -= 1;
        }
    }
    let i_start = curr_i;
    
    let mut alignment_cost = d[i_star][n];
    if i_start > 0 {
        alignment_cost += 1.0 + 0.1 * i_start as f32;
    }
    if i_star < m {
        alignment_cost += 1.0 + 0.1 * (m - i_star) as f32;
    }
    
    alignment_cost
}

// UPGMA Tree node representation
enum Node {
    Leaf {
        id: usize,
        _tale_name: String,
        length: usize,
    },
    Internal {
        left: Box<Node>,
        right: Box<Node>,
        merge_dist: f32,
        max_tale_len: usize,
        members: Vec<usize>,
    },
}

impl Node {
    fn members(&self) -> &[usize] {
        match self {
            Node::Leaf { id, .. } => std::slice::from_ref(id),
            Node::Internal { members, .. } => members,
        }
    }
    
    fn max_tale_len(&self) -> usize {
        match self {
            Node::Leaf { length, .. } => *length,
            Node::Internal { max_tale_len, .. } => *max_tale_len,
        }
    }
}

// Calculate the average linkage distance between two nodes
fn cluster_distance(c1: &Node, c2: &Node, dist_matrix: &[Vec<f32>]) -> f32 {
    let m1 = c1.members();
    let m2 = c2.members();
    let mut sum = 0.0;
    for &u in m1 {
        for &v in m2 {
            sum += dist_matrix[u][v];
        }
    }
    sum / (m1.len() * m2.len()) as f32
}

// Helper to recursively collect leaf IDs from a node
fn collect_leaves(node: &Node, leaves: &mut Vec<usize>) {
    match node {
        Node::Leaf { id, .. } => {
            leaves.push(*id);
        }
        Node::Internal { left, right, .. } => {
            collect_leaves(left, leaves);
            collect_leaves(right, leaves);
        }
    }
}

// Process the UPGMA dendrogram using the cut threshold and length mismatch splitter
fn process_node(node: &Node, cut: f32, families: &mut Vec<Vec<usize>>) {
    match node {
        Node::Leaf { id, .. } => {
            families.push(vec![*id]);
        }
        Node::Internal { left, right, merge_dist, .. } => {
            // 1. If merge distance is greater than the cut threshold, we must split
            if *merge_dist > cut {
                process_node(left, cut, families);
                process_node(right, cut, families);
                return;
            }
            
            // 2. Apply Length Mismatch Splitting Filter
            let max_l1 = left.max_tale_len() as f32;
            let max_l2 = right.max_tale_len() as f32;
            let min_max_len = max_l1.min(max_l2);
            
            let ratio = if min_max_len > 0.0 {
                *merge_dist / min_max_len
            } else {
                0.0
            };
            
            if ratio >= 0.3 {
                // Split clusters recursively
                process_node(left, cut, families);
                process_node(right, cut, families);
            } else {
                // Cohesive family
                let mut family = Vec::new();
                collect_leaves(node, &mut family);
                families.push(family);
            }
        }
    }
}

// Helper to compute a simple majority-vote consensus RVD sequence for a family
fn get_consensus(rvds_list: &[&Vec<String>]) -> String {
    if rvds_list.is_empty() {
        return String::new();
    }
    let max_len = rvds_list.iter().map(|r| r.len()).max().unwrap_or(0);
    let mut consensus = Vec::new();
    for pos in 0..max_len {
        let mut counts = std::collections::HashMap::new();
        for seq in rvds_list {
            if let Some(rvd) = seq.get(pos) {
                *counts.entry(rvd).or_insert(0) += 1;
            }
        }
        if let Some((most_frequent, _)) = counts.into_iter().max_by_key(|&(_, count)| count) {
            consensus.push(most_frequent.clone());
        }
    }
    consensus.join("-")
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !args.outdir.exists() {
        std::fs::create_dir_all(&args.outdir)
            .with_context(|| format!("Failed to create output directory {:?}", args.outdir))?;
    }

    let reader = fasta::Reader::from_file(&args.input)
        .with_context(|| format!("Failed to read input FASTA file {:?}", args.input))?;

    let mut tales: Vec<(String, Vec<String>)> = Vec::new();

    for record in reader.records() {
        let rec = record?;
        let id = rec.id().to_string();
        let seq = rec.seq();
        if is_dna(seq) {
            let rvds = dna_to_rvds(seq);
            tales.push((id, rvds));
        } else {
            let seq_str = String::from_utf8_lossy(seq);
            let rvds = parse_rvd_sequence(&seq_str);
            tales.push((id, rvds));
        }
    }

    if tales.is_empty() {
        println!("No TALE sequences found in input.");
        return Ok(());
    }

    println!("Loaded {} TALE sequences.", tales.len());

    if tales.len() == 1 {
        // Handle single TALE edge case
        let (id, rvds) = &tales[0];
        let tsv_path = args.outdir.join("family_assignments.tsv");
        let mut tsv_writer = BufWriter::new(File::create(&tsv_path)?);
        writeln!(tsv_writer, "TALE\tFamily\tRVDs")?;
        writeln!(tsv_writer, "{}\tFamily_1\t{}", id, rvds.join("-"))?;
        
        let summary_path = args.outdir.join("family_summary.txt");
        let mut summary_writer = BufWriter::new(File::create(&summary_path)?);
        writeln!(summary_writer, "================================================================================")?;
        writeln!(summary_writer, "Family 1 (1 members)")?;
        writeln!(summary_writer, "================================================================================")?;
        writeln!(summary_writer, "Consensus: {}", rvds.join("-"))?;
        writeln!(summary_writer, "Members:")?;
        writeln!(summary_writer, "  - {} ({} repeats): {}", id, rvds.len(), rvds.join("-"))?;
        
        println!("Clustering complete. Assignments and summary written to {:?}", args.outdir);
        return Ok(());
    }

    println!("Calculating pairwise glocal dynamic programming distance matrix in parallel...");
    let n_tales = tales.len();
    let mut dist_matrix = vec![vec![0.0; n_tales]; n_tales];
    
    let mut pairs = Vec::new();
    for i in 0..n_tales {
        for j in (i + 1)..n_tales {
            pairs.push((i, j));
        }
    }
    
    let pair_distances: Vec<f32> = pairs
        .par_iter()
        .map(|&(i, j)| glocal_distance(&tales[i].1, &tales[j].1))
        .collect();
        
    for (idx, &(i, j)) in pairs.iter().enumerate() {
        let dist = pair_distances[idx];
        dist_matrix[i][j] = dist;
        dist_matrix[j][i] = dist;
    }

    println!("Building UPGMA average-linkage hierarchical dendrogram...");
    let mut active_clusters = Vec::new();
    for (i, (id, rvds)) in tales.iter().enumerate() {
        active_clusters.push(Node::Leaf {
            id: i,
            _tale_name: id.clone(),
            length: rvds.len(),
        });
    }

    while active_clusters.len() > 1 {
        let mut min_dist = f32::INFINITY;
        let mut merge_pair = (0, 0);

        let n_active = active_clusters.len();
        for i in 0..n_active {
            for j in (i + 1)..n_active {
                let dist = cluster_distance(&active_clusters[i], &active_clusters[j], &dist_matrix);
                if dist < min_dist {
                    min_dist = dist;
                    merge_pair = (i, j);
                }
            }
        }

        let (idx1, idx2) = merge_pair;
        let right_node = active_clusters.remove(idx2);
        let left_node = active_clusters.remove(idx1);

        let mut members = Vec::with_capacity(left_node.members().len() + right_node.members().len());
        members.extend_from_slice(left_node.members());
        members.extend_from_slice(right_node.members());

        let max_tale_len = std::cmp::max(left_node.max_tale_len(), right_node.max_tale_len());

        let parent = Node::Internal {
            merge_dist: min_dist,
            max_tale_len,
            members,
            left: Box::new(left_node),
            right: Box::new(right_node),
        };

        active_clusters.push(parent);
    }

    let root = active_clusters.remove(0);

    println!("Traversing tree with cut threshold {} and normalized length mismatch splitting...", args.cut);
    let mut families = Vec::new();
    process_node(&root, args.cut, &mut families);

    // Sort families by member count in descending order
    families.sort_by(|f1, f2| f2.len().cmp(&f1.len()));

    println!("Formed {} homologous TALE families.", families.len());

    // Generate output files
    let tsv_path = args.outdir.join("family_assignments.tsv");
    let mut tsv_writer = BufWriter::new(File::create(&tsv_path)?);
    writeln!(tsv_writer, "TALE\tFamily\tRVDs")?;

    for (fam_idx, family) in families.iter().enumerate() {
        let family_id = format!("Family_{}", fam_idx + 1);
        for &member_idx in family {
            let (id, rvds) = &tales[member_idx];
            writeln!(tsv_writer, "{}\t{}\t{}", id, family_id, rvds.join("-"))?;
        }
    }

    let summary_path = args.outdir.join("family_summary.txt");
    let mut summary_writer = BufWriter::new(File::create(&summary_path)?);

    for (fam_idx, family) in families.iter().enumerate() {
        let family_id = format!("Family_{}", fam_idx + 1);
        let rvds_list: Vec<&Vec<String>> = family.iter().map(|&idx| &tales[idx].1).collect();
        let consensus = get_consensus(&rvds_list);

        writeln!(summary_writer, "================================================================================")?;
        writeln!(summary_writer, "{} ({} members)", family_id, family.len())?;
        writeln!(summary_writer, "================================================================================")?;
        writeln!(summary_writer, "Consensus: {}", consensus)?;
        writeln!(summary_writer, "Members:")?;
        for &member_idx in family {
            let (id, rvds) = &tales[member_idx];
            writeln!(summary_writer, "  - {} ({} repeats): {}", id, rvds.len(), rvds.join("-"))?;
        }
        writeln!(summary_writer)?;
    }

    println!("Clustering complete. Outputs written to {:?}", args.outdir);
    Ok(())
}
