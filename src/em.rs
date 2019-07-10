use std::collections::{HashMap,HashSet};
use bam::EqClass;

use failure::Error;
use crate::{utils, bam, config, pseudoaligner};

#[derive(Default)]
pub struct EqClassCounts {
    pub nitems: usize,
    pub counts: HashMap<EqClass, u32>,
    pub counts_reads: HashMap<EqClass, u32>,
    pub nreads: u32,
}

impl EqClassCounts {
    pub fn new() -> EqClassCounts {
        EqClassCounts {
            nitems: 0,
            counts: HashMap::new(),
            counts_reads: HashMap::new(),
            nreads: 0,
        }
    }
    pub fn pair_reads_explained(&self, a: u32, b: u32) -> u32 {
        let mut exp = 0;
        for (cls, count) in self.counts_reads.iter() {
            if cls.contains(&a) || cls.contains(&b) {
                exp += *count;
            }
        }
        exp
    }
}

#[allow(clippy::needless_range_loop)]
impl EmProblem for EqClassCounts {
    fn init(&self) -> Vec<f64> {
        vec![1.0/(self.nitems as f64); self.nitems]
    }

    fn reg(&self, theta: &mut [f64]) {
        let mut norm = 0.0;

        for i in 0 .. theta.len() {

            // Clamp weight
            let mut v = theta[i];
            if v > 1.0 { v = 1.0 };
            if v < 0.0 { v = 1e-15 };
            theta[i] = v;
            norm += v;
        }

        let inv_norm = 1.0 / norm;
        for i in 0 .. theta.len() {
            theta[i] *= inv_norm;
        }
    }

    fn F(&self, theta1: &[f64], theta2: &mut [f64]) {
        let nitems = self.nitems;

        let mut total_counts = 0.0;

        for i in 0..theta2.len() {
            theta2[i] = 0.0;
        }

        for (class, count) in &self.counts {

            let mut norm = 0.0;
            for tx in class {
                norm += theta1[*tx as usize];
            }

            for tx in class {
                let tx_count = theta1[*tx as usize] / norm * f64::from(*count);
                theta2[*tx as usize] += tx_count;
                total_counts += tx_count;
            }
        }

        let mut max_abs_diff = 0.0;
        let mut max_rel_diff = 0.0;

        for i in 0 .. nitems {
            let old_weights = theta1[i];
            let new_weights = theta2[i] / total_counts;

            let abs_diff = (old_weights - new_weights).abs();
            let rel_diff = abs_diff / old_weights;

            if abs_diff > max_abs_diff {
                max_abs_diff = abs_diff;
            }

            if new_weights > 1e-2 && rel_diff > max_rel_diff {
                max_rel_diff = rel_diff
            }

            theta2[i] = new_weights;
        }
    }

    fn L(&self, theta: &[f64]) -> f64 {
        let mut ll = 0.0;

        for (class, count) in &self.counts {

            // Exclude empty equivalence classes
            if class.is_empty() {
                continue;
            }

            let mut theta_sum = 0.0;
            for tx in class {
                theta_sum += theta[*tx as usize];
            }

            ll += f64::from(*count) * theta_sum.ln();
        }

        ll
    }
}

pub fn em(eqclasses: &EqClassCounts) -> Vec<f64> {
    let nitems = eqclasses.nitems;

    // initialize weights
    let mut weights = vec![1.0/(nitems as f64); nitems];
    
    let mut iters = 0;

    // Abundance required to 'care' about a relative change
    //let rel_care_thresh = 1e-3 / (nitems as f64);

    loop {

        let mut pseudocounts = vec![0.0; nitems];
        let mut total_counts = 0.0;

        for (class, count) in &eqclasses.counts {

            let mut norm = 0.0;
            for tx in class {
                norm += weights[*tx as usize];
            }

            for tx in class {
                let tx_count = weights[*tx as usize] / norm * f64::from(*count);
                pseudocounts[*tx as usize] += tx_count;
                total_counts += tx_count;
            }
        }

        let mut max_abs_diff = 0.0;
        let mut max_rel_diff = 0.0;
        let mut simpsons = 0.0;

        for i in 0 .. nitems {
            let old_weights = weights[i];
            let new_weights = pseudocounts[i] / total_counts;

            let abs_diff = (old_weights - new_weights).abs();
            let rel_diff = abs_diff / old_weights;

            if abs_diff > max_abs_diff {
                max_abs_diff = abs_diff;
            }

            if new_weights > 1e-2 && rel_diff > max_rel_diff {
                max_rel_diff = rel_diff
            }

            weights[i] = new_weights;
            simpsons += new_weights*new_weights;
        }

        let ll = eqclasses.L(&weights);
        debug!("iter: {}, ll: {}, div: {}, rel_diff: {}, abs_diff: {}", iters, ll, 1.0/simpsons, max_rel_diff, max_abs_diff);
        iters += 1;
        if (max_abs_diff < 0.00005 && max_rel_diff < 0.0001) || iters > 5000 {
            break;
        }

        
    }

    weights
}



/// Encapsulate an EM optimization problem so that it can run through an accelerated EM loop (SquareM).
pub trait EmProblem {

    // Make an initial estimate of the parameter vector. May be a naive estimate.
    fn init(&self) -> Vec<f64>;
    
    // Regularize a parameter vector -- fixup an inconsistencies in the parameters.
    // E.g. legalize values or enforce global constrains.
    fn reg(&self, theta: &mut[f64]);

    // Update the parameters -- one EM step
    fn F(&self, theta1: &[f64], theta_2: &mut [f64]);

    // Compute the likelihood of a parameter set
    fn L(&self, theta: &[f64]) -> f64;
}



/// SquareM EM acceleration method. 
/// As described in:
/// Varadhan, Ravi, and Christophe Roland. 
/// "Simple and globally convergent methods for accelerating the convergence of any EM algorithm." 
/// Scandinavian Journal of Statistics 35.2 (2008): 335-353.
/// Takes an implementation of `EmProblem` and applies the accelerated EM algorithm.
pub fn squarem<T: EmProblem>(p: &T) -> Vec<f64> {


    // Array for holding theta
    let mut theta0 = p.init();
    let mut theta1 = theta0.clone();
    let mut theta2 = theta0.clone();
    let mut theta_sq = theta0.clone();

    let mut r = theta0.clone();
    let mut v = theta0.clone();
    
    let n = theta0.len();
    //let prevLL = std::f64::MIN;
    let mut iters = 0; 

    let rel_care_thresh = Some(1e-5);

    loop {

        // Get theta1
        p.F(&theta0, &mut theta1);
        p.F(&theta1, &mut theta2);

        let mut rsq: f64 = 0.0;
        let mut vsq: f64 = 0.0;

        for i in 0..n {
            r[i] = theta1[i] - theta0[i];
            rsq += r[i].powi(2);

            v[i] = theta2[i] - theta1[i] - r[i];
            vsq += v[i].powi(2);
        }

        let mut alpha = -rsq.sqrt() / vsq.sqrt();
        let mut alpha_tries = 1;
        let mut lsq = 0.0;
        let mut l2 = 0.0;

        loop {
            
            let alpha_sq = alpha.powi(2);

            for i in 0..n {
                theta_sq[i] = theta0[i] - 2.0 * alpha * r[i] + alpha_sq * v[i]
            }

            p.reg(&mut theta_sq);

            lsq = p.L(&theta_sq);
            l2 = p.L(&theta2);

            if lsq > l2 || alpha_tries > 5 { 
                break;
            } else {
                alpha = (alpha + -1.0) / 2.0;
            }

            alpha_tries += 1;
        }


        let (max_rel_diff, max_abs_diff) = 
            if lsq > l2 {
                let diff = diffs(&theta0, &theta_sq, rel_care_thresh);
                std::mem::swap(&mut theta0, &mut theta_sq);
                diff
            } else {
                let diff = diffs(&theta0, &theta2, rel_care_thresh);
                std::mem::swap(&mut theta0, &mut theta2);
                diff
            };


        debug!("iter: {}, ll2: {}, llsq: {}, alpha_tries: {}, rel_diff: {}, abs_diff: {}", iters, l2, lsq, alpha_tries, max_rel_diff, max_abs_diff);
        iters += 1;

        if (max_abs_diff < 5e-4 && max_rel_diff < 5e-3) || iters > 2000 {
            break;
        }
    }

    theta0
}


/// Compute the change in the parameter vectors, returning the largest relative and absolute change, respectively.
/// Only parameters with a value greater than rel_thresh (if set), are counted in the relative change check.
fn diffs(t1: &[f64], t2: &[f64], rel_thresh: Option<f64>) -> (f64, f64) {

    let mut max_abs_diff = 0.0;
    let mut max_rel_diff = 0.0;

    for i in 0 .. t1.len() {
        let old_weights = t1[i];
        let new_weights = t2[i];

        let abs_diff = (old_weights - new_weights).abs();
        let rel_diff = abs_diff / old_weights;

        if abs_diff > max_abs_diff {
            max_abs_diff = abs_diff;
        }

        if rel_thresh.map_or(true, |thresh| new_weights > thresh) && rel_diff > max_rel_diff {
            max_rel_diff = rel_diff
        }
    }

    (max_rel_diff, max_abs_diff)
}

pub fn hla_em(hla_index: String, hla_counts: String, cds_db: String, gen_db: String) -> Result<(&'static str,&'static str), Error> {
    let index: pseudoaligner::Pseudoaligner<config::KmerType> = utils::read_obj(&hla_index)?;
    let mut eq_counts: bam::EqClassDb = utils::read_obj(&hla_counts)?;
    let mut genes_to_eq_classes : HashMap<&str,Vec<u32>> = HashMap::new();
    let mut single_gene_eq_classes : Vec<usize> = Vec::new();
    
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;
    let mut weights_file = BufWriter::new(File::create("weights.tsv")?);
    let mut pairs_file = BufWriter::new(File::create("pairs.tsv")?);
    let mut alleles_called: HashSet<String> = HashSet::new();
    
    // identify eq classes that have all transcripts from the same gene
    for (i, tx_ids) in index.eq_classes.clone().iter().enumerate() {
        let mut gene_names : HashSet<&str> = HashSet::new();
        for t in tx_ids {
            let s = &index.tx_names[*t as usize];
            let mut star_split = s.split('*');
            let gene = star_split.next().ok_or_else(|| format_err!("no split: {}", s))?;
            gene_names.insert(gene);
        }
        if gene_names.len() == 1 { 
            single_gene_eq_classes.push(i); 
            for g in gene_names.iter() {
                let c = genes_to_eq_classes.entry(g).or_insert_with(Vec::new);
                c.push(i as u32);
            }
        }
    }
    
    // count the number of transcripts for each gene
    let mut gene_tx_counts : HashMap<&str,Vec<&str>> = HashMap::new();
    for s in &index.tx_names {
        let mut star_split = s.split('*');
        let gene = star_split.next().ok_or_else(|| format_err!("no split: {}", s))?;
        let c = gene_tx_counts.entry(gene).or_insert_with(Vec::new);
        c.push(s.as_str());
    }
    debug!("{} out of {} equivalence classes correspond to a single gene", single_gene_eq_classes.len(), index.eq_classes.len());
    
    for (gene, eqs) in genes_to_eq_classes.iter() {
        let txs = &gene_tx_counts[gene];
        info!("Computing EM for gene {} with {} gene-specific equivalence classes and {} transcripts",gene,eqs.len(),txs.len());
        /*let eq_class_counts = eq_counts.eq_class_counts(index.tx_names.len());
        let weights = squarem(&eq_class_counts);
        let mut weight_names: Vec<(f64, &String, usize)> = 
                weights.
                into_iter().
                enumerate().
                map(|(i,w)| (w, &index.tx_names[i], i)).
                collect();
            weight_names.sort_by(|(wa,_, _), (wb, _, _)| (-wa).partial_cmp(&-wb).unwrap());
        for (w, name, _) in &weight_names {
            if *w > 1e-2 {
                println!("{}\t{}", name, w);
            } else {break;}
        }*/
        
        if let Some((eqclass_counts, tx_order, reads_explained)) = eq_counts.eq_class_counts_from_vec(eqs, &index.eq_classes) {
            
            /*for (cls, n) in eqclass_counts.counts_reads.iter() {
                for i in cls { print!("{} ", &index.tx_names[tx_order[*i as usize] as usize]); }
                println!("{}",n);
            }
            println!("");
            for (cls, n) in eqclass_counts.counts.iter() {
                for i in cls { print!("{} ", &index.tx_names[tx_order[*i as usize] as usize]); }
                println!("{}",n);
            }*/
            
            let weights = squarem(&eqclass_counts);
            debug!("Calculated weights for {} transcripts",weights.len());

            let mut weight_names: Vec<(f64, &String, usize, u32)> = 
                weights.
                into_iter().
                enumerate().
                map(|(i,w)| (w, &index.tx_names[tx_order[i] as usize], reads_explained[i], i as u32)).
                collect();
            weight_names.sort_by(|(wa,_, _,_), (wb, _, _,_)| (-wa).partial_cmp(&-wb).unwrap());
        
            let mut weights_written = 0;
            for (w, name, exp, _) in &weight_names {
                if *w > 1e-2 {
                    weights_written += 1;
                    writeln!(weights_file, "{}\t{}\t{}", name, w, f64::from(*exp as u32) / f64::from(eqclass_counts.nreads))?;
                } else {
                    break;
                }
            }
            
            let mut pairs : Vec<(&String,&String,f64)> = Vec::new();
            
            for i in 0..weights_written-1 {
                for j in i+1..weights_written {
                    let exp = eqclass_counts.pair_reads_explained(weight_names[i].3, weight_names[j].3);
                    pairs.push((weight_names[i].1, weight_names[j].1,f64::from(exp) / f64::from(eqclass_counts.nreads)));
                }
            }
            pairs.sort_by(|(_,_,wa), (_,_,wb)| (-wa).partial_cmp(&-wb).unwrap());
            for i in 0..std::cmp::min(5,pairs.len()) {
                writeln!(pairs_file, "{}\t{}\t{}", pairs[i].0, pairs[i].1, pairs[i].2 )?;
            }
            // Write the sequences
            if !pairs.is_empty() {
                alleles_called.insert(pairs[0].0.clone());
                alleles_called.insert(pairs[0].1.clone());
            } else if !weight_names.is_empty() {
                alleles_called.insert(weight_names[0].1.clone());
            }
        }
    }
    info!("Writing CDS of {} alleles to file", alleles_called.len());
    use bio::io::fasta;
    let mut gen_file = fasta::Writer::to_file("gen_pseudoaln.fasta")?;
    let mut cds_file = fasta::Writer::to_file("cds_pseudoaln.fasta")?;
    let cds_db_file = fasta::Reader::from_file(&cds_db)?;
    let gen_db_file = fasta::Reader::from_file(&gen_db)?;
    
    for result in cds_db_file.records() {
        let record = result?;
        let allele_str = record.desc().ok_or_else(|| format_err!("no HLA allele"))?;
        let allele_str = allele_str.split(' ').next().ok_or_else(||format_err!("no HLA allele"))?;
        if alleles_called.contains(allele_str) {
            println!("Writing CDS of {}",allele_str);
            cds_file.write_record(&record)?;
        }
    }
    cds_file.flush()?;
    info!("Writing genomic sequence of {} alleles to file", alleles_called.len());
    for result in gen_db_file.records() {
        let record = result?;
        let allele_str = record.desc().ok_or_else(|| format_err!("no HLA allele"))?;
        let allele_str = allele_str.split(' ').next().ok_or_else(||format_err!("no HLA allele"))?;
        if alleles_called.contains(allele_str) {
            println!("Writing genomic sequence of {}",allele_str);
            gen_file.write_record(&record)?;
        }
    }
    gen_file.flush()?;
    Ok(("gen_pseudoaln.fasta","cds_pseudoaln.fasta"))
}

#[cfg(test)]
mod test_em {
    use super::*;
    use bam::EqClass;

    fn test_ds() -> EqClassCounts {
        let mut counts = HashMap::new();

        let eqA = EqClass::from(vec![0]);
        let eqAB = EqClass::from(vec![0,1]);

        let eqC = EqClass::from(vec![2]);
        let eqD = EqClass::from(vec![3]);

        counts.insert(eqA, 1);
        counts.insert(eqAB, 19);
        counts.insert(eqC, 10);
        counts.insert(eqD, 10);
        
        EqClassCounts { counts, nitems: 4, counts_reads: HashMap::new(), nreads: 0 }
    }


    fn test2_ds() -> EqClassCounts {
        let mut counts = HashMap::new();

        let eqA = EqClass::from(vec![0]);
        let eqAB = EqClass::from(vec![0,1]);

        let eqC = EqClass::from(vec![2]);
        let eqD = EqClass::from(vec![3]);

        let eqE = EqClass::from(vec![4,5]);

        counts.insert(eqA, 1);
        counts.insert(eqAB, 19);
        counts.insert(eqC, 10);
        counts.insert(eqD, 10);
        counts.insert(eqE, 20);
        
        EqClassCounts { counts, nitems: 6, counts_reads: HashMap::new(), nreads: 0 }
    }


    #[test]
    fn simple_inf() {

        let eqc = test_ds();
        let res = em(&eqc);

        println!("{:?}", res);
    }

    #[test]
    fn med_inf() {

        let eqc = test2_ds();
        let res = em(&eqc);

        println!("{:?}", res);
    }


    #[test]
    fn accel_inf() {

        let eqc = test_ds();
        let res = squarem(&eqc);

        println!("{:?}", res);
    }


}
