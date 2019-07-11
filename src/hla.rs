// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.
use std::collections::{HashMap,HashSet};
use std::io::{self, Read, Write};
use failure::Error;
use serde::{Serialize};
use bio::io::fasta;
use debruijn::dna_string::DnaString;
use itertools::Itertools;
use std::fmt;
use regex::Regex;
use std::str::FromStr;
use std::string::String;

#[derive(Clone, Debug, Serialize, Deserialize, Eq, PartialEq, Ord, PartialOrd)]
pub struct Allele {
    pub gene: Vec<u8>,
    pub f1: u16,
    pub f2: Option<u16>,
    pub f3: Option<u16>,
    pub f4: Option<u16>,
    pub name : Vec<u8>,
}

impl fmt::Display for Allele {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", String::from_utf8(self.name.clone()).unwrap())
    }
}

pub struct AlleleDb {
    alleles: Vec<Allele>,
}

pub fn all_same<T: Eq>(mut items: impl Iterator<Item=T>) -> Option<T> {
    let first_item = items.next();
    loop {
        let next_item = items.next();
        if next_item.is_none() {
            return first_item;
        }

        if next_item != first_item {
            return None
        }
    }
}

impl AlleleDb {
    pub fn lowest_common_allele(&self, eq_classes: &[usize]) -> Option<Allele> {

        if eq_classes.is_empty() {
            return None
        }

        // Unique allele hit.
        if eq_classes.len() == 1 {
            return Some(self.alleles[eq_classes[0]].clone())
        }

        let gene = all_same(eq_classes.iter().map(|c| &self.alleles[*c].gene));
        let f1 = all_same(eq_classes.iter().map(|c| &self.alleles[*c].f1));

        None
    }
}

pub struct AlleleParser {
    valid_regex: Regex,
    field_regex: Regex,
}

impl AlleleParser {
    pub fn new() -> AlleleParser {
        let valid_regex = Regex::new("^[A-Z0-9]+[*][0-9]+(:[0-9]+)*[A-Z]?$").unwrap();
        let field_regex = Regex::new("[0-9]+(:[0-9]+)*").unwrap();
        AlleleParser { valid_regex, field_regex }
    }

    pub fn parse(&self, s: &str) -> Result<Allele, Error> {
        
        if !self.valid_regex.is_match(s) {
            return Err(format_err!("invalid allele string: {}", s));
        }

        let mut star_split = s.split('*');
        let gene = star_split.next().ok_or_else(|| format_err!("no split: {}", s))?;
        let suffix = star_split.next().ok_or_else(|| format_err!("invalid allele no star separator: {}", s))?;

        let flds = self.field_regex.find(suffix).ok_or_else(|| format_err!("no alleles found {}", s))?;

        let fld_str = flds.as_str();
        let mut flds = fld_str.split(':');
        let f1 = u16::from_str(flds.next().unwrap()).unwrap();
        let f2 = flds.next().map(|f| u16::from_str(f).unwrap());
        let f3 = flds.next().map(|f| u16::from_str(f).unwrap());
        let f4 = flds.next().map(|f| u16::from_str(f).unwrap());
        
        Ok(Allele {
            gene: gene.as_bytes().to_vec(),
            f1, f2, f3, f4,
            name: s.as_bytes().to_vec(),
        })
    }
}


// Parse headers of the form:
// >HLA:HLA01534 A*02:53N 1098 bp
// Get HLA allele sequences from:
// ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta
pub fn read_hla_cds(
    reader: fasta::Reader<impl Read>,
    allele_set: HashSet<String>,
) -> Result<(Vec<DnaString>, Vec<String>, HashMap<String, Allele>), Error> {
    let mut seqs = Vec::new();
    let mut transcript_counter = 0;
    let mut tx_ids = Vec::new();
    let mut tx_to_allele_map = HashMap::new();

    let allele_parser = AlleleParser::new();

    info!("Starting reading the Fasta file\n");
    let mut hlas = Vec::new();
    for result in reader.records() {
        // obtain record or fail with error
        let record = result?;

        // Sequence
        //let dna_string = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
        let dna_string = DnaString::from_acgt_bytes(record.seq());
        let allele_str = record.desc().ok_or_else(|| format_err!("no HLA allele"))?;
        let allele_str = allele_str.split(' ').next().ok_or_else(||format_err!("no HLA allele"))?;
        if !allele_set.contains(&allele_str.to_string()) { continue; }
        let allele = allele_parser.parse(allele_str)?;
        let tx_id = record.id().to_string();
        let data = (allele, tx_id, allele_str.to_string(), dna_string);
        hlas.push(data);
    }

    hlas.sort();

    let mut lengths = HashMap::new();

    // All the alleles with a common 2-digit prefix must have the same length -- they can only differ by an synonymous mutation.
    // Collate these lengths so that we can filter out non-full length sequences.
    for (two_digit, alleles) in &hlas.iter().group_by(|v| (v.0.gene.clone(), v.0.f1, v.0.f2)) {

        let mut ma: Vec<_> = alleles.collect();
        //println!("td: {:?}, alleles: {:?}", two_digit, ma.len());

        // Pick the longest representative
        ma.sort_by_key(|v| v.3.len());
        let longest = ma.pop().unwrap();
        let (_, _, _, dna_string) = longest;

        lengths.insert(two_digit.clone(), dna_string.len());
    }

    for (three_digit, alleles) in &hlas.iter().group_by(|v| (v.0.gene.clone(), v.0.f1, v.0.f2, v.0.f3)) {

        let mut ma: Vec<_> = alleles.collect();
        //println!("td: {:?}, alleles: {:?}", three_digit, ma.len());
        let nalleles = ma.len();

        // Pick the longest representative
        ma.sort_by_key(|v| v.3.len());

        let longest = ma.pop().unwrap();
        let (allele, tx_id, allele_str, dna_string) = longest;

        // Get the length of longest 2-digit entry
        let req_len = lengths[&(three_digit.0.clone(), three_digit.1, three_digit.2)];

        let mylen = dna_string.len();

        //println!("td: {:?}, alleles: {:?}, max_len: {}, req_len: {}", three_digit, nalleles, mylen, req_len);

        if mylen >= req_len {
            //the CDS only has three-digit resolution so having more than that in allele_str is misleading
            //when that's reported as the HLA type
            seqs.push(dna_string.clone());
            tx_ids.push(allele_str.to_string());
            tx_to_allele_map.insert(tx_id.to_string(), allele.clone());
            transcript_counter += 1;
        }
    }

    println!( "Read {} Alleles, deduped into {} full-length 3-digit alleles", hlas.len(), transcript_counter);
    Ok((seqs, tx_ids, tx_to_allele_map))
}

#[cfg(test)]
mod test {
    use super::*;

    const T1: &str = "A*01:01:01:01";

    #[test]
    fn test_parse1() {
        let parser = AlleleParser::new();
        let al = parser.parse(T1).unwrap();
        assert_eq!(String::from_utf8(al.gene).unwrap(), "A");
        assert_eq!(al.f1, 1);
        assert_eq!(al.f2, Some(1));
        assert_eq!(al.f3, Some(1));
        assert_eq!(al.f4, Some(1));
    }


    const T2: &str = "A*01:01:38L";

    #[test]
    fn test_parse2() {
        let parser = AlleleParser::new();
        let al = parser.parse(T2).unwrap();
        assert_eq!(String::from_utf8(al.gene).unwrap(), "A");
        assert_eq!(al.f1, 1);
        assert_eq!(al.f2, Some(1));
        assert_eq!(al.f3, Some(38));
        assert_eq!(al.f4, None);
    }

    const T3: &str = "MICB*012";

    #[test]
    fn test_parse3() {
        let parser = AlleleParser::new();
        let al = parser.parse(T3).unwrap();
        assert_eq!(String::from_utf8(al.gene).unwrap(), "MICB");
        assert_eq!(al.f1, 12);
        assert_eq!(al.f2, None);
        assert_eq!(al.f3, None);
        assert_eq!(al.f4, None);
    }

    const T4: &str = "MICB*012,5";

    #[test]
    fn test_parse4() {
        let parser = AlleleParser::new();
        let al = parser.parse(T4);
        assert!(al.is_err());
    }
}