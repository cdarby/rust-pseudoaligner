
use rust_htslib::bam::Read;
use rust_htslib::bam::{Record, IndexedReader};
use failure::Error;
use debruijn::dna_string::DnaString;
use pseudoaligner::Pseudoaligner;
use config::KmerType;
use std::path::Path;
use std::str::FromStr;
use std::collections::{HashMap};
//use shardio::{ShardWriter, SortKey};
//use std::borrow::Cow;
use std::path::PathBuf;
use smallvec::SmallVec;

use crate::locus::Locus;
use crate::pseudoaligner::intersect;
use crate::utils;

pub struct BamSeqReader {
    reader: IndexedReader,
    tmp_record: Record,
}

impl BamSeqReader  {
    pub fn new(reader: IndexedReader) -> BamSeqReader {

        BamSeqReader {
            reader,
            tmp_record: Record::new(),
        }
    }

    pub fn fetch(&mut self, locus: &Locus) {
        let tid = self.reader.header().tid(locus.chrom.as_bytes()).unwrap();
        self.reader.fetch(tid, locus.start, locus.end);
    }
}


pub type Barcode = SmallVec<[u8; 24]>;
pub type Umi = SmallVec<[u8; 16]>;
pub type EqClass = SmallVec<[u32; 4]>;

pub const PROC_BC_SEQ_TAG: &[u8]  = b"CB";
pub const PROC_UMI_SEQ_TAG: &[u8] = b"UB";


#[derive(Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BamCrRead {
    key: BamCrReadKey,
    sequence: DnaString,
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BamCrReadKey {
    barcode: Barcode,
    umi: Umi,
}

/*struct Key;
impl SortKey<BamCrRead> for Key {
    type Key = BamCrReadKey;
    fn sort_key(item: &BamCrRead) -> Cow<BamCrReadKey> {
        Cow::Borrowed(&item.key)
    }
}*/


#[derive(Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BusCount {
    barcode_id: u32,
    umi_id: u32,
    eq_class_id: u32,
}

#[derive(Serialize, Deserialize)]
pub struct EqClassDb {
    barcodes: HashMap<Barcode, u32>,
    umis: HashMap<Umi, u32>,
    pub eq_classes: HashMap<EqClass, u32>,
    counts: Vec<BusCount>,
}

impl<'a> EqClassDb {
    pub fn new(_nitems: usize) -> EqClassDb {
        EqClassDb {
            barcodes: HashMap::new(),
            umis: HashMap::new(),
            eq_classes: HashMap::new(),
            counts: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }
    
    pub fn is_empty(&self) -> bool {
        self.counts.len() == 0
    }

    pub fn count(&mut self, barcode: &Barcode, umi: &Umi, eq_class: &EqClass) {
        let barcode_id = match self.barcodes.get(barcode).cloned() {
            Some(id) => id,
            None => {
                let id = self.barcodes.len() as u32;
                self.barcodes.insert(barcode.clone(), id);
                id 
            }
        };

        let umi_id = match self.umis.get(umi).cloned() {
            Some(id) => id,
            None => {
                let id = self.umis.len() as u32;
                self.umis.insert(umi.clone(), id);
                id
            }
        };

        let eq_class_id = match self.eq_classes.get(eq_class).cloned() {
            Some(id) => id,
            None => {
                let id = self.eq_classes.len() as u32;
                self.eq_classes.insert(eq_class.clone(), id);
                id 
            },
        };

        let count = BusCount {
            barcode_id,
            umi_id,
            eq_class_id,
        };

        self.counts.push(count);
    }

    /*
    pub fn get(&self, i: usize) -> (&Barcode, &Umi, &EqClass) {

        let ids = &self.counts[i];

        let bc = self.barcodes.get_by_right(&ids.barcode_id).unwrap();
        let umi = self.umis.get_by_right(&ids.umi_id).unwrap();
        let eq_class = self.eq_classes.get_by_right(&ids.eq_class_id).unwrap();

        (bc, umi, eq_class)
    }
    */

    pub fn sort(&mut self) {
        self.counts.sort();
    }

    pub fn eq_class_counts(&mut self, nitems: usize) -> crate::em::EqClassCounts {
        let mut rev_map = HashMap::<u32, &EqClass>::new();
        let mut counts: HashMap<EqClass, u32> = HashMap::new();
        let mut counts_reads: HashMap<EqClass, u32> = HashMap::new();

        self.sort();

        for (cls, id) in &self.eq_classes {
            rev_map.insert(*id, cls);
        }

        let mut uniq = 0;
        let mut total_reads = 0;
        let mut buf = Vec::new();
        
        use itertools::Itertools;
        for ((_bc, _umi), mut hits) in &self.counts.iter().group_by(|c| (c.barcode_id, c.umi_id)) {
            
            let c = hits.next().unwrap();
            buf.clear();
            buf.extend(rev_map[&c.eq_class_id]);
            total_reads += 1;

            for c in hits {
                intersect(&mut buf, rev_map[&c.eq_class_id]);
                total_reads += 1;
            }

            let eqclass = EqClass::from_slice(&buf);
            let count = counts.entry(eqclass).or_default();
            *count += 1;

            uniq += 1;
        }

        println!("mean umis/read: {}", f64::from(uniq) / f64::from(total_reads));

        let empty = EqClass::new();
        println!("empty eq_class counts: {:?} of {}", counts.get(&empty), self.counts.len());

        crate::em::EqClassCounts {
            nitems,
            counts_umi : counts,
            counts_reads,
            nreads: total_reads,
        }
    }
    
    pub fn eq_class_counts_updated(&mut self, nitems: usize) -> (crate::em::EqClassCounts, Vec<usize>) {
        let mut rev_map = HashMap::<u32, &EqClass>::new();
        let mut counts_reads: HashMap<EqClass, u32> = HashMap::new();
        let mut counts_umi: HashMap<EqClass, u32> = HashMap::new();
        let mut reads_explained = vec![0; nitems];
        self.sort();

        for (cls, id) in &self.eq_classes {
            rev_map.insert(*id, cls);
        }

        let mut uniq = 0;
        let mut total_reads = 0;
        let mut buf = Vec::new();
        
        use itertools::Itertools;
        for ((_bc, _umi), mut hits) in &self.counts.iter().group_by(|c| (c.barcode_id, c.umi_id)) {
            buf.clear();
            let mut flag : bool = false;
            for c in hits {
                if rev_map.contains_key(&c.eq_class_id) {
                    if flag { intersect(&mut buf, rev_map[&c.eq_class_id].as_ref()); }
                    else { buf.extend(rev_map[&c.eq_class_id].as_ref()); flag = true; }
                    total_reads += 1;
                    for t in rev_map[&c.eq_class_id].as_ref() {
                        reads_explained[*t as usize] += 1;
                    }
                    let eqclass = EqClass::from_slice(rev_map[&c.eq_class_id].as_ref());
                    let count = counts_reads.entry(eqclass).or_default();
                    *count += 1;
                }
            }
            if flag {
                buf.sort();
                buf.dedup();
                let eqclass = EqClass::from_slice(&buf);
                let count = counts_umi.entry(eqclass).or_default();
                *count += 1;
                uniq += 1;
            }
        }

        debug!("mean reads per UMI: {}", f64::from(total_reads) / f64::from(uniq));

        //let empty = EqClass::new();
        //println!("empty eq_class counts: {:?} of {}", counts.get(&empty), self.counts.len());

        (crate::em::EqClassCounts {
            nitems,
            counts_reads,
            counts_umi,
            nreads: total_reads,
        }, reads_explained)
    }
    
    // Takes only reads aligning to the subset of the equivalence classes as indicated in the vector included_eq_classes
    /*pub fn eq_class_counts_from_vec(&mut self, included_eqs: &Vec<u32>, eq_classes: &Vec<Vec<u32>>) -> Option<(crate::em::EqClassCounts, Vec<u32>, Vec<usize>)> {
        let mut rev_map = HashMap::<u32, EqClass>::new();
        let mut counts_reads: HashMap<EqClass, u32> = HashMap::new();
        let mut counts_umi: HashMap<EqClass, u32> = HashMap::new();
        let mut tx_order : Vec<u32> = Vec::new();
        
        self.sort();

        for (i, tx_ids) in eq_classes.iter().enumerate() {
            if included_eqs.contains(&(i as u32)) { 
                tx_order.extend(tx_ids);
            }
        }
        tx_order.sort();
        tx_order.dedup();
        //re-number the transcripts that occur
        for (i, tx_ids) in eq_classes.iter().enumerate() {
            if included_eqs.contains(&(i as u32)) { 
                let mut new_cls : Vec<u32> = Vec::new();
                for c in tx_ids {
                    if let Ok(newnumber) = tx_order.binary_search(c) { new_cls.push(newnumber as u32); }
                }
                rev_map.insert(i as u32, EqClass::from_slice(&new_cls)); 
            }
        }
        
        let mut uniq = 0;
        let mut total_reads = 0;
        let mut buf = Vec::new();
        let mut reads_explained = vec![0; tx_order.len()];
        
        use itertools::Itertools;
        for ((_bc, _umi), mut hits) in &self.counts.iter().group_by(|c| (c.barcode_id, c.umi_id)) {
            buf.clear();
            let mut flag : bool = false;
            for c in hits {
                if rev_map.contains_key(&c.eq_class_id) {
                    //if flag { intersect(&mut buf, rev_map[&c.eq_class_id].as_ref()); }
                    //else { buf.extend(rev_map[&c.eq_class_id].as_ref()); flag = true; }
                    flag = true;
                    buf.extend(rev_map[&c.eq_class_id].as_ref());
                    total_reads += 1;
                    for t in rev_map[&c.eq_class_id].as_ref() {
                        reads_explained[*t as usize] += 1;
                    }
                    let eqclass = EqClass::from_slice(rev_map[&c.eq_class_id].as_ref());
                    let count = counts_reads.entry(eqclass).or_default();
                    *count += 1;
                }
            }
            if flag {
                buf.sort();
                buf.dedup();
                let eqclass = EqClass::from_slice(&buf);
                let count = counts_umi.entry(eqclass).or_default();
                *count += 1;
                uniq += 1;
            }
        }
        if total_reads < 100 { return None }
        info!("{} reads mapped",total_reads);
        info!("mean reads per UMI: {}", f64::from(total_reads) / f64::from(uniq)) ;
        //let empty = EqClass::new();
        //println!("empty eq_class counts: {:?} of {}", counts.get(&empty), self.counts.len());
        //println!("{:?}",counts_reads);
        //println!("{:?}",counts_umi);
        Some((crate::em::EqClassCounts {
            nitems : tx_order.len(), //number of transcripts over all the eq classes 
            counts : counts_umi, // counts_reads.clone(),
            counts_reads,
            nreads : total_reads,
        }, tx_order, reads_explained))
    }*/
}


impl Iterator for BamSeqReader {
    type Item = Result<BamCrRead, Error>;

    fn next(&mut self) -> Option<Self::Item> {

        loop {
            let r = self.reader.read(&mut self.tmp_record);

            if let Err(e) = r {
                if e.is_eof() { return None } else { return Some(Err(e.into())) }
            };
            if self.tmp_record.is_secondary() || self.tmp_record.is_supplementary() {
                continue; //TODO want to skip these?
            }

            // Use reads that have no GN tag, or a GN tag that matches "HLA-*"
            /*let gene_filter =
                match self.tmp_record.aux(GENE_TAG) {
                    Some(gn_aux) => { 
                        
                        let gn_bytes = gn_aux.string();

                        let gn_iter = gn_bytes.split(|x| *x == b';');
                        let mut result = false;
                        for gn in gn_iter {
                            let gn_str = std::str::from_utf8(gn).unwrap();

                            if self.gene_regex.is_match(gn_str) {
                                result = true;
                                break;
                            }
                        }
                        result
                    },
                    None => true,
                };

            if !gene_filter { continue };*/



            // Get original read sequence from record.
            //let mut sequence = DnaString::from_acgt_bytes_hashn(&self.tmp_record.seq().as_bytes(), self.tmp_record.qname());
            let mut sequence = DnaString::from_acgt_bytes(&self.tmp_record.seq().as_bytes());
            if !self.tmp_record.is_reverse() {
                sequence = sequence.reverse();
            }
            //println!("{:?}",sequence);
            let barcode = match self.tmp_record.aux(PROC_BC_SEQ_TAG).map(|x| Barcode::from_slice(x.string())) {
                Some(bc) => bc,
                None => continue,
            };

            let umi = match self.tmp_record.aux(PROC_UMI_SEQ_TAG).map(|x| Umi::from_slice(x.string())) {
                Some(bc) => bc,
                None => continue,
            };

            return Some(Ok(BamCrRead { sequence, key: BamCrReadKey { barcode, umi }}));
        }
    }
}

pub fn map_bam(bam: impl AsRef<Path>, align: Pseudoaligner<KmerType>, locus_string: &Option<String>, outs: &Path) -> Result<(), Error> {
    let mut rdr = IndexedReader::from_path(bam)?;
    rdr.set_threads(4).unwrap();

    let mut itr = BamSeqReader::new(rdr);
    
    if let Some(l) = locus_string {
        let locus = Locus::from_str(l)?;
        debug!("Locus: {:?}", locus);
        itr.fetch(&locus);
    }
    

    let mut hits_file = outs.to_path_buf();
    hits_file.set_extension("counts.bin");

    //let mut reads_file = outs.to_path_buf();
    //reads_file.set_extension("shards.bin");

    //let sw = ShardWriter::<BamCrRead, Key>::new(&reads_file, 16, 64, 1<<20);

    let mut rid = 0;
    let mut some_aln = 0;
    let mut long_aln = 0;

    let mut eq_counts = EqClassDb::new(align.tx_names.len());
    let mut nodes = Vec::new();
    let mut eq_class = Vec::new();

    for _rec in itr {
        rid += 1;
        let rec = _rec?;

        let aln = align.map_read_to_nodes(&rec.sequence, &mut nodes);
        
        if let Some(cov) = aln {
            some_aln += 1;
            
            align.nodes_to_eq_class(&nodes, &mut eq_class);
            if cov > 40 {
                long_aln += 1;

                eq_counts.count(&rec.key.barcode, &rec.key.umi, &EqClass::from_slice(&eq_class));

                /*
                for n in &nodes {
                    let bc = std::str::from_utf8(&rec.barcode).unwrap();
                    let umi = std::str::from_utf8(&rec.umi).unwrap();
                    writeln!(hits, "{}\t{}\t{}\t{}\t{}", bc, umi, rid, n, cov)?;
                }
                */
            }
        }

        if rid % 100_000 == 0 {
            info!("analyzed {} reads. Mapped {}, long {}", rid, some_aln, long_aln);
        }
    }

    info!("analyzed {} reads. Mapped {}, long {}", rid, some_aln, long_aln);

    /*
    for node in align.dbg.iter_nodes() {
        let eq_class_id = node.data();
        use debruijn::Mer;
        let seq_len = node.sequence().len();
        let eq_class = &align.eq_classes[*eq_class_id as usize];

        for tx_id in eq_class.iter() {
            writeln!(idx, "{}\t{}\t{}\t{}\t{}", node.node_id, eq_class.len(), seq_len, tx_id, align.tx_names[*tx_id as usize]);
        }
    }
    */

    crate::utils::write_obj(&eq_counts, &hits_file)?;
    Ok(())
}


/*fn filter_rec(rec: &Record, locus: &Locus) -> bool {

    if rec.mapq() < 30 {
        debug!("skipping read {} due to low mapping quality", 
                String::from_utf8(rec.qname().to_vec()).unwrap());

        false
    }
    else if useful_alignment(locus, &rec).unwrap() == false {
        debug!("skipping read {} due to not being useful", 
                String::from_utf8(rec.qname().to_vec()).unwrap());
        false
    } else {
        true
    }
}


pub fn useful_alignment(locus: &Locus, rec: &Record) -> Result<bool, Error> {
    // filter alignments to ensure that they truly overlap the region of interest
    // for now, overlap will be defined as having an aligned base anywhere in the locus
        let cigar = rec.cigar();
        for i in locus.start..=locus.end {
            // Don't include soft-clips but do include deletions
            let t = cigar.read_pos(i, false, true)?; 
            if t.is_some() {
                return Ok(true)
            }
        }
        Ok(false)
}*/


pub fn hla_map_bam(hla_index: String, outdir: String, bam: String, locus: Option<String>) -> Result<String, Error> {
    info!("Reading index from disk");
    let (index, _tx_allele_map) : (Pseudoaligner<KmerType>, HashMap<String, crate::hla::Allele>) = utils::read_obj(&hla_index)?;
    // tx_allele_map is: "HLA:HLA08579": Allele { gene: [66], f1: 40, f2: Some(223), f3: None, f4: None }, ...
    // index.tx_names is: "B*51:154", "B*51:155", ...
    // index.eq_classes.len() is: 22910
    info!("Finished reading index!");
    info!("Mapping reads from BAM");
    let path = PathBuf::from(outdir.clone());
    map_bam(&bam, index, &locus, &path)?;
    Ok(outdir)
}
