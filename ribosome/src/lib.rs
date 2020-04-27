#[macro_use]
use amino::{AminoCat, AminoCell};
use category::Cat;
use rna::{RNA, RNACat, RNACell};
use polymer::{Polymer, Strand};

type RNACodon = (RNACell, RNACell, RNACell);

pub struct Ribosome {
    pub amino_c: AminoCat,
    pub rna_c: RNACat,
}

pub struct Segment {
    pub junk: bool,
    pub strand: Strand<AminoCell>,
}

impl Ribosome {
    pub fn new(a: AminoCat, r: RNACat) -> Self {
        Ribosome {
            amino_c: a,
            rna_c: r,
        }
    }

    pub fn plain_translate(&self, r: Strand<RNACell>) -> Option<Strand<AminoCell>> {
        match r.contents.len() % 3 == 0 {
            false => None,
            _ => {
                let mut x = r.contents.into_iter().peekable();
                let mut y = Strand::<AminoCell>::new();

                while x.peek().is_some() {
                    let chunk: Vec<RNACell> = x.by_ref().take(3).collect();
                    let z = (chunk[0].clone(), chunk[1].clone(), chunk[2].clone());
                    let a = self.codon_to_amino(&z);
                    y.push(a);
                }

                Some(y)
            }
        }
    }

    pub fn translate(&self, r: Strand<RNACell>) -> Option<Vec<Segment>> {
        match r.contents.len() % 3 == 0 {
            false => None,
            _ => {
                let mut is_junk = true;

                let mut rna_contents = r.contents.into_iter().peekable();
                let mut segments = Vec::<Segment>::new();
                let mut current_seg = Strand::<AminoCell>::new();

                while rna_contents.peek().is_some() {
                    let chunk: Vec<RNACell> = rna_contents.by_ref().take(3).collect();
                    let codon = (chunk[0].clone(), chunk[1].clone(), chunk[2].clone());
                    let amino = self.codon_to_amino(&codon);
                    match (is_junk, amino == self.amino_c.morphisms.met, amino == self.amino_c.morphisms.stop)  {
                        // In a stopped state, found start
                        (true, true, _) =>  {
                            is_junk = false;
                            if current_seg.contents.len() > 0 {
                                segments.push(
                                    Segment {
                                        junk: true,
                                        strand: current_seg,
                                    }
                                );
                                current_seg = Strand::<AminoCell>::new();
                            }
                        }
                        
                        // More Junk:
                        (true, false, _) => {
                            current_seg.push(amino);
                        }

                        // Continuing after start 
                        (false, _, false) => {
                            current_seg.push(amino);
                        }

                        // Stopping after start 
                        (false, _, true) => {
                            is_junk = true;
                            if current_seg.contents.len() > 0 {
                                segments.push(
                                    Segment {
                                        junk: false,
                                        strand: current_seg,
                                    }
                                );
                                current_seg = Strand::<AminoCell>::new();
                            }
                        }
                    }
                }

                Some(segments)
            }
        }
    }

    pub fn codon_eq(&self, fst: &RNACodon, snd: &RNACodon) -> bool {
        self.codon_to_amino(fst) == self.codon_to_amino(snd)
    }

    pub fn amino_codon_eq(&self, a: &AminoCell, c: &RNACodon) -> bool {
        a == &self.codon_to_amino(c)
    }

    pub fn codon_to_amino(&self, r: &RNACodon) -> AminoCell {
        let (fst, snd, thd) = r;

        let x = fst.read();
        let y = snd.read();
        let z = thd.read();

        let plain_codon = (
            x.as_ref().unwrap(),
            y.as_ref().unwrap(),
            z.as_ref().unwrap(),
        );

        match plain_codon {
            // Reference: https://cnx.org/contents/GFy_h8cu@9.87:QEibhJMi@8/The-Genetic-Code
            // Square 1:
            (&RNA::U, &RNA::U, &RNA::A) | (&RNA::U, &RNA::U, &RNA::G) => self.amino_c.morphisms.leu.clone(),
            (&RNA::U, &RNA::U, &RNA::U) | (&RNA::U, &RNA::U, &RNA::C) => self.amino_c.morphisms.phe.clone(),

            // Square 2:
            (&RNA::U, &RNA::C, &RNA::U)
            | (&RNA::U, &RNA::C, &RNA::C)
            | (&RNA::U, &RNA::C, &RNA::A)
            | (&RNA::U, &RNA::C, &RNA::G) => self.amino_c.morphisms.ser.clone(),

            // Square 3:
            (&RNA::U, &RNA::A, &RNA::U) | (&RNA::U, &RNA::A, &RNA::C) => self.amino_c.morphisms.tyr.clone(),
            (&RNA::U, &RNA::A, &RNA::A) | (&RNA::U, &RNA::A, &RNA::G) => self.amino_c.morphisms.stop.clone(),

            // Square 4:
            (&RNA::U, &RNA::G, &RNA::U) | (&RNA::U, &RNA::G, &RNA::C) => self.amino_c.morphisms.cys.clone(),
            (&RNA::U, &RNA::G, &RNA::A) => self.amino_c.morphisms.stop.clone(),
            (&RNA::U, &RNA::G, &RNA::G) => self.amino_c.morphisms.trp.clone(),

            // Square 5:
            (&RNA::C, &RNA::U, &RNA::U)
            | (&RNA::C, &RNA::U, &RNA::C)
            | (&RNA::C, &RNA::U, &RNA::A)
            | (&RNA::C, &RNA::U, &RNA::G) => self.amino_c.morphisms.leu.clone(),

            // Square 6:
            (&RNA::C, &RNA::C, &RNA::U)
            | (&RNA::C, &RNA::C, &RNA::C)
            | (&RNA::C, &RNA::C, &RNA::A)
            | (&RNA::C, &RNA::C, &RNA::G) => self.amino_c.morphisms.pro.clone(),

            // Square 7:
            (&RNA::C, &RNA::A, &RNA::U) | (&RNA::C, &RNA::A, &RNA::C) => self.amino_c.morphisms.his.clone(),
            (&RNA::C, &RNA::A, &RNA::A) | (&RNA::C, &RNA::A, &RNA::G) => self.amino_c.morphisms.gln.clone(),

            // Square 8:
            (&RNA::C, &RNA::G, &RNA::U)
            | (&RNA::C, &RNA::G, &RNA::C)
            | (&RNA::C, &RNA::G, &RNA::A)
            | (&RNA::C, &RNA::G, &RNA::G) => self.amino_c.morphisms.arg.clone(),

            // Square 9:
            (&RNA::A, &RNA::U, &RNA::U)
            | (&RNA::A, &RNA::U, &RNA::C)
            | (&RNA::A, &RNA::U, &RNA::A) => self.amino_c.morphisms.ile.clone(),
            // NOTE: Not start, assumed start has been dealt with
            (&RNA::A, &RNA::U, &RNA::G) => self.amino_c.morphisms.met.clone(),

            // Square 10:
            (&RNA::A, &RNA::C, &RNA::U)
            | (&RNA::A, &RNA::C, &RNA::C)
            | (&RNA::A, &RNA::C, &RNA::A)
            | (&RNA::A, &RNA::C, &RNA::G) => self.amino_c.morphisms.thr.clone(),

            // Square 11:
            (&RNA::A, &RNA::A, &RNA::U) | (&RNA::A, &RNA::A, &RNA::C) => self.amino_c.morphisms.asn.clone(),
            (&RNA::A, &RNA::A, &RNA::A) | (&RNA::A, &RNA::A, &RNA::G) => self.amino_c.morphisms.lys.clone(),

            // Square 12:
            (&RNA::A, &RNA::G, &RNA::U) | (&RNA::A, &RNA::G, &RNA::C) => self.amino_c.morphisms.ser.clone(),
            (&RNA::A, &RNA::G, &RNA::A) | (&RNA::A, &RNA::G, &RNA::G) => self.amino_c.morphisms.arg.clone(),

            // Square 13:
            (&RNA::G, &RNA::U, &RNA::U)
            | (&RNA::G, &RNA::U, &RNA::C)
            | (&RNA::G, &RNA::U, &RNA::A)
            | (&RNA::G, &RNA::U, &RNA::G) => self.amino_c.morphisms.val.clone(),

            // Square 14:
            (&RNA::G, &RNA::C, &RNA::U)
            | (&RNA::G, &RNA::C, &RNA::C)
            | (&RNA::G, &RNA::C, &RNA::A)
            | (&RNA::G, &RNA::C, &RNA::G) => self.amino_c.morphisms.ala.clone(),

            // Square 15:
            (&RNA::G, &RNA::A, &RNA::U) | (&RNA::G, &RNA::A, &RNA::C) => self.amino_c.morphisms.asp.clone(),
            (&RNA::G, &RNA::A, &RNA::A) | (&RNA::G, &RNA::A, &RNA::G) => self.amino_c.morphisms.glu.clone(),

            // Square 16:
            (&RNA::G, &RNA::G, &RNA::U)
            | (&RNA::G, &RNA::G, &RNA::C)
            | (&RNA::G, &RNA::G, &RNA::A)
            | (&RNA::G, &RNA::G, &RNA::G) => self.amino_c.morphisms.gly.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_translate_test() {
        let ac = AminoCat::new();
        let rc = RNACat::new();
        let ribo = Ribosome::new(ac, rc);

        let t_string = String::from("augacggaucagccgcaagcggaauuggcguuuacguacgaugcgccguaa");
        let strand = ribo.rna_c.polymer_from_string(t_string).unwrap();
        let protien = ribo.plain_translate(strand).unwrap();

        assert_eq!(protien.contents.len(), 17);

        let result = vec![
            ribo.amino_c.monomer_from_string(String::from("met")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("thr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("asp")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("gln")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("pro")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("gln")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("glu")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("leu")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("phe")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("thr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("tyr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("asp")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("pro")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("stop")).unwrap(),
        ];

        assert_eq!(protien.contents, result);
    }

    #[test]
    fn translate_test() {
        let ac = AminoCat::new();
        let rc = RNACat::new();
        let ribo = Ribosome::new(ac, rc);

        let t_string = String::from("augaugacggaucagccgcaagcggaauuggcguuuacguacgaugcgccguaa");
        let strand = ribo.rna_c.polymer_from_string(t_string).unwrap();
        let protien = ribo.translate(strand).unwrap();

        assert_eq!(protien[0].strand.contents.len(), 16);

        let result = vec![
            ribo.amino_c.monomer_from_string(String::from("met")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("thr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("asp")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("gln")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("pro")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("gln")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("glu")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("leu")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("phe")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("thr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("tyr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("asp")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("pro")).unwrap(),
        ];

        assert_eq!(protien[0].strand.contents, result);
    }

    #[test]
    fn junk_data_test() {
        let ac = AminoCat::new();
        let rc = RNACat::new();
        let ribo = Ribosome::new(ac, rc);

        let mut seg_string = String::from("augaugacggaucagccgcaagcggaauuggcguuuacguacgaugcgccguaa");
        let junk_string = String::from("guagauuaguag");

        seg_string.push_str(&junk_string);
        seg_string.push_str(&seg_string.clone());

        let strand = ribo.rna_c.polymer_from_string(seg_string).unwrap();
        let protien = ribo.translate(strand).unwrap();

        assert_eq!(protien[0].strand.contents.len(), 16);
        assert_eq!(protien[1].strand.contents.len(), 4);
        assert_eq!(protien[2].strand.contents.len(), 16);

        assert_eq!(protien[0].junk, false);
        assert!(protien[1].junk);
        assert_eq!(protien[2].junk, false);

        let result = vec![
            ribo.amino_c.monomer_from_string(String::from("met")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("thr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("asp")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("gln")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("pro")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("gln")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("glu")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("leu")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("phe")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("thr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("tyr")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("asp")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("ala")).unwrap(),
            ribo.amino_c.monomer_from_string(String::from("pro")).unwrap(),
        ];

        assert_eq!(protien[0].strand.contents, result);
        assert_eq!(protien[2].strand.contents, result);
    }
}