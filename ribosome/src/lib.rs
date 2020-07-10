#[macro_use]
use amino::{AminoCat, AminoCell, Amino};
use category::Cat;
use polymer::{Polymer, Strand};
use rna::{RNACat, RNACell, RNA};

type RNACodon = (RNACell, RNACell, RNACell);

pub struct Ribosome {
    pub amino_c: AminoCat,
    pub rna_c: RNACat,
}

#[derive(Debug)]
pub enum Segment {
    Protien(Strand<AminoCell>),
    Junk(Strand<RNACell>),
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
        let mut is_junk = true;
        let mut segments = Vec::<Segment>::new();

        let mut current_protien = Strand::<AminoCell>::new();
        let mut current_rna = Strand::<RNACell>::new();
        let mut current_codon = Strand::<RNACell>::new();

        for nucl in r.contents {
            match is_junk {
                true => {
                    current_rna.push(nucl.clone());
                    current_codon.push(nucl.clone());

                    if current_codon.contents.len() == 3 {
                        let x = self.codon_to_amino(&(
                            current_codon.contents[0].clone(),
                            current_codon.contents[1].clone(),
                            current_codon.contents[2].clone(),
                        ));

                        if x == self.amino_c.morphisms.met {
                            is_junk = false;
                            current_rna
                                .contents
                                .truncate(current_rna.contents.len() - 3);

                            if current_rna.contents.len() > 0 {
                                segments.push(Segment::Junk(current_rna));
                            }

                            current_codon = Strand::<RNACell>::new();
                            current_rna = Strand::<RNACell>::new();
                        } else {
                            // TODO: GET FANCIER.

                            let fst = current_codon.contents[1].clone();
                            let snd = current_codon.contents[2].clone();

                            current_codon = Strand::<RNACell>::new();

                            current_codon.push(fst);
                            current_codon.push(snd);
                        }
                    }
                }

                _ => {
                    current_codon.push(nucl.clone());

                    if current_codon.contents.len() == 3 {
                        let x = self.codon_to_amino(&(
                            current_codon.contents[0].clone(),
                            current_codon.contents[1].clone(),
                            current_codon.contents[2].clone(),
                        ));

                        if x == self.amino_c.morphisms.stop {
                            if current_protien.contents.len() > 0 {
                                segments.push(Segment::Protien(current_protien));
                                current_protien = Strand::<AminoCell>::new();
                            }

                            is_junk = true;
                        } else {
                            current_protien.push(x);
                        }

                        current_codon = Strand::<RNACell>::new();
                    }
                }
            }
        }
        Some(segments)
    }

    pub fn codon_eq(&self, fst: &RNACodon, snd: &RNACodon) -> bool {
        self.codon_to_amino(fst) == self.codon_to_amino(snd)
    }

    pub fn amino_codon_eq(&self, a: &AminoCell, c: &RNACodon) -> bool {
        a == &self.codon_to_amino(c)
    }

    pub fn amino_to_condon_vec(&self, a: &AminoCell) -> Option<Vec<RNACodon>> {
        let x = a.read();
        match x.as_ref().unwrap() {
            &Amino::Ala => Some(vec![
                (self.rna_c.g.clone(), self.rna_c.c.clone(), self.rna_c.u.clone()),
                (self.rna_c.g.clone(), self.rna_c.c.clone(), self.rna_c.c.clone()),
                (self.rna_c.g.clone(), self.rna_c.c.clone(), self.rna_c.a.clone()),
                (self.rna_c.g.clone(), self.rna_c.c.clone(), self.rna_c.g.clone()),
            ]),
            &Amino::Arg => Some(vec![
                (self.rna_c.c.clone(), self.rna_c.g.clone(), self.rna_c.u.clone()),
                (self.rna_c.c.clone(), self.rna_c.g.clone(), self.rna_c.c.clone()),
                (self.rna_c.c.clone(), self.rna_c.g.clone(), self.rna_c.a.clone()),
                (self.rna_c.c.clone(), self.rna_c.g.clone(), self.rna_c.g.clone()),
            ]),
            &Amino::Asn => Some(vec![
                (self.rna_c.a.clone(), self.rna_c.a.clone(), self.rna_c.u.clone()),
                (self.rna_c.a.clone(), self.rna_c.a.clone(), self.rna_c.c.clone()),
            ]),
            // TODO: WORK FROM HERE:
            &Amino::Asp => Some(vec![
                
            ]),
            &Amino::Cys => Some(vec![
                
            ]),
            &Amino::Gln => Some(vec![
                
            ]),
            &Amino::Glu => Some(vec![
                
            ]),
            &Amino::Gly => Some(vec![
                
            ]),
            &Amino::His => Some(vec![
                
            ]),
            &Amino::Ile => Some(vec![
                
            ]),
            &Amino::Leu => Some(vec![
                
            ]),
            &Amino::Lys => Some(vec![
                
            ]),
            &Amino::Met => Some(vec![
                
            ]),
            &Amino::Phe => Some(vec![
                
            ]),
            &Amino::Pro => Some(vec![
                
            ]),
            &Amino::Ser => Some(vec![
                
            ]),
            &Amino::Thr => Some(vec![
                
            ]),
            &Amino::Trp => Some(vec![
                
            ]),
            &Amino::Tyr => Some(vec![
                
            ]),
            &Amino::Val => Some(vec![
                
            ]),
            _ => None
        }
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
            (&RNA::U, &RNA::U, &RNA::A) | (&RNA::U, &RNA::U, &RNA::G) => {
                self.amino_c.morphisms.leu.clone()
            }
            (&RNA::U, &RNA::U, &RNA::U) | (&RNA::U, &RNA::U, &RNA::C) => {
                self.amino_c.morphisms.phe.clone()
            }

            // Square 2:
            (&RNA::U, &RNA::C, &RNA::U)
            | (&RNA::U, &RNA::C, &RNA::C)
            | (&RNA::U, &RNA::C, &RNA::A)
            | (&RNA::U, &RNA::C, &RNA::G) => self.amino_c.morphisms.ser.clone(),

            // Square 3:
            (&RNA::U, &RNA::A, &RNA::U) | (&RNA::U, &RNA::A, &RNA::C) => {
                self.amino_c.morphisms.tyr.clone()
            }
            (&RNA::U, &RNA::A, &RNA::A) | (&RNA::U, &RNA::A, &RNA::G) => {
                self.amino_c.morphisms.stop.clone()
            }

            // Square 4:
            (&RNA::U, &RNA::G, &RNA::U) | (&RNA::U, &RNA::G, &RNA::C) => {
                self.amino_c.morphisms.cys.clone()
            }
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
            (&RNA::C, &RNA::A, &RNA::U) | (&RNA::C, &RNA::A, &RNA::C) => {
                self.amino_c.morphisms.his.clone()
            }
            (&RNA::C, &RNA::A, &RNA::A) | (&RNA::C, &RNA::A, &RNA::G) => {
                self.amino_c.morphisms.gln.clone()
            }

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
            (&RNA::A, &RNA::A, &RNA::U) | (&RNA::A, &RNA::A, &RNA::C) => {
                self.amino_c.morphisms.asn.clone()
            }
            (&RNA::A, &RNA::A, &RNA::A) | (&RNA::A, &RNA::A, &RNA::G) => {
                self.amino_c.morphisms.lys.clone()
            }

            // Square 12:
            (&RNA::A, &RNA::G, &RNA::U) | (&RNA::A, &RNA::G, &RNA::C) => {
                self.amino_c.morphisms.ser.clone()
            }
            (&RNA::A, &RNA::G, &RNA::A) | (&RNA::A, &RNA::G, &RNA::G) => {
                self.amino_c.morphisms.arg.clone()
            }

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
            (&RNA::G, &RNA::A, &RNA::U) | (&RNA::G, &RNA::A, &RNA::C) => {
                self.amino_c.morphisms.asp.clone()
            }
            (&RNA::G, &RNA::A, &RNA::A) | (&RNA::G, &RNA::A, &RNA::G) => {
                self.amino_c.morphisms.glu.clone()
            }

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
            ribo.amino_c
                .monomer_from_string(String::from("met"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("thr"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("asp"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("gln"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("pro"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("gln"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("ala"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("glu"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("leu"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("ala"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("phe"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("thr"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("tyr"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("asp"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("ala"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("pro"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("stop"))
                .unwrap(),
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

        match &protien[0] {
            Segment::Protien(x) => {
                assert_eq!(x.contents.len(), 16);
                let result = vec![
                    ribo.amino_c
                        .monomer_from_string(String::from("met"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("thr"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("asp"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("gln"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("pro"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("gln"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("ala"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("glu"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("leu"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("ala"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("phe"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("thr"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("tyr"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("asp"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("ala"))
                        .unwrap(),
                    ribo.amino_c
                        .monomer_from_string(String::from("pro"))
                        .unwrap(),
                ];

                assert_eq!(x.contents, result);
            }
            _ => panic!("Expected Protien, got Junk in test_translate!"),
        }
    }

    #[test]
    fn junk_data_test() {
        let ac = AminoCat::new();
        let rc = RNACat::new();
        let ribo = Ribosome::new(ac, rc);

        let mut seg_string = String::from("augaugacggaucagccgcaagcggaauuggcguuuacguacgaugcgccguaa");
        let junk_string = String::from("gccgccgccgcc");

        seg_string.push_str(&junk_string);
        seg_string.push_str(&seg_string.clone());

        let strand = ribo.rna_c.polymer_from_string(seg_string).unwrap();
        let protien = ribo.translate(strand).unwrap();

        let result = vec![
            ribo.amino_c
                .monomer_from_string(String::from("met"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("thr"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("asp"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("gln"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("pro"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("gln"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("ala"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("glu"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("leu"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("ala"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("phe"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("thr"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("tyr"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("asp"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("ala"))
                .unwrap(),
            ribo.amino_c
                .monomer_from_string(String::from("pro"))
                .unwrap(),
        ];

        match &protien[0] {
            Segment::Junk(_) => panic!("Expected Protien, got Junk in Segment 0 in Junk Test"),
            Segment::Protien(x) => {
                assert_eq!(x.contents.len(), 16);
                assert_eq!(x.contents, result);
            }
        }

        match &protien[1] {
            Segment::Protien(_) => panic!("Expected Junk, got Protien in Segment 1 in Junk Test"),
            Segment::Junk(x) => {
                println!("{:?}", x.contents);
                assert_eq!(x.contents.len(), 12);
            }
        }

        match &protien[2] {
            Segment::Junk(_) => panic!("Expected Protien, got Junk in Segment 2 in Junk Test"),
            Segment::Protien(x) => {
                assert_eq!(x.contents.len(), 16);
                assert_eq!(x.contents, result);
            }
        }
    }
}
