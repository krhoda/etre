#[macro_use]
extern crate enum_display_derive;

use category::{Cat, ICat, NCat};
use monomer::{IMono, Mono, NucleicAcid};
use once_mono::{IMonomer};
use polymer::{Helix, Polymer, Strand};

use std::fmt::Display;

// TODO: Support the other DNA alphabets.
#[derive(Debug, Display, PartialEq, Clone)]
pub enum DNA {
    A,
    T,
    G,
    C,
}

impl NucleicAcid for DNA {
    fn is_g_or_c(&self) -> bool {
        match self {
            DNA::G => true,
            DNA::C => true,
            _ => false,
        }
    }
}

impl Mono for DNA {
    fn from_string(s: String) -> Option<DNA> {
        match s.to_uppercase().as_ref() {
            "A" => Some(DNA::A),
            "T" => Some(DNA::T),
            "C" => Some(DNA::C),
            "G" => Some(DNA::G),
            _ => None,
        }
    }
}

impl IMono for DNA {
    fn inverse(c: &Self) -> Self {
        match c {
            DNA::A => DNA::T,
            DNA::T => DNA::A,
            DNA::C => DNA::G,
            DNA::G => DNA::C,
        }
    }
}

type DNACell = IMonomer<DNA>;

pub struct DNACat {
    a: DNACell,
    t: DNACell,
    c: DNACell,
    g: DNACell,
}

impl Cat<DNACell, Helix<DNACell>> for DNACat {
    fn new() -> Self {
        DNACat {
            a: IMonomer::from_string(String::from("a")).unwrap(),
            t: IMonomer::from_string(String::from("t")).unwrap(),
            c: IMonomer::from_string(String::from("c")).unwrap(),
            g: IMonomer::from_string(String::from("g")).unwrap(),
        }
    }

    fn monomer_from_string(&self, s: String) -> Option<DNACell> {
        match DNA::from_string(s) {
            None => None,
            Some(x) => match x {
                DNA::A => Some(self.a.clone()),
                DNA::T => Some(self.t.clone()),
                DNA::C => Some(self.c.clone()),
                DNA::G => Some(self.g.clone()),
            },
        }
    }

    fn polymer_from_string(&self, s: String) -> Option<Helix<DNACell>> {
        let mut x = Helix::<DNACell>::new();
        let mut ok = true;

        for c in s.chars() {
            match self.monomer_from_string(c.to_string()) {
                Some(y) => x.push(y),
                None => {
                    ok = false;
                    break;
                }
            }
        }

        match ok {
            true => Some(x),
            _ => None,
        }
    }
}

impl ICat<DNACell, Helix<DNACell>> for DNACat {
    fn inverse_m(&self, t: &DNACell) -> DNACell {
        match t.read().as_ref().unwrap() {
            DNA::A => self.t.clone(),
            DNA::T => self.a.clone(),
            DNA::C => self.g.clone(),
            DNA::G => self.c.clone(),
        }
    }

    fn inverse_p(&self, u: &Helix<DNACell>) -> Helix<DNACell> {
        let mut y = Helix::<DNACell>::new();
        for z in &u.strand.contents {
            y.push(self.inverse_m(&z));
        }
        y.strand.contents =  y.strand.contents.into_iter().rev().collect();
        y
    }
}

impl NCat<DNACell, Helix<DNACell>> for DNACat {
    fn gc_content(u: Helix<DNACell>) -> (u64, u64) {
        let mut d = 0;
        let mut n = 0;
        for x in u.strand.contents {
            d = d + 1;
            if DNACell::is_g_or_c(&x) {
                n = n + 1;
            }
        }

        (n, d)
    }
}

impl DNACat {
    pub fn pairs(&self, u: Helix<DNACell>) -> Vec<(DNACell, DNACell)> {
        let mut next = Vec::<(DNACell, DNACell)>::new();
        for x in u.strand.contents {
            next.push((x.clone(), self.inverse_m(&x)));
        }
        next
    }

    pub fn strands(&self, u: Helix<DNACell>) -> (Strand<DNACell>, Strand<DNACell>) {
        let x = self.inverse_p(&u);
        (u.strand, x.strand)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_from_char() {
        let c = DNACat::new();

        assert_eq!(
            c.monomer_from_string(String::from("A"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::A
        );
        assert_eq!(
            c.monomer_from_string(String::from("a"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::A
        );

        assert_eq!(
            c.monomer_from_string(String::from("T"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::T
        );

        assert_eq!(
            c.monomer_from_string(String::from("t"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::T
        );

        assert_eq!(
            c.monomer_from_string(String::from("C"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::C
        );
        assert_eq!(
            c.monomer_from_string(String::from("c"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::C
        );

        assert_eq!(
            c.monomer_from_string(String::from("G"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::G
        );

        assert_eq!(
            c.monomer_from_string(String::from("g"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::G
        );

        assert_eq!(c.monomer_from_string(String::from("d")), None);
    }

    // TODO Add Inverse:
    // Now for polynomial tests.

    #[test]
    fn new_helix_cat() {
        let c = DNACat::new();

        let helix_str = String::from("atcg");
        let h = c.polymer_from_string(helix_str).unwrap();
        let mut h2 = Helix::<DNACell>::new();

        h2.push(c.a.clone());
        h2.push(c.t.clone());
        h2.push(c.c.clone());
        h2.push(c.g.clone());

        assert_eq!(h, h2);
    }

    #[test]
    fn helix_from_string() {
        let c = DNACat::new();

        let helix_str = String::from("gattaca");
        let bad_str = String::from("bad");

        let maybe_h = c.polymer_from_string(helix_str);
        let maybe_none = c.polymer_from_string(bad_str);

        match maybe_h {
            Some(_) => assert!(true),
            None => panic!("Recieved 'None' on good input"),
        }

        match maybe_none {
            None => assert!(true),
            Some(x) => panic!("Recieved 'Some({:?})' on bad input", x),
        }
    }

    #[test]
    fn helix_concat() {
        let c = DNACat::new();

        let str1 = String::from("at");
        let str2 = String::from("cg");
        let str3 = String::from("atcg");

        let mut h1 = c.polymer_from_string(str1).unwrap();
        let mut h2 = c.polymer_from_string(str2).unwrap();

        let h3 = c.polymer_from_string(str3).unwrap();
        h1.concat(&mut h2);

        assert_eq!(h1, h3)
    }

    #[test]
    fn test_inverse() {
        let c = DNACat::new();

        let str1 = String::from("atcg");
        let str2 = String::from("cgat");

        let h1 = c.polymer_from_string(str1).unwrap();
        let h2 = c.polymer_from_string(str2).unwrap();

        let i1 = c.inverse_p(&h1);
        let i2 = c.inverse_p(&i1);

        assert_eq!(i1, h2);
        assert_eq!(i2, h1);
    }

    #[test]
    fn test_pairs() {
        let c = DNACat::new();

        let str1 = String::from("atcg");
        let h1 = c.polymer_from_string(str1).unwrap();
        let pairs = c.pairs(h1);

        assert_eq!(
            vec![
                (c.a.clone(), c.t.clone()),
                (c.t.clone(), c.a.clone()),
                (c.c.clone(), c.g.clone()),
                (c.g.clone(), c.c.clone()),
            ],
            pairs
        );
    }

    #[test]
    fn test_strands() {
        let c = DNACat::new();

        let str1 = String::from("atcg");
        let str2 = String::from("cgat");

        let h1 = c.polymer_from_string(str1).unwrap();
        let h2 = c.polymer_from_string(str2).unwrap();

        let (fst1, snd1) = c.strands(h1);
        let (fst2, snd2) = c.strands(h2);

        assert_eq!(fst1, snd2);
        assert_eq!(fst2, snd1);
    }

    #[test]
    fn test_gc_content() {
        let c = DNACat::new();

        let str1 = String::from("atcg");
        let str2 = String::from("atata");

        let h1 = c.polymer_from_string(str1).unwrap();
        let h2 = c.polymer_from_string(str2).unwrap();

        let (maybe_2, maybe_4) = h1.gc_content();
        let (maybe_0, maybe_5) = h2.gc_content();
        assert_eq!(0, maybe_0);
        assert_eq!(2, maybe_2);
        assert_eq!(4, maybe_4);
        assert_eq!(5, maybe_5);
    }
}
