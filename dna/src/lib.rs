#[macro_use]
extern crate enum_display_derive;
// use once_cell::{OnceCell, OnceVal};

// TODO RESTORE:
use category::{Cat, ICat, NCat};
use monomer::{IMono, Mono, NucleicAcid};
use once_mono::{IMonomer, Monomer};
use polymer::{Helix, Polymer, Strand};

use std::fmt::Display;
// use std::sync::Arc;

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
    fn from_char(c: char) -> Option<DNA> {
        match c.to_uppercase().to_string().as_ref() {
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
            a: IMonomer::from_char(String::from("a").chars().next().unwrap()).unwrap(),
            t: IMonomer::from_char(String::from("t").chars().next().unwrap()).unwrap(),
            c: IMonomer::from_char(String::from("c").chars().next().unwrap()).unwrap(),
            g: IMonomer::from_char(String::from("g").chars().next().unwrap()).unwrap(),
        }
    }

    fn from_char(&self, c: char) -> Option<DNACell> {
        match DNA::from_char(c) {
            None => None,
            Some(x) => match x {
                DNA::A => Some(self.a.clone()),
                DNA::T => Some(self.t.clone()),
                DNA::C => Some(self.c.clone()),
                DNA::G => Some(self.g.clone()),
            },
        }
    }

    fn from_string(&self, s: String) -> Option<Helix<DNACell>> {
        let mut x = Helix::<DNACell>::new();
        let mut ok = true;

        for c in s.chars() {
            match self.from_char(c) {
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
    use polymer::{Helix, IPolymer, Polymer};

    #[test]
    fn dna_from_char() {
        let c = DNACat::new();

        assert_eq!(
            c.from_char("A".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::A
        );
        assert_eq!(
            c.from_char("a".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::A
        );

        assert_eq!(
            c.from_char("T".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::T
        );

        assert_eq!(
            c.from_char("t".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::T
        );

        assert_eq!(
            c.from_char("C".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::C
        );
        assert_eq!(
            c.from_char("c".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::C
        );

        assert_eq!(
            c.from_char("G".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::G
        );

        assert_eq!(
            c.from_char("g".chars().next().unwrap())
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &DNA::G
        );

        assert_eq!(c.from_char("d".chars().next().unwrap()), None);
    }

    // TODO Add Inverse:

    // Now for polynomial tests.

    #[test]
    fn helix_new_push_eq() {
        let helix_str = String::from("atcg");
        let h = Helix::<DNA>::from_string(helix_str).unwrap();
        let mut h2 = Helix::new();

        h2.push(DNA::A);
        h2.push(DNA::T);
        h2.push(DNA::C);
        h2.push(DNA::G);

        assert_eq!(h, h2);
    }

    #[test]
    fn helix_from_string() {
        let helix_str = String::from("gattaca");
        let bad_str = String::from("bad");

        let maybe_h = Helix::<DNA>::from_string(helix_str);
        let maybe_none = Helix::<DNA>::from_string(bad_str);

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
        let str1 = String::from("at");
        let str2 = String::from("cg");
        let str3 = String::from("atcg");

        let mut h1 = Helix::<DNA>::from_string(str1).unwrap();
        let mut h2 = Helix::<DNA>::from_string(str2).unwrap();
        let h3 = Helix::<DNA>::from_string(str3).unwrap();
        h1.concat(&mut h2);
        assert_eq!(h1, h3)
    }

    #[test]
    fn test_inverse() {
        let str1 = String::from("atcg");
        let str2 = String::from("cgat");

        let h1 = Helix::<DNA>::from_string(str1).unwrap();
        let h2 = Helix::<DNA>::from_string(str2).unwrap();

        let i1 = h1.inverse();
        let i2 = i1.inverse();
        assert_eq!(i1, h2);
        assert_eq!(i2, h1);
    }

    #[test]
    fn test_pairs() {
        let str1 = String::from("atcg");
        let h1 = Helix::<DNA>::from_string(str1).unwrap();
        let pairs = h1.pairs();

        assert_eq!(
            vec![
                (DNA::A, DNA::T),
                (DNA::T, DNA::A),
                (DNA::C, DNA::G),
                (DNA::G, DNA::C),
            ],
            pairs
        );
    }

    #[test]
    fn test_strands() {
        let str1 = String::from("atcg");
        let str2 = String::from("cgat");

        let h1 = Helix::<DNA>::from_string(str1).unwrap();
        let h2 = Helix::<DNA>::from_string(str2).unwrap();

        let (fst1, snd1) = h1.strands();
        let (fst2, snd2) = h2.strands();

        assert_eq!(fst1, snd2);
        assert_eq!(fst2, snd1);
    }

    #[test]
    fn test_gc_content() {
        let str1 = String::from("atcg");
        let str2 = String::from("atata");

        let h1 = Helix::<DNA>::from_string(str1).unwrap();
        let h2 = Helix::<DNA>::from_string(str2).unwrap();

        let (maybe_2, maybe_4) = h1.gc_content();
        let (maybe_0, maybe_5) = h2.gc_content();
        assert_eq!(0, maybe_0);
        assert_eq!(2, maybe_2);
        assert_eq!(4, maybe_4);
        assert_eq!(5, maybe_5);
    }
}
