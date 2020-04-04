#[macro_use]
extern crate enum_display_derive;

use once_cell::{OnceCell, OnceVal};

use category::Cat;
use monomer::Monomer;
use polymer::Polymer;

use std::fmt::Display;
use std::sync::Arc;

// TODO: Support the other RNA alphabets.
#[derive(Debug, Display, PartialEq)]
pub enum RNucleotide {
    A,
    U,
    G,
    C,
}

impl Monomer for RNucleotide {
    fn from_char(c: char) -> Option<RNucleotide> {
        match c.to_uppercase().to_string().as_ref() {
            "A" => Some(RNucleotide::A),
            "U" => Some(RNucleotide::U),
            "C" => Some(RNucleotide::C),
            "G" => Some(RNucleotide::G),
            _ => None,
        }
    }
}

pub type RNucl = OnceVal<RNucleotide>;

#[derive(Debug)]
pub struct Strand(Vec<RNucl>);

impl Polymer<RNucleotide> for Strand {
    fn new() -> Self {
        Strand(Vec::<RNucl>::new())
    }

    fn push(&mut self, n: RNucl) {
        self.0.push(n);
    }

    fn concat(&mut self, other: &mut Self) {
        self.0.append(&mut other.0);
    }
}

#[derive(Debug)]
struct CategoryRNAMachine {
    a: OnceVal<RNucleotide>,
    u: OnceVal<RNucleotide>,
    g: OnceVal<RNucleotide>,
    c: OnceVal<RNucleotide>,
}

#[derive(Clone, Debug)]
pub struct CategoryRNA(Arc<CategoryRNAMachine>);

impl CategoryRNA {
    pub fn compliment(&self, n: &RNucl) -> RNucl {
        match n.read().as_ref().unwrap() {
            RNucleotide::A => self.0.u.clone(),
            RNucleotide::U => self.0.a.clone(),
            RNucleotide::C => self.0.g.clone(),
            RNucleotide::G => self.0.c.clone(),
        }
    }

    pub fn inverse(&self, h: &Strand) -> Strand {
        let mut next = Strand::new();

        for x in &h.0 {
            next.push(self.compliment(&x.clone()))
        }

        next.0 = next.0.into_iter().rev().collect();
        next
    }

    // TODO: Better result than tuple fraction.
    pub fn gc_content(&self, h: &Strand) -> (u64, u64) {
        let mut n: u64 = 0;
        let mut d: u64 = 0;
        for x in &h.0 {
            d = d + 1;
            if CategoryRNA::is_g_or_c(x.read().as_ref().unwrap()) {
                n = n + 1;
            }
        }

        (n, d)
    }

    fn is_g_or_c(n: &RNucleotide) -> bool {
        match n {
            RNucleotide::G => true,
            RNucleotide::C => true,
            _ => false,
        }
    }
}

impl Cat<RNucleotide, Strand> for CategoryRNA {
    fn new() -> CategoryRNA {
        let mut ac = OnceCell::<RNucleotide>::new();
        let mut uc = OnceCell::<RNucleotide>::new();
        let mut gc = OnceCell::<RNucleotide>::new();
        let mut cc = OnceCell::<RNucleotide>::new();

        ac.write(RNucleotide::A).unwrap();
        uc.write(RNucleotide::U).unwrap();
        gc.write(RNucleotide::G).unwrap();
        cc.write(RNucleotide::C).unwrap();

        let a = ac.read().unwrap();
        let u = uc.read().unwrap();
        let g = gc.read().unwrap();
        let c = cc.read().unwrap();

        CategoryRNA(Arc::new(CategoryRNAMachine {
            a: a,
            u: u,
            g: g,
            c: c,
        }))
    }

    fn read(&self, n: RNucleotide) -> RNucl {
        match n {
            RNucleotide::A => self.0.a.clone(),
            RNucleotide::U => self.0.u.clone(),
            RNucleotide::C => self.0.c.clone(),
            RNucleotide::G => self.0.g.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: ADD READ TEST.

    #[test]
    fn rna_cat_new() {
        let cat = CategoryRNA::new();

        assert_eq!(&RNucleotide::A, cat.0.a.read().as_ref().unwrap());
        assert_eq!(&RNucleotide::U, cat.0.u.read().as_ref().unwrap());
        assert_eq!(&RNucleotide::C, cat.0.c.read().as_ref().unwrap());
        assert_eq!(&RNucleotide::G, cat.0.g.read().as_ref().unwrap());
    }

    #[test]
    fn dna_complement() {
        let cat = CategoryRNA::new();

        assert_eq!(
            cat.0.a.read().as_ref().unwrap(),
            cat.compliment(&cat.0.u.clone()).read().as_ref().unwrap()
        );
        assert_eq!(
            cat.0.u.read().as_ref().unwrap(),
            cat.compliment(&cat.0.a.clone()).read().as_ref().unwrap()
        );
        assert_eq!(
            cat.0.c.read().as_ref().unwrap(),
            cat.compliment(&cat.0.g.clone()).read().as_ref().unwrap()
        );
        assert_eq!(
            cat.0.g.read().as_ref().unwrap(),
            cat.compliment(&cat.0.c.clone()).read().as_ref().unwrap()
        );
    }

    #[test]
    fn rna_from_char() {
        assert_eq!(
            RNucleotide::from_char("A".chars().next().unwrap()),
            Some(RNucleotide::A)
        );
        assert_eq!(
            RNucleotide::from_char("a".chars().next().unwrap()),
            Some(RNucleotide::A)
        );

        assert_eq!(
            RNucleotide::from_char("U".chars().next().unwrap()),
            Some(RNucleotide::U)
        );
        assert_eq!(
            RNucleotide::from_char("u".chars().next().unwrap()),
            Some(RNucleotide::U)
        );

        assert_eq!(
            RNucleotide::from_char("C".chars().next().unwrap()),
            Some(RNucleotide::C)
        );
        assert_eq!(
            RNucleotide::from_char("c".chars().next().unwrap()),
            Some(RNucleotide::C)
        );

        assert_eq!(
            RNucleotide::from_char("G".chars().next().unwrap()),
            Some(RNucleotide::G)
        );
        assert_eq!(
            RNucleotide::from_char("g".chars().next().unwrap()),
            Some(RNucleotide::G)
        );

        assert_eq!(RNucleotide::from_char("d".chars().next().unwrap()), None);
    }

    #[test]
    fn strand_new() {
        let cat = CategoryRNA::new();
        let mut s = Strand::new();
        let mut tvec = Vec::<RNucl>::new();
        assert_eq!(&tvec, &s.0);
        s.push(cat.read(RNucleotide::A));
        tvec.push(cat.0.a.clone());

        assert_eq!(&tvec, &s.0);

        s.push(cat.read(RNucleotide::U));
        tvec.push(cat.0.u.clone());

        assert_eq!(&tvec, &s.0);

        s.push(cat.read(RNucleotide::C));
        tvec.push(cat.0.c.clone());

        assert_eq!(&tvec, &s.0);

        s.push(cat.read(RNucleotide::G));
        tvec.push(cat.0.g.clone());

        assert_eq!(&tvec, &s.0);
    }

    #[test]
    fn helix_from_string() {
        let cat = CategoryRNA::new();

        let strand = String::from("gauuaca");
        let badstrand = String::from("a123g");

        let mut control_strand = Strand::new();

        control_strand.push(cat.read(RNucleotide::G));
        control_strand.push(cat.read(RNucleotide::A));
        control_strand.push(cat.read(RNucleotide::U));
        control_strand.push(cat.read(RNucleotide::U));
        control_strand.push(cat.read(RNucleotide::A));
        control_strand.push(cat.read(RNucleotide::C));
        control_strand.push(cat.read(RNucleotide::A));

        let maybe_none = cat.from_string(badstrand);
        match maybe_none {
            None => assert!(true),
            Some(x) => panic!("Should've received nothing, got: {:?}", x),
        }

        let maybe_strand = cat.from_string(strand);
        match maybe_strand {
            None => panic!("Failed in Strand from_string with good string"),
            Some(h) => assert_eq!(h.0, control_strand.0),
        }
    }

    #[test]
    fn test_inverse_and_strands() {
        let cat = CategoryRNA::new();

        let strandstr = String::from("gauuaca");
        let invstr = String::from("uguaauc");

        let strand = cat.from_string(strandstr).unwrap();
        let strand2 = cat.from_string(invstr).unwrap();

        assert_eq!(cat.inverse(&strand).0, strand2.0);
        assert_eq!(cat.inverse(&cat.inverse(&strand)).0, strand.0);

    }

    #[test]
    fn test_gc_content() {
        let cat = CategoryRNA::new();

        let strandstr = String::from("gauuaca");
        let nogc = String::from("uuaauu");

        let strand = cat.from_string(strandstr).unwrap();
        let strand2 = cat.from_string(nogc).unwrap();
        let (should_be_2, should_be_7) = cat.gc_content(&strand);
        let (should_be_0, should_be_6) = cat.gc_content(&strand2);

        assert_eq!(should_be_0, 0);
        assert_eq!(should_be_2, 2);
        assert_eq!(should_be_6, 6);
        assert_eq!(should_be_7, 7);
    }
}
