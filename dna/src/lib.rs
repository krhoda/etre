#[macro_use]
extern crate enum_display_derive;

use once_cell::{OnceCell, OnceVal};

use category::Cat;
use monomer::Monomer;
use polymer::Polymer;

use std::fmt::Display;
use std::sync::Arc;

// TODO: Support the other DNA alphabets.

#[derive(Debug, Display, PartialEq)]
pub enum DNucleotide {
    A,
    T,
    G,
    C,
}

impl Monomer for DNucleotide {
    fn from_char(c: char) -> Option<DNucleotide> {
        match c.to_uppercase().to_string().as_ref() {
            "A" => Some(DNucleotide::A),
            "T" => Some(DNucleotide::T),
            "C" => Some(DNucleotide::C),
            "G" => Some(DNucleotide::G),
            _ => None,
        }
    }
}

pub type DNucl = OnceVal<DNucleotide>;

#[derive(Debug)]
pub struct Helix(Vec<DNucl>);

impl Polymer<DNucleotide> for Helix {
    fn new() -> Self {
        Helix(Vec::<DNucl>::new())
    }

    fn push(&mut self, n: DNucl) {
        self.0.push(n);
    }

    fn concat(&mut self, other: &mut Self) {
        self.0.append(&mut other.0);
    }
}

#[derive(Debug)]
struct CategoryDNAMachine {
    a: OnceVal<DNucleotide>,
    t: OnceVal<DNucleotide>,
    g: OnceVal<DNucleotide>,
    c: OnceVal<DNucleotide>,
}

// The below is a memory/speed(?) optimization I'm testing.
// It's also a mathematical-categorical view of DNA.
// For this reason, it has it's name which will likely change
// As it's structure becomes clearer.
#[derive(Clone, Debug)]
pub struct CategoryDNA(Arc<CategoryDNAMachine>);

impl CategoryDNA {
    pub fn compliment(&self, n: &DNucl) -> DNucl {
        match n.read().as_ref().unwrap() {
            DNucleotide::A => self.0.t.clone(),
            DNucleotide::T => self.0.a.clone(),
            DNucleotide::C => self.0.g.clone(),
            DNucleotide::G => self.0.c.clone(),
        }
    }

    pub fn inverse(&self, h: &Helix) -> Helix {
        let mut next = Helix::new();

        for x in &h.0 {
            next.push(self.compliment(&x.clone()))
        }

        next.0 = next.0.into_iter().rev().collect();
        next
    }

    pub fn pairs(&self, h: &Helix) -> Vec<(DNucl, DNucl)> {
        let mut next = Vec::<(DNucl, DNucl)>::new();
        for x in &h.0 {
            next.push((x.clone(), self.compliment(&x.clone())));
        }

        next
    }

    pub fn strands(&self, h: Helix) -> (Helix, Helix) {
        let x = self.inverse(&h);
        (h, x)
    }

    // TODO: Better result than tuple fraction.
    pub fn gc_content(&self, h: &Helix) -> (u64, u64) {
        let mut n: u64 = 0;
        let mut d: u64 = 0;
        for x in &h.0 {
            d = d + 1;
            if CategoryDNA::is_g_or_c(x.read().as_ref().unwrap()) {
                n = n + 1;
            }
        }

        (n, d)
    }

    fn is_g_or_c(n: &DNucleotide) -> bool {
        match n {
            DNucleotide::G => true,
            DNucleotide::C => true,
            _ => false,
        }
    }
}

impl Cat<DNucleotide, Helix> for CategoryDNA {
    fn new() -> CategoryDNA {
        let mut ac = OnceCell::<DNucleotide>::new();
        let mut tc = OnceCell::<DNucleotide>::new();
        let mut gc = OnceCell::<DNucleotide>::new();
        let mut cc = OnceCell::<DNucleotide>::new();

        ac.write(DNucleotide::A).unwrap();
        tc.write(DNucleotide::T).unwrap();
        gc.write(DNucleotide::G).unwrap();
        cc.write(DNucleotide::C).unwrap();

        let a = ac.read().unwrap();
        let t = tc.read().unwrap();
        let g = gc.read().unwrap();
        let c = cc.read().unwrap();

        CategoryDNA(Arc::new(CategoryDNAMachine {
            a: a,
            t: t,
            g: g,
            c: c,
        }))
    }

    fn read(&self, n: DNucleotide) -> DNucl {
        match n {
            DNucleotide::A => self.0.a.clone(),
            DNucleotide::T => self.0.t.clone(),
            DNucleotide::C => self.0.c.clone(),
            DNucleotide::G => self.0.g.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: ADD READ TEST.

    #[test]
    fn dna_cat_new() {
        let cat = CategoryDNA::new();

        assert_eq!(&DNucleotide::A, cat.0.a.read().as_ref().unwrap());
        assert_eq!(&DNucleotide::T, cat.0.t.read().as_ref().unwrap());
        assert_eq!(&DNucleotide::C, cat.0.c.read().as_ref().unwrap());
        assert_eq!(&DNucleotide::G, cat.0.g.read().as_ref().unwrap());
    }

    #[test]
    fn dna_complement() {
        let cat = CategoryDNA::new();

        assert_eq!(
            cat.0.a.read().as_ref().unwrap(),
            cat.compliment(&cat.0.t.clone()).read().as_ref().unwrap()
        );
        assert_eq!(
            cat.0.t.read().as_ref().unwrap(),
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
    fn dna_from_char() {
        assert_eq!(
            DNucleotide::from_char("A".chars().next().unwrap()),
            Some(DNucleotide::A)
        );
        assert_eq!(
            DNucleotide::from_char("a".chars().next().unwrap()),
            Some(DNucleotide::A)
        );

        assert_eq!(
            DNucleotide::from_char("T".chars().next().unwrap()),
            Some(DNucleotide::T)
        );
        assert_eq!(
            DNucleotide::from_char("t".chars().next().unwrap()),
            Some(DNucleotide::T)
        );

        assert_eq!(
            DNucleotide::from_char("C".chars().next().unwrap()),
            Some(DNucleotide::C)
        );
        assert_eq!(
            DNucleotide::from_char("c".chars().next().unwrap()),
            Some(DNucleotide::C)
        );

        assert_eq!(
            DNucleotide::from_char("G".chars().next().unwrap()),
            Some(DNucleotide::G)
        );
        assert_eq!(
            DNucleotide::from_char("g".chars().next().unwrap()),
            Some(DNucleotide::G)
        );

        assert_eq!(DNucleotide::from_char("d".chars().next().unwrap()), None);
    }

    #[test]
    fn helix_new() {
        let cat = CategoryDNA::new();
        let mut h = Helix::new();
        let mut tvec = Vec::<DNucl>::new();
        assert_eq!(&tvec, &h.0);
        h.push(cat.read(DNucleotide::A));
        tvec.push(cat.0.a.clone());

        assert_eq!(&tvec, &h.0);

        h.push(cat.read(DNucleotide::T));
        tvec.push(cat.0.t.clone());

        assert_eq!(&tvec, &h.0);

        h.push(cat.read(DNucleotide::C));
        tvec.push(cat.0.c.clone());

        assert_eq!(&tvec, &h.0);

        h.push(cat.read(DNucleotide::G));
        tvec.push(cat.0.g.clone());

        assert_eq!(&tvec, &h.0);
    }

    #[test]
    fn helix_from_string() {
        let cat = CategoryDNA::new();

        let hxstr = String::from("gattaca");
        let badstr = String::from("a123g");

        let mut control_h = Helix::new();

        control_h.push(cat.read(DNucleotide::G));
        control_h.push(cat.read(DNucleotide::A));
        control_h.push(cat.read(DNucleotide::T));
        control_h.push(cat.read(DNucleotide::T));
        control_h.push(cat.read(DNucleotide::A));
        control_h.push(cat.read(DNucleotide::C));
        control_h.push(cat.read(DNucleotide::A));

        let maybe_none = cat.from_string(badstr);
        match maybe_none {
            None => assert!(true),
            Some(x) => panic!("Should've recieved nothing, got: {:?}", x),
        }

        let maybe_h = cat.from_string(hxstr);
        match maybe_h {
            None => panic!("Failed in Helix from_string with good string"),
            Some(h) => assert_eq!(h.0, control_h.0),
        }
    }

    #[test]
    fn test_inverse_and_strands() {
        let cat = CategoryDNA::new();

        let hxstr = String::from("gattaca");
        let invhxstr = String::from("tgtaatc");
        let h = cat.from_string(hxstr).unwrap();
        let h2 = cat.from_string(invhxstr).unwrap();
        assert_eq!(cat.inverse(&h).0, h2.0);
        assert_eq!(cat.inverse(&cat.inverse(&h)).0, h.0);

        let (fst, snd) = cat.strands(h);
        assert_eq!(fst.0, cat.inverse(&h2).0);
        assert_eq!(snd.0, h2.0);
    }

    #[test]
    fn test_gc_content() {
        let cat = CategoryDNA::new();

        let hxstr = String::from("gattaca");
        let nogc = String::from("ttaatt");

        let h = cat.from_string(hxstr).unwrap();
        let h2 = cat.from_string(nogc).unwrap();
        let (should_be_2, should_be_7) = cat.gc_content(&h);
        let (should_be_0, should_be_6) = cat.gc_content(&h2);

        assert_eq!(should_be_0, 0);
        assert_eq!(should_be_2, 2);
        assert_eq!(should_be_6, 6);
        assert_eq!(should_be_7, 7);
    }
}
