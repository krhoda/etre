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
pub enum Nucleotide {
    A,
    T,
    G,
    C,
}

impl Monomer for Nucleotide {
    fn from_char(c: char) -> Option<Nucleotide> {
        match c.to_uppercase().to_string().as_ref() {
            "A" => Some(Nucleotide::A),
            "T" => Some(Nucleotide::T),
            "C" => Some(Nucleotide::C),
            "G" => Some(Nucleotide::G),
            _ => None,
        }
    }
}

pub type Nucl = OnceVal<Nucleotide>;

#[derive(Debug)]
pub struct Helix(Vec<Nucl>);

impl Polymer<Nucleotide> for Helix {
    fn new() -> Self {
        Helix(Vec::<Nucl>::new())
    }

    fn push(&mut self, n: Nucl) {
        self.0.push(n);
    }

    fn concat(&mut self, other: &mut Self) {
        self.0.append(&mut other.0);
    }
}

#[derive(Debug)]
struct CategoryDNAMachine {
    a: OnceVal<Nucleotide>,
    t: OnceVal<Nucleotide>,
    g: OnceVal<Nucleotide>,
    c: OnceVal<Nucleotide>,
}

// The below is a memory/speed(?) optimization I'm testing.
// It's also a mathematical-categorical view of DNA.
// For this reason, it has it's name which will likely change
// As it's structure becomes clearer.
#[derive(Clone, Debug)]
pub struct CategoryDNA(Arc<CategoryDNAMachine>);

impl CategoryDNA {
    pub fn compliment(&self, n: &Nucl) -> Nucl {
        match n.read().as_ref().unwrap() {
            Nucleotide::A => self.0.t.clone(),
            Nucleotide::T => self.0.a.clone(),
            Nucleotide::C => self.0.g.clone(),
            Nucleotide::G => self.0.c.clone(),
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

    pub fn pairs(&self, h: &Helix) -> Vec<(Nucl, Nucl)> {
        let mut next = Vec::<(Nucl, Nucl)>::new();
        for x in &h.0 {
            next.push((x.clone(), self.compliment(&x.clone())));
        }

        next
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

    fn is_g_or_c(n: &Nucleotide) -> bool {
        match n {
            Nucleotide::G => true,
            Nucleotide::C => true,
            _ => false,
        }
    }
}

impl Cat<Nucleotide, Helix> for CategoryDNA {
    fn new() -> CategoryDNA {
        let mut ac = OnceCell::<Nucleotide>::new();
        let mut tc = OnceCell::<Nucleotide>::new();
        let mut gc = OnceCell::<Nucleotide>::new();
        let mut cc = OnceCell::<Nucleotide>::new();

        ac.write(Nucleotide::A).unwrap();
        tc.write(Nucleotide::T).unwrap();
        gc.write(Nucleotide::G).unwrap();
        cc.write(Nucleotide::C).unwrap();

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

    fn read(&self, n: Nucleotide) -> Nucl {
        match n {
            Nucleotide::A => self.0.a.clone(),
            Nucleotide::T => self.0.t.clone(),
            Nucleotide::C => self.0.c.clone(),
            Nucleotide::G => self.0.g.clone(),
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

        assert_eq!(&Nucleotide::A, cat.0.a.read().as_ref().unwrap());
        assert_eq!(&Nucleotide::T, cat.0.t.read().as_ref().unwrap());
        assert_eq!(&Nucleotide::C, cat.0.c.read().as_ref().unwrap());
        assert_eq!(&Nucleotide::G, cat.0.g.read().as_ref().unwrap());
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
            Nucleotide::from_char("A".chars().next().unwrap()),
            Some(Nucleotide::A)
        );
        assert_eq!(
            Nucleotide::from_char("a".chars().next().unwrap()),
            Some(Nucleotide::A)
        );

        assert_eq!(
            Nucleotide::from_char("T".chars().next().unwrap()),
            Some(Nucleotide::T)
        );
        assert_eq!(
            Nucleotide::from_char("t".chars().next().unwrap()),
            Some(Nucleotide::T)
        );

        assert_eq!(
            Nucleotide::from_char("C".chars().next().unwrap()),
            Some(Nucleotide::C)
        );
        assert_eq!(
            Nucleotide::from_char("c".chars().next().unwrap()),
            Some(Nucleotide::C)
        );

        assert_eq!(
            Nucleotide::from_char("G".chars().next().unwrap()),
            Some(Nucleotide::G)
        );
        assert_eq!(
            Nucleotide::from_char("g".chars().next().unwrap()),
            Some(Nucleotide::G)
        );

        assert_eq!(Nucleotide::from_char("d".chars().next().unwrap()), None);
    }

    #[test]
    fn helix_new() {
        let cat = CategoryDNA::new();
        let mut h = Helix::new();
        let mut tvec = Vec::<Nucl>::new();
        assert_eq!(&tvec, &h.0);
        h.push(cat.read(Nucleotide::A));
        tvec.push(cat.0.a.clone());

        assert_eq!(&tvec, &h.0);

        h.push(cat.read(Nucleotide::T));
        tvec.push(cat.0.t.clone());

        assert_eq!(&tvec, &h.0);

        h.push(cat.read(Nucleotide::C));
        tvec.push(cat.0.c.clone());

        assert_eq!(&tvec, &h.0);

        h.push(cat.read(Nucleotide::G));
        tvec.push(cat.0.g.clone());

        assert_eq!(&tvec, &h.0);
    }

    #[test]
    fn helix_from_string() {
        let cat = CategoryDNA::new();

        let hxstr = String::from("gattaca");
        let badstr = String::from("a123g");

        let mut control_h = Helix::new();

        control_h.push(cat.read(Nucleotide::G));
        control_h.push(cat.read(Nucleotide::A));
        control_h.push(cat.read(Nucleotide::T));
        control_h.push(cat.read(Nucleotide::T));
        control_h.push(cat.read(Nucleotide::A));
        control_h.push(cat.read(Nucleotide::C));
        control_h.push(cat.read(Nucleotide::A));

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
}
