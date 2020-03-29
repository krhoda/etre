#[macro_use]
extern crate enum_display_derive;

use once_cell::{OnceCell, OnceVal};
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

pub fn from_char(c: char) -> Option<Nucleotide> {
    match c.to_uppercase().to_string().as_ref() {
        "A" => Some(Nucleotide::A),
        "T" => Some(Nucleotide::T),
        "C" => Some(Nucleotide::C),
        "G" => Some(Nucleotide::G),
        _ => None,
    }
}

pub type Nucl = OnceVal<Nucleotide>;

#[derive(Debug)]
struct CatMachine {
    a: OnceCell<Nucleotide>,
    t: OnceCell<Nucleotide>,
    g: OnceCell<Nucleotide>,
    c: OnceCell<Nucleotide>,
}

// The below is a memory/speed(?) optimization I'm testing.
// It's also a mathematical-categorical view of DNA.
// For this reason, it has it's name which will likely change
// As it's structure becomes clearer.
#[derive(Clone, Debug)]
pub struct Cat(Arc<CatMachine>);

impl Cat {
    pub fn new() -> Cat {
        let mut a = OnceCell::<Nucleotide>::new();
        let mut t = OnceCell::<Nucleotide>::new();
        let mut g = OnceCell::<Nucleotide>::new();
        let mut c = OnceCell::<Nucleotide>::new();

        a.write(Nucleotide::A).unwrap();
        t.write(Nucleotide::T).unwrap();
        g.write(Nucleotide::G).unwrap();
        c.write(Nucleotide::C).unwrap();

        Cat(Arc::new(CatMachine {
            a: a,
            t: t,
            g: g,
            c: c,
        }))
    }

    pub fn A(&self) -> Nucl {
        self.0.a.read().unwrap()
    }

    pub fn T(&self) -> Nucl {
        self.0.t.read().unwrap()
    }

    pub fn C(&self) -> Nucl {
        self.0.c.read().unwrap()
    }

    pub fn G(&self) -> Nucl {
        self.0.g.read().unwrap()
    }

    pub fn compliment(&self, n: Nucl) -> Nucl {
        match n.read().as_ref().unwrap() {
            Nucleotide::A => self.T(),
            Nucleotide::T => self.A(),
            Nucleotide::C => self.G(),
            Nucleotide::G => self.C(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_cat_new() {
        let cat = Cat::new();

        assert_eq!(&Nucleotide::A, cat.A().read().as_ref().unwrap());
        assert_eq!(&Nucleotide::T, cat.T().read().as_ref().unwrap());
        assert_eq!(&Nucleotide::C, cat.C().read().as_ref().unwrap());
        assert_eq!(&Nucleotide::G, cat.G().read().as_ref().unwrap());
    }

    #[test]
    fn dna_complement() {
        let cat = Cat::new();
        let a = cat.A();
        let t = cat.T();
        let c = cat.C();
        let g = cat.G();

        assert_eq!(
            a.read().as_ref().unwrap(),
            cat.compliment(t.clone()).read().as_ref().unwrap()
        );
        assert_eq!(
            t.read().as_ref().unwrap(),
            cat.compliment(a.clone()).read().as_ref().unwrap()
        );
        assert_eq!(
            c.read().as_ref().unwrap(),
            cat.compliment(g.clone()).read().as_ref().unwrap()
        );
        assert_eq!(
            g.read().as_ref().unwrap(),
            cat.compliment(c.clone()).read().as_ref().unwrap()
        );
    }
}
