#[macro_use]
extern crate enum_display_derive;

use category::{Cat, ICat, NCat};
use monomer::{IMono, Mono, NucleicAcid};
use once_mono::{IMonomer};
use polymer::{Polymer, Strand};

use std::fmt::Display;

// TODO: Support the other RNA alphabets.
#[derive(Debug, Display, PartialEq, Clone)]
pub enum RNA {
    A,
    U,
    G,
    C,
}

impl NucleicAcid for RNA {
    fn is_g_or_c(&self) -> bool {
        match self {
            RNA::G => true,
            RNA::C => true,
            _ => false,
        }
    }
}

impl Mono for RNA {
    fn from_string(s: String) -> Option<RNA> {
        match s.to_uppercase().as_ref() {
            "A" => Some(RNA::A),
            "U" => Some(RNA::U),
            "C" => Some(RNA::C),
            "G" => Some(RNA::G),
            _ => None,
        }
    }
}

impl IMono for RNA {
    fn inverse(c: &Self) -> Self {
        match c {
            RNA::A => RNA::U,
            RNA::U => RNA::A,
            RNA::C => RNA::G,
            RNA::G => RNA::C,
        }
    }
}

pub type RNACell = IMonomer<RNA>;

pub struct RNACat {
    a: RNACell,
    u: RNACell,
    c: RNACell,
    g: RNACell,
}

impl Cat<RNACell, Strand<RNACell>> for RNACat {
    fn new() -> Self {
        RNACat {
            a: IMonomer::from_string(String::from("a")).unwrap(),
            u: IMonomer::from_string(String::from("u")).unwrap(),
            c: IMonomer::from_string(String::from("c")).unwrap(),
            g: IMonomer::from_string(String::from("g")).unwrap(),
        }
    }

    fn monomer_from_string(&self, s: String) -> Option<RNACell> {
        match RNA::from_string(s) {
            None => None,
            Some(x) => match x {
                RNA::A => Some(self.a.clone()),
                RNA::U => Some(self.u.clone()),
                RNA::C => Some(self.c.clone()),
                RNA::G => Some(self.g.clone()),
            },
        }
    }

    fn polymer_from_string(&self, s: String) -> Option<Strand<RNACell>> {
        let mut x = Strand::<RNACell>::new();
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

impl ICat<RNACell, Strand<RNACell>> for RNACat {
    fn inverse_m(&self, t: &RNACell) -> RNACell {
        match t.read().as_ref().unwrap() {
            RNA::A => self.u.clone(),
            RNA::U => self.a.clone(),
            RNA::C => self.g.clone(),
            RNA::G => self.c.clone(),
        }
    }

    fn inverse_p(&self, u: &Strand<RNACell>) -> Strand<RNACell> {
        let mut y = Strand::<RNACell>::new();
        for z in &u.contents {
            y.push(self.inverse_m(&z));
        }
        y.contents =  y.contents.into_iter().rev().collect();
        y
    }
}

impl NCat<RNACell, Strand<RNACell>> for RNACat {
    fn gc_content(u: Strand<RNACell>) -> (u64, u64) {
        let mut d = 0;
        let mut n = 0;
        for x in u.contents {
            d = d + 1;
            if RNACell::is_g_or_c(&x) {
                n = n + 1;
            }
        }

        (n, d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rna_from_char() {
        let c = RNACat::new();

        assert_eq!(
            c.monomer_from_string(String::from("A"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::A
        );
        assert_eq!(
            c.monomer_from_string(String::from("a"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::A
        );

        assert_eq!(
            c.monomer_from_string(String::from("U"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::U
        );

        assert_eq!(
            c.monomer_from_string(String::from("u"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::U
        );

        assert_eq!(
            c.monomer_from_string(String::from("C"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::C
        );
        assert_eq!(
            c.monomer_from_string(String::from("c"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::C
        );

        assert_eq!(
            c.monomer_from_string(String::from("G"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::G
        );

        assert_eq!(
            c.monomer_from_string(String::from("g"))
                .unwrap()
                .read()
                .as_ref()
                .unwrap(),
            &RNA::G
        );

        assert_eq!(c.monomer_from_string(String::from("d")), None);
    }

    // TODO Add Inverse:
    // Now for polynomial tests.

    #[test]
    fn new_strand_cat() {
        let c = RNACat::new();

        let strand_str = String::from("aucg");
        let h = c.polymer_from_string(strand_str).unwrap();
        let mut h2 = Strand::<RNACell>::new();

        h2.push(c.a.clone());
        h2.push(c.u.clone());
        h2.push(c.c.clone());
        h2.push(c.g.clone());

        assert_eq!(h, h2);
    }

    #[test]
    fn strand_from_string() {
        let c = RNACat::new();

        let strand_str = String::from("gauuaca");
        let bad_str = String::from("bad");

        let maybe_h = c.polymer_from_string(strand_str);
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
    fn strand_concat() {
        let c = RNACat::new();

        let str1 = String::from("au");
        let str2 = String::from("cg");
        let str3 = String::from("aucg");

        let mut h1 = c.polymer_from_string(str1).unwrap();
        let mut h2 = c.polymer_from_string(str2).unwrap();

        let h3 = c.polymer_from_string(str3).unwrap();
        h1.concat(&mut h2);

        assert_eq!(h1, h3)
    }

    #[test]
    fn test_inverse() {
        let c = RNACat::new();

        let str1 = String::from("aucg");
        let str2 = String::from("cgau");

        let h1 = c.polymer_from_string(str1).unwrap();
        let h2 = c.polymer_from_string(str2).unwrap();

        let i1 = c.inverse_p(&h1);
        let i2 = c.inverse_p(&i1);

        assert_eq!(i1, h2);
        assert_eq!(i2, h1);
    }

    #[test]
    fn test_gc_content() {
        let c = RNACat::new();

        let str1 = String::from("aucg");
        let str2 = String::from("auaua");

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