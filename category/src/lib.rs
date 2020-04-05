use monomer::{IMono, Mono, NucleicAcid};
use once_cell::{OnceCell, OnceVal};
use polymer::{IPolymer, Polymer};

pub trait Cat<T, U>
where
    T: Mono,
    U: Polymer<T>,
{
    fn new() -> Self;
    // fn new_poly() -> U;
    fn from_char(&self, c: char) -> Option<T>;

    fn from_string(&self, s: String) -> Option<U>;
    // fn push(&self, t: T, u: &mut U);
    // fn concat(&self, fst: &mut U, snd: &mut U);
}

pub trait ICat<T, U>: Cat<T, U>
where
    T: IMono,
    U: IPolymer<T>,
{
    fn inverse_m(&self, t: &T) -> T;
    fn inverse_p(&self, u: &U) -> U;
}

pub trait NCat<T, U>: ICat<T, U>
where
    T: NucleicAcid,
    U: IPolymer<T>,
{
    fn gc_content(u: U) -> (u64, u64);
}


// pub trait Cat<T, U>
// where
//     T: Monomer + PartialEq,
//     U: Polymer<T>,
// {
//     fn new() -> Self;
//     fn read(&self, t: T) -> OnceVal<T>;

//     // NOTE: Ultra optimization here would be to demand a direct from_char implementation here.
//     fn from_char(&self, c: char) -> Option<OnceVal<T>> {
//         match T::from_char(c) {
//             Some(x) => Some(self.read(x)),
//             None => None,
//         }
//     }

//     fn from_string(&self, s: String) -> Option<U>
//     where
//         U: std::marker::Sized,
//     {
//         let mut x = U::new();
//         let mut ok = true;
//         for c in s.chars() {
//             match self.from_char(c) {
//                 Some(y) => x.push(y),
//                 None => {
//                     ok = false;
//                     break;
//                 }
//             }
//         }

//         match ok {
//             true => Some(x),
//             _ => None,
//         }
//     }
// }

// ---- FROM EARLIER:
// #[derive(Debug)]
// struct CategoryDNAMachine {
//     a: OnceVal<DNucleotide>,
//     t: OnceVal<DNucleotide>,
//     g: OnceVal<DNucleotide>,
//     c: OnceVal<DNucleotide>,
// }

// The below is a memory/speed(?) optimization I'm testing.
// It's also a mathematical-categorical view of DNA.
// For this reason, it has it's name which will likely change
// As it's structure becomes clearer.
// #[derive(Clone, Debug)]
// pub struct CategoryDNA(Arc<CategoryDNAMachine>);

// impl CategoryDNA {
//     pub fn compliment(&self, n: &DNucl) -> DNucl {
//         match n.read().as_ref().unwrap() {
//             DNucleotide::A => self.0.t.clone(),
//             DNucleotide::T => self.0.a.clone(),
//             DNucleotide::C => self.0.g.clone(),
//             DNucleotide::G => self.0.c.clone(),
//         }
//     }

//     pub fn inverse(&self, h: &Helix) -> Helix {
//         let mut next = Helix::new();

//         for x in &h.0 {
//             next.push(self.compliment(&x.clone()))
//         }

//         next.0 = next.0.into_iter().rev().collect();
//         next
//     }

//     pub fn pairs(&self, h: &Helix) -> Vec<(DNucl, DNucl)> {
//         let mut next = Vec::<(DNucl, DNucl)>::new();
//         for x in &h.0 {
//             next.push((x.clone(), self.compliment(&x.clone())));
//         }

//         next
//     }

//     pub fn strands(&self, h: Helix) -> (Helix, Helix) {
//         let x = self.inverse(&h);
//         (h, x)
//     }

//     // TODO: Better result than tuple fraction.
//     pub fn gc_content(&self, h: &Helix) -> (u64, u64) {
//         let mut n: u64 = 0;
//         let mut d: u64 = 0;
//         for x in &h.0 {
//             d = d + 1;
//             if CategoryDNA::is_g_or_c(x.read().as_ref().unwrap()) {
//                 n = n + 1;
//             }
//         }

//         (n, d)
//     }

//     fn is_g_or_c(n: &DNucleotide) -> bool {
//         match n {
//             DNucleotide::G => true,
//             DNucleotide::C => true,
//             _ => false,
//         }
//     }
// }

// impl Cat<DNucleotide, Helix> for CategoryDNA {
//     fn new() -> CategoryDNA {
//         let mut ac = OnceCell::<DNucleotide>::new();
//         let mut tc = OnceCell::<DNucleotide>::new();
//         let mut gc = OnceCell::<DNucleotide>::new();
//         let mut cc = OnceCell::<DNucleotide>::new();

//         ac.write(DNucleotide::A).unwrap();
//         tc.write(DNucleotide::T).unwrap();
//         gc.write(DNucleotide::G).unwrap();
//         cc.write(DNucleotide::C).unwrap();

//         let a = ac.read().unwrap();
//         let t = tc.read().unwrap();
//         let g = gc.read().unwrap();
//         let c = cc.read().unwrap();

//         CategoryDNA(Arc::new(CategoryDNAMachine {
//             a: a,
//             t: t,
//             g: g,
//             c: c,
//         }))
//     }

//     fn read(&self, n: DNucleotide) -> DNucl {
//         match n {
//             DNucleotide::A => self.0.a.clone(),
//             DNucleotide::T => self.0.t.clone(),
//             DNucleotide::C => self.0.c.clone(),
//             DNucleotide::G => self.0.g.clone(),
//         }
//     }
// }
