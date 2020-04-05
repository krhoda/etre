use monomer::{IMono, Mono, NucleicAcid};
use polymer::{IPolymer, Polymer};

pub trait Cat<T, U>
where
    T: Mono,
    U: Polymer<T>,
{
    fn new() -> Self;
    fn from_char(&self, c: char) -> Option<T>;
    fn from_string(&self, s: String) -> Option<U>;
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