use monomer::{IMono, Mono, NucleicAcid};

// TODO: See how to implement this thing as an iterator for read only
pub trait Polymer<T>: IntoIterator + std::marker::Sized + PartialEq
where
    T: Mono,
{
    fn new() -> Self;
    fn push(&mut self, t: T);
    fn concat(&mut self, other: &mut Self);

    fn from_string(s: String) -> Option<Self> {
        let mut x = Self::new();
        let mut ok = true;

        for c in s.chars() {
            match T::from_char(c) {
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

pub trait IPolymer<T>: Polymer<T>
where
    T: IMono,
{
    fn inverse(&self) -> Self;
}

// NOTE: No implementation for Mono.
// Means No protien strucutres as of yet.
#[derive(Debug)]
pub struct Strand<T>(Vec<T>)
where
    T: Mono;

impl<T: Mono> IntoIterator for Strand<T> {
    type Item = T;
    type IntoIter = std::vec::IntoIter<Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<T: Mono> PartialEq for Strand<T> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<T: Mono> Polymer<T> for Strand<T> {
    fn new() -> Self {
        Strand(Vec::<T>::new())
    }

    fn push(&mut self, t: T) {
        self.0.push(t);
    }

    fn concat(&mut self, other: &mut Self) {
        self.0.append(&mut other.0);
    }
}

impl<T: IMono> IPolymer<T> for Strand<T> {
    fn inverse(&self) -> Self {
        let mut next = Strand::new();

        for x in &self.0 {
            // NOTE: Using clone here, because I'm counting on OnceVals
            // worth remembering.
            next.push(T::inverse(&x.clone()));

            // TODO: TEST:
            // next.push(T::inverse(&x))
        }

        next.0 = next.0.into_iter().rev().collect();
        next
    }
}

#[derive(Debug)]
pub struct Helix<T>(Strand<T>)
where
    T: IMono;

impl<T: IMono> IntoIterator for Helix<T> {
    type Item = T;
    type IntoIter = std::vec::IntoIter<Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<T: IMono> PartialEq for Helix<T> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<T: IMono> Polymer<T> for Helix<T> {
    fn new() -> Self {
        Helix(Strand::<T>::new())
    }

    fn push(&mut self, t: T) {
        self.0.push(t);
    }

    fn concat(&mut self, other: &mut Self) {
        self.0.concat(&mut other.0);
    }
}

impl<T: IMono> IPolymer<T> for Helix<T> {
    fn inverse(&self) -> Self {
        let mut x = Helix::<T>::new();
        x.0 = self.0.inverse();
        x
    }
}

impl<T: IMono> Helix<T> {
    pub fn pairs(&self) -> Vec<(T, T)> {
        let mut next = Vec::<(T, T)>::new();
        let a = &self.0;
        for x in &a.0 {
            next.push((x.clone(), T::inverse(&x)));
        }

        next
    }

    pub fn strands(self) -> (Strand<T>, Strand<T>) {
        let x = self.inverse();
        (self.0, x.0)
    }
}

// TODO: MOVE TO STRAND!
impl<T: NucleicAcid> Helix<T> {
    pub fn gc_content(&self) -> (u64, u64) {
        let mut n = 0;
        let mut d = 0;

        let a = &self.0;
        for x in &a.0 {
            d = d + 1;
            if x.is_g_or_c() {
                n = n + 1;
            }
        }

        (n, d)
    }
}
