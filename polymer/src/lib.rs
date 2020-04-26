use monomer::{IMono, Mono, NucleicAcid};

// TODO: See how to implement this thing as an iterator for read only
pub trait Polymer<T>: std::marker::Sized + PartialEq
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
            match T::from_string(c.to_string()) {
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
pub struct Strand<T>
where
    T: Mono,
{
    pub contents: Vec<T>,
}

impl<T: Mono> PartialEq for Strand<T> {
    fn eq(&self, other: &Self) -> bool {
        self.contents == other.contents
    }
}

impl<T: Mono> Polymer<T> for Strand<T> {
    fn new() -> Self {
        Strand {
            contents: Vec::<T>::new(),
        }
    }

    fn push(&mut self, t: T) {
        self.contents.push(t);
    }

    fn concat(&mut self, other: &mut Self) {
        self.contents.append(&mut other.contents);
    }
}

impl<T: IMono> IPolymer<T> for Strand<T> {
    fn inverse(&self) -> Self {
        let mut next = Strand::new();

        for x in &self.contents {
            // NOTE: Using clone here, because I'm counting on OnceVals
            // worth remembering.
            next.push(T::inverse(&x));

            // TODO: TEST:
            // next.push(T::inverse(&x))
        }

        next.contents = next.contents.into_iter().rev().collect();
        next
    }
}

impl<T: NucleicAcid> Strand<T> {
    pub fn gc_content(&self) -> (u64, u64) {
        let mut n = 0;
        let mut d = 0;

        for x in &self.contents {
            d = d + 1;
            if x.is_g_or_c() {
                n = n + 1;
            }
        }

        (n, d)
    }
}

#[derive(Debug)]
pub struct Helix<T>
where
    T: IMono,
{
    pub strand: Strand<T>,
}

impl<T: IMono> PartialEq for Helix<T> {
    fn eq(&self, other: &Self) -> bool {
        self.strand == other.strand
    }
}

impl<T: IMono> Polymer<T> for Helix<T> {
    fn new() -> Self {
        Helix {
            strand: Strand::<T>::new(),
        }
    }

    fn push(&mut self, t: T) {
        self.strand.push(t);
    }

    fn concat(&mut self, other: &mut Self) {
        self.strand.concat(&mut other.strand);
    }
}

impl<T: IMono> IPolymer<T> for Helix<T> {
    fn inverse(&self) -> Self {
        let mut x = Helix::<T>::new();
        x.strand = self.strand.inverse();
        x
    }
}

impl<T: IMono> Helix<T> {
    pub fn pairs(&self) -> Vec<(T, T)> {
        let mut next = Vec::<(T, T)>::new();
        let a = &self.strand;
        for x in &a.contents {
            next.push((x.clone(), T::inverse(&x)));
        }

        next
    }

    pub fn strands(self) -> (Strand<T>, Strand<T>) {
        let x = self.inverse();
        (self.strand, x.strand)
    }
}

// TODO: MOVE TO STRAND!
impl<T: NucleicAcid> Helix<T> {
    pub fn gc_content(&self) -> (u64, u64) {
        self.strand.gc_content()
    }
}
