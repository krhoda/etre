use once_cell::{OnceCell, OnceVal};

pub trait Mono: std::marker::Sized + PartialEq + Clone {
    fn from_char(c: char) -> Option<Self>;
}

pub trait IMono: Mono {
    fn inverse(c: &Self) -> Self;
}

pub trait NucleicAcid: IMono {
    fn is_g_or_c(&self) -> bool;
}

fn wrap<T>(t: T) -> OnceVal<T>
where
    T: Mono,
{
    let mut x = OnceCell::<T>::new();
    x.write(t).unwrap();
    x.read().unwrap()
}

#[derive(Debug)]
pub struct Monomer<T: Mono>(OnceVal<T>);
impl<T: Mono> PartialEq for Monomer<T> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<T: Mono> Clone for Monomer<T> {
    fn clone(&self) -> Self {
        Monomer(self.0.clone())
    }
}

impl<T: Mono> Mono for Monomer<T> {
    fn from_char(c: char) -> Option<Self> {
        let x = T::from_char(c);
        match x {
            Some(y) => Some(Monomer(wrap(y))),
            None => None,
        }
    }
}

#[derive(Debug)]
pub struct IMonomer<T: IMono>(OnceVal<T>);

impl<T: IMono> PartialEq for IMonomer<T> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<T: IMono> Clone for IMonomer<T> {
    fn clone(&self) -> Self {
        IMonomer(self.0.clone())
    }
}

impl<T: IMono> Mono for IMonomer<T> {
    fn from_char(c: char) -> Option<Self> {
        let x = T::from_char(c);
        match x {
            Some(y) => Some(IMonomer(wrap(y))),
            None => None,
        }
    }
}

impl<T: IMono> IMono for IMonomer<T> {
    fn inverse(c: &Self) -> Self {
        IMonomer(wrap(T::inverse(c.0.read().as_ref().unwrap())))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
