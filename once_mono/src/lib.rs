use monomer::{IMono, Mono, NucleicAcid};
use once_cell::{OnceCell, OnceVal};

pub fn create<T>(t: T) -> OnceVal<T>
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
    fn from_string(c: String) -> Option<Self> {
        let x = T::from_string(c);
        match x {
            Some(y) => Some(Monomer(create(y))),
            None => None,
        }
    }
}

impl<T: Mono> Monomer<T> {
    pub fn read(&self) -> std::sync::RwLockReadGuard<'_, Option<T>> {
        self.0.read()
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
    fn from_string(s: String) -> Option<Self> {
        let x = T::from_string(s);
        match x {
            Some(y) => Some(IMonomer(create(y))),
            None => None,
        }
    }
}

impl<T: IMono> IMono for IMonomer<T> {
    fn inverse(c: &Self) -> Self {
        IMonomer(create(T::inverse(c.0.read().as_ref().unwrap())))
    }
}

impl<T: NucleicAcid> NucleicAcid for IMonomer<T> {
    fn is_g_or_c(&self) -> bool {
        self.0.read().as_ref().unwrap().is_g_or_c()
    }
}

impl<T: IMono> IMonomer<T> {
    pub fn read(&self) -> std::sync::RwLockReadGuard<'_, Option<T>> {
        self.0.read()
    }
}
