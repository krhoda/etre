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

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
