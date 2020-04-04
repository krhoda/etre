pub trait Monomer {
    fn from_char(c: char) -> Option<Self> where Self: std::marker::Sized;
}

pub trait IMonomer: Monomer {
    fn inverse(c: Self) -> Self;
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
