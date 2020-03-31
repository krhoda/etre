use monomer::Monomer;
use once_cell::OnceVal;
use polymer::Polymer;

pub trait Cat<T, U>
where
    T: Monomer + PartialEq,
    U: Polymer<T>,
{
    fn new() -> Self;
    fn read(&self, t: T) -> OnceVal<T>;

    // NOTE: Ultra optimization here would be to demand a direct from_char implementation here.
    fn from_char(&self, c: char) -> Option<OnceVal<T>> {
        match T::from_char(c) {
            Some(x) => Some(self.read(x)),
            None => None,
        }
    }

    fn from_string(&self, s: String) -> Option<U>
    where
        U: std::marker::Sized,
    {
        let mut x = U::new();
        let mut ok = true;
        for c in s.chars() {
            match self.from_char(c) {
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
