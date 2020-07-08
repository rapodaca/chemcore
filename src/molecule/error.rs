use purr::read::Error as PurrError;

#[derive(Debug, PartialEq)]
pub enum Error {
    Hypervalent(usize),
    ImpossibleIsotope(usize),
    AtomParityNotAllowed(usize),
    UnsupportedBond(usize, usize),
    InvalidConformation(usize, usize),
    CanNotKekulize,
    Purr(PurrError)
}

impl std::convert::From<PurrError> for Error {
    fn from(err: PurrError) -> Error {
        Error::Purr(err)
    }
}