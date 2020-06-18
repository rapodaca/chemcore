use purr::read::Error as PurrError;

// See: https://doc.rust-lang.org/stable/rust-by-example/error/multiple_error_types/wrap_error.html

#[derive(PartialEq, Debug)]
pub enum Error {
    HypervalentAtom,
    DuplicateBond,
    MisplacedBond,
    ImpossibleIsotope,
    CanNotKekulize,
    Purr(PurrError)
}

impl std::convert::From<PurrError> for Error {
    fn from(err: PurrError) -> Error {
        Error::Purr(err)
    }
}