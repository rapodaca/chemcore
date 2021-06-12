use purr::read::Error as PurrError;

#[derive(Debug, PartialEq)]
pub enum Error {
    Character(usize),
    BondKind(usize),
    Valence(usize),
    Isotope(usize),
    Parity(usize),
    ChargedStar(usize),
    IncompatibleJoin(usize, usize),
    Kekulization,
    EndOfLine,
}

impl From<PurrError> for Error {
    fn from(error: PurrError) -> Self {
        match error {
            PurrError::InvalidCharacter(cursor) => Error::Character(cursor),
            PurrError::EndOfLine => Error::EndOfLine,
        }
    }
}
