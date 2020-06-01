#[derive(Clone, Copy, Eq, Hash, PartialEq, Debug)]
pub enum Parity {
    Positive,
    Negative
}

impl Parity {
    pub fn negate(&self) -> Parity {
        match self {
            Parity::Positive => Parity::Negative,
            Parity::Negative => Parity::Positive
        }
    }

    pub fn multiply(&self, other: &Parity) -> Self {
        if self == other {
            Self::Positive
        } else {
            Self::Negative
        }
    }
}