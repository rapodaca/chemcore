use purr::parts;

#[derive(Debug,PartialEq,Clone)]
pub enum Parity {
    Positive,
    Negative
}

impl Into<Parity> for &parts::Parity {
    fn into(self) -> Parity {
        match self {
            parts::Parity::Clockwise => Parity::Positive,
            parts::Parity::Counterclockwise => Parity::Negative,
            _ => unimplemented!("unsupported parity")
        }
    }
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