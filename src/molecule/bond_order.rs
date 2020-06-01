#[derive(Clone, Copy, Eq, Hash, PartialEq, Debug)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple
}

impl Default for BondOrder {
    fn default() -> Self {
        Self::Single
    }
}