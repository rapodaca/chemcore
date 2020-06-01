#[derive(PartialEq, Eq, Debug)]
pub enum Error {
    HypervalentAtom,
    DuplicateBond,
    MisplacedBond,
    ImpossibleIsotope
}