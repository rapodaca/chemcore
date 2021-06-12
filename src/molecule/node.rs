use super::{Atom, Bond};

#[derive(Debug, PartialEq)]
pub struct Node {
    pub atom: Atom,
    pub bonds: Vec<Bond>,
}
