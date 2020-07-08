use purr::mol::{ Atom, Parity as PurrParity };
use crate::molecule::{ Parity };

pub fn atom_parity(id: usize, atom: &Atom) -> Option<Parity> {
    let mut result = match atom.nub.parity {
        Some(PurrParity::Counterclockwise) => Parity::Negative,
        Some(PurrParity::Clockwise) => Parity::Positive,
        None => return None
    };

    if let Some(hcount) = atom.nub.hcount {
        if id > 0 && hcount > 0 {
            result = result.negate();
        }
    }

    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use purr::read::read;

    #[test]
    fn none() {
        let atoms = read(&"CCC").unwrap();
        let parity = atom_parity(1, &atoms[0]);

        assert_eq!(parity, None);
    }

    #[test]
    fn counterclockwise_leading_with_hydrogen() {
        let atoms = read(&"[C@H](N)(O)C").unwrap();
        let parity = atom_parity(0, &atoms[0]);

        assert_eq!(parity, Some(Parity::Negative));
    }

    #[test]
    fn clockwise_leading_with_hydrogen() {
        let atoms = read(&"[C@@H](N)(O)C").unwrap();
        let parity = atom_parity(0, &atoms[0]);

        assert_eq!(parity, Some(Parity::Positive));
    }

    #[test]
    fn counterclockwise_trailing_with_hydrogen() {
        let atoms = read(&"N[C@H](O)C").unwrap();
        let parity = atom_parity(1, &atoms[1]);

        assert_eq!(parity, Some(Parity::Positive));
    }

    #[test]
    fn clockwise_trailing_with_hydrogen() {
        let atoms = read(&"N[C@@H](O)C").unwrap();
        let parity = atom_parity(1, &atoms[1]);

        assert_eq!(parity, Some(Parity::Negative));
    }

    #[test]
    fn clockwise_trailing_without_hydrogen() {
        let atoms = read(&"N[C@@]([H])(O)C").unwrap();
        let parity = atom_parity(1, &atoms[1]);

        assert_eq!(parity, Some(Parity::Positive));
    }

    #[test]
    fn counterclockwise_leading_with_hydrogen_and_cut() {
        let atoms = read(&"[C@H](C)1OC1").unwrap();
        let parity = atom_parity(0, &atoms[0]);

        assert_eq!(parity, Some(Parity::Negative));
    }

    #[test]
    fn clockwise_leading_with_hydrogen_and_cut() {
        let atoms =read(&"[C@@H](C)1OC1").unwrap();
        let parity = atom_parity(0, &atoms[0]);

        assert_eq!(parity, Some(Parity::Positive));
    }

    #[test]
    fn counterclockwise_trailing_without_hydrogen_with_cut() {
        let atoms = read(&"[H][C@]1(N)OC1").unwrap();
        let parity = atom_parity(1, &atoms[1]);

        assert_eq!(parity, Some(Parity::Negative));
    }

    #[test]
    fn counterclockwise_trailing_without_hydrogen_with_distal_cut() {
        let atoms = read(&"[H][C@](N)1OC1").unwrap();
        let parity = atom_parity(1, &atoms[1]);

        assert_eq!(parity, Some(Parity::Negative));
    }
}