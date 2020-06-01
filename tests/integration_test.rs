use gamma::traversal::{ depth_first, breadth_first };
use chemcore::daylight;

#[test]
fn depth_first_propane_from_smiles() {
    let molecule = daylight::Molecule::build(&"CCC").unwrap();
    let traversal = depth_first(&molecule, &0).unwrap();

    assert_eq!(traversal.collect::<Vec<_>>(), vec![
        (&0, &1, false),
        (&1, &2, false)
    ]);
}

#[test]
fn depth_first_propane_from_smiles_inside() {
    let molecule = daylight::Molecule::build(&"C(C)C").unwrap();
    let traversal = depth_first(&molecule, &0).unwrap();

    assert_eq!(traversal.collect::<Vec<_>>(), vec![
        (&0, &1, false),
        (&0, &2, false)
    ]);
}

#[test]
fn depth_first_cyclopropane() {
    let molecule = daylight::Molecule::build(&"C1CC1").unwrap();
    let traversal = depth_first(&molecule, &0).unwrap();

    assert_eq!(traversal.collect::<Vec<_>>(), vec![
        (&0, &2, false),
        (&2, &1, false),
        (&1, &0, true)
    ]);
}

#[test]
fn breadth_first_cyclopropane() {
    let molecule = daylight::Molecule::build(&"C1CC1").unwrap();
    let traversal = breadth_first(&molecule, &0).unwrap();

    assert_eq!(traversal.collect::<Vec<_>>(), vec![
        (&0, &2, false),
        (&0, &1, false),
        (&2, &1, true)
    ]);
}