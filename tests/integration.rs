use chemcore::daylight;

#[test]
fn test() {
    let _ = daylight::read_smiles("C", None);
}