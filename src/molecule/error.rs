use gamma::graph::Error as GraphError;

#[derive(Debug, PartialEq)]
pub enum Error {
    Valence(usize),
    Isotope(usize),
    Parity(usize),
    Graph(GraphError),
}

impl From<GraphError> for Error {
    fn from(err: GraphError) -> Error {
        Error::Graph(err)
    }
}
