use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug)]
pub struct GraphDB {
    pub nodes: HashMap<String, Node>,
    pub edges: HashSet<Edge>,
}

#[derive(Clone, Debug)]
pub struct Node {
    pub id: String,
    edges: HashSet<String>,
    deleted: bool,
}

#[derive(Clone, Debug)]
pub struct Edge {
    id: String,
    pub attributes: HashMap<String, String>,
    deleted: bool,
}

pub fn get_edge_id(id_a: &String, id_b: &String) -> String {
    format!("{}:{}", id_a, id_b)
}

impl GraphDB {
    pub fn new() -> Self {
        GraphDB {
            nodes: HashMap::new(),
            edges: HashSet::new(),
        }
    }

    pub fn add_node(&mut self, id: String) {
        self.nodes.entry(id.clone()).or_insert(Node{
            id: id.clone(),
            edges: HashSet::new(),
            deleted: false,
        });
    }

    pub fn add_edge(&mut self, id_a: &String, id_b: &String, attributes: HashMap<String, String>) -> Edge {
        let id = get_edge_id(id_a, id_b);
        self.edges.insert(Edge {
            id: id.clone(),
            attributes,
            deleted: false,
        });

        self.add_node(id_a.clone());
        self.add_node(id_b.clone());
        self.nodes.get_mut(id_a).unwrap().edges.insert(id.clone());
        self.nodes.get_mut(id_b).unwrap().edges.insert(id.clone());

        self.get_edge(&id).unwrap().clone()
    }

    pub fn get_node(&self, id: &String) -> Option<&Node> {
        self.nodes.get(id)
    }

    pub fn get_edge(&self, id: &String) -> Option<&Edge> {
        self.edges.iter().find(|edge| !edge.deleted && edge.id.eq(id))
    }

    pub fn get_edges_for_node(&self, node_id: &String) -> Vec<&Edge> {
        match self.nodes.get(node_id) {
            Some(node) => node
                .edges
                .iter()
                .map(|edge_id| self.get_edge(edge_id))
                .filter(|edge| edge.is_some())
                .map(|edge| edge.unwrap())
                .collect(),
            None => Vec::new(),
        }
    }

    pub fn delete_node(&mut self, node_id: &String) {
        let node = self.nodes.get_mut(node_id).unwrap();
        let edge_ids = node.edges.iter();
        for edge_id in edge_ids {
            match self.edges.iter().find(|edge| edge.id.eq(edge_id)) {
                Some(e) => {
                    self.edges.replace(Edge{
                        id: e.id.clone(),
                        attributes: e.attributes.clone(),
                        deleted: true,
                    });
                },
                None => (),
            }
        }
        node.deleted = true;
    }

    pub fn get_other_end(&self, edge: &Edge, node_id: &String) -> Option<&Node> {
        let edge_id = &edge.id;
        let ids: Vec<&str> = edge_id.split(":").collect();
        let other_end = ids.iter().find(|id| id == &node_id);
        match other_end {
            Some(id) => self.get_node(&id.to_string()),
            None => None,
        }
    }

    pub fn get_edge_nodes(&self, edge: &Edge) -> (&Node, &Node) {
        let edge_id = &edge.id;
        let ids: Vec<&str> = edge_id.split(":").collect();
        let node_a = self.get_node(&ids[0].to_string());
        let node_b = self.get_node(&ids[1].to_string());
        (node_a.unwrap(), node_b.unwrap())
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for Edge {}

impl std::hash::Hash for Edge {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use crate::graphdb::{get_edge_id, GraphDB};

    #[test]
    pub fn test_add_node() {
        let mut graph = GraphDB::new();
        graph.add_node("A".to_string());
        graph.add_node("B".to_string());
        assert!(graph.nodes["A"].id.eq(&"A".to_string()));
        assert!(graph.nodes["B"].id.eq(&"B".to_string()));
    }

    #[test]
    pub fn test_add_edge() {
        let mut graph = GraphDB::new();
        graph.add_node("A".to_string());
        graph.add_node("B".to_string());
        let mut attrs = HashMap::new();
        attrs.insert("dg".to_string(), "-12.00".to_string());
        graph.add_edge(&"A".to_string(), &"B".to_string(), attrs);

        let a = graph.get_node(&"A".to_string()).unwrap();
        let b = graph.get_node(&"B".to_string()).unwrap();
        assert_eq!(graph.edges.len(), 1);
        assert_eq!(a.edges.len(), 1);
        assert_eq!(b.edges.len(), 1);
    }

    #[test]
    pub fn test_get_node() {
        let mut graph = GraphDB::new();
        graph.add_node("A".to_string());

        let a = graph.get_node(&"A".to_string());
        assert!(a.is_some());

        let b = graph.get_node(&"B".to_string());
        assert!(b.is_none());
    }

    #[test]
    pub fn test_get_edge() {
        let mut graph = GraphDB::new();
        graph.add_node("A".to_string());
        graph.add_node("B".to_string());
        graph.add_edge(&"A".to_string(), &"B".to_string(), HashMap::new());

        let a = graph.get_edge(&"A:B".to_string());
        assert!(a.is_some());
    }

    #[test]
    pub fn test_get_edges_for_node() {
        let mut graph = GraphDB::new();
        graph.add_node("A".to_string());
        graph.add_node("B".to_string());
        graph.add_node("C".to_string());
        graph.add_edge(&"A".to_string(), &"B".to_string(), HashMap::new());
        graph.add_edge(&"A".to_string(), &"C".to_string(), HashMap::new());

        let a_edges = graph.get_edges_for_node(&"A".to_string());
        assert_eq!(a_edges.len(), 2);

        let b_edges = graph.get_edges_for_node(&"B".to_string());
        assert_eq!(b_edges.len(), 1);
    }

    #[test]
    pub fn test_delete_node() {
        let mut graph = GraphDB::new();
        let id_a = "A".to_string();
        let id_b = "B".to_string();
        let id_c = "C".to_string();
        graph.add_node(id_a.clone());
        graph.add_node(id_b.clone());
        graph.add_node(id_c.clone());
        graph.add_edge(&"A".to_string(), &"B".to_string(), HashMap::new());
        graph.add_edge(&"A".to_string(), &"C".to_string(), HashMap::new());
        graph.delete_node(&id_a);

        let node_a = graph.get_node(&id_a).unwrap();
        assert!(node_a.deleted);
        assert_eq!(graph.get_edges_for_node(&id_a).len(), 0);

        let edge_ab_id = get_edge_id(&id_a, &id_b);
        let edge_a = graph.edges
            .iter()
            .find(|e| e.id.eq(&edge_ab_id))
            .unwrap();
        assert_eq!(edge_a.deleted, true);
    }
}