// Copyright 2020 Franklin Delehelle
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A phylogenetic tree is represented as a directed graph.
//! Each node is a taxon, identified as a string.
//! The edges are weighted by the phylogenetic distance if it was defined, or f32::NAN otherwise.

use std::{ffi::OsStr, fs, io, path::Path, str::FromStr};

use nom::bytes::complete::take_till;
use petgraph::{
    graph::{Graph, NodeIndex},
    visit::EdgeRef,
    Directed,
    EdgeDirection::Outgoing,
};

pub type Taxon = String;
pub type Proximity = f32;

pub type TreeGraph = Graph<Taxon, Proximity, Directed>;

/// Representation of a phylogenetic tree.
///
/// The root is at NodeIndex 0.
///
/// String conversions and file IO are to/from the [Newick format](https://en.wikipedia.org/wiki/Newick_format).
/// Extra whitespace, quoting strings, and comments are currently not supported.
pub struct Tree {
    pub g: TreeGraph,
}

impl Tree {
    /// Create a new empty Tree.
    pub fn new() -> Self {
        Tree {
            g: TreeGraph::new(),
        }
    }
}

impl ToString for Tree {
    /// Convert the Tree to the Newick format.
    fn to_string(&self) -> String {
        fn subtree_to_string(i: NodeIndex, g: &TreeGraph, mut s: String) -> String {
            let mut iter = g.edges_directed(i, Outgoing).peekable();
            if iter.peek().is_some() {
                s += "(";
                let mut first = true;
                for edge in iter {
                    println!("{:?}", edge.target());
                    if first {
                        first = false;
                    } else {
                        s += ",";
                    }
                    s = subtree_to_string(edge.target(), g, s);
                    if !edge.weight().is_nan() {
                        s += ":";
                        s += &edge.weight().to_string();
                    }
                }
                s += ")";
            }
            if let Some(name) = g.node_weight(i) {
                s += name;
            }
            s
        }
        subtree_to_string(0.into(), &self.g, String::new()) + ";"
    }
}

impl FromStr for Tree {
    type Err = String;

    /// Parse a string in Newick format.
    fn from_str(s: &str) -> Result<Tree, String> {
        use nom::{
            branch::alt,
            bytes::complete::tag,
            combinator::{map, opt},
            multi::separated_list1,
            number::complete::float,
            sequence::{delimited, pair, preceded, terminated, tuple},
            IResult,
        };

        type Result<'a, O> = IResult<&'a str, O>;
        enum ParseTree<'a> {
            Leaf(&'a str),
            Internal((Vec<(ParseTree<'a>, f32)>, &'a str)),
        }

        impl ParseTree<'_> {
            fn to_tree(&self, t: &mut Tree) -> NodeIndex {
                match self {
                    ParseTree::Leaf(name) => t.g.add_node(name.to_string()),
                    ParseTree::Internal((children, name)) => {
                        let node = t.g.add_node(name.to_string());
                        // Add the children in reverse order, so that PetGraph iterates them in the normal order.
                        // Useful for tests.
                        children.iter().rev().for_each(|(pt, d)| {
                            let child_node = pt.to_tree(t);
                            t.g.add_edge(node, child_node, *d);
                        });
                        node
                    }
                }
            }
        }

        // Grammar taken from https://en.wikipedia.org/wiki/Newick_format#Grammar.
        fn length(s: &str) -> Result<f32> {
            map(opt(preceded(tag(":"), float)), |o| o.unwrap_or(f32::NAN))(s)
        }
        fn name(s: &str) -> Result<&str> {
            take_till(|c| ";(),:".find(c).is_some())(s)
        }
        fn leaf(s: &str) -> Result<ParseTree> {
            map(name, ParseTree::Leaf)(s)
        }
        fn branch(s: &str) -> Result<(ParseTree, f32)> {
            tuple((subtree, length))(s)
        }
        fn branchset(s: &str) -> Result<Vec<(ParseTree, f32)>> {
            separated_list1(tag(","), branch)(s)
        }
        fn internal(s: &str) -> Result<ParseTree> {
            map(
                pair(delimited(tag("("), branchset, tag(")")), name),
                ParseTree::Internal,
            )(s)
        }
        fn subtree(s: &str) -> Result<ParseTree> {
            alt((internal, leaf))(s)
        }
        fn tree(s: &str) -> Result<ParseTree> {
            terminated(subtree, tag(";"))(s)
        }
        let mut t = Tree::new();
        map(tree, |pt| pt.to_tree(&mut t))(s).map_err(|x| x.to_string())?;
        Ok(t)
    }
}

impl Tree {
    /// Read from a `.tree` file in Newick format.
    pub fn from_file(p: &Path) -> Result<Self, String> {
        assert!(p.extension() == Some(OsStr::new("tree")));
        fs::read_to_string(p).map_err(|e| e.to_string())?.parse()
    }

    /// Write to a `.tree` file in Newick format.
    pub fn to_file(self: &Self, p: &Path) -> io::Result<()> {
        assert!(p.extension() == Some(OsStr::new("tree")));
        Ok(fs::write(p, self.to_string())?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phylogeny_from_to_string() {
        // From Wikipedia: https://en.wikipedia.org/wiki/Newick_format
        let strings = vec![
            // This is not supported currently
            //"(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;", //all have a distance to parent, including the root.
            ";",
            "A;",
            "(A,B);",
            "(,,(,));",                            //no nodes are named
            "(A,B,(C,D));",                        //leaf nodes are named
            "(A,B,(C,D)E)F;",                      //all nodes are named
            "(:0.1,:0.2,(:0.3,:0.4):0.5);",        //all but root node have a distance to parent
            "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",    //distances and leaf names (popular)
            "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",  //distances and all names
            "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;", //a tree rooted on a leaf node (rare)
        ];
        for s in strings {
            let t = s.parse::<Tree>().unwrap();
            assert_eq!(t.to_string(), s);
        }
    }
}
