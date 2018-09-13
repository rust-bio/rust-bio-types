#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate quick_error;
#[macro_use]
extern crate lazy_static;

extern crate regex;

pub mod alignment;
pub mod annot;
pub mod genome;
pub mod sequence;
pub mod strand;
pub mod variant;
