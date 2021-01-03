#[cfg(feature = "serde")]
extern crate serde;
#[macro_use]
extern crate quick_error;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate derive_new;

extern crate regex;

pub mod alignment;
pub mod annot;
pub mod genome;
pub mod sequence;
pub mod strand;
pub mod variant;
