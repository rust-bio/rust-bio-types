#[cfg_attr(feature = "serde_feature", macro_use)]
#[cfg(feature = "serde_feature")]
extern crate serde_derive;
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
