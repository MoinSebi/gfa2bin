use std::collections::{HashMap};
use std::hash::Hash;


#[allow(dead_code)]
pub fn counting<T: PartialEq + Eq + Hash>(vecci: & Vec<T>, what: & mut HashMap<T, usize>) {

    for x in vecci.iter(){
        *what.get_mut(x).unwrap() += 1;
    }
}