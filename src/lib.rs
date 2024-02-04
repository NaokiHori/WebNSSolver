#![deny(missing_docs)]

//! A Navier-Stokes solver which runs in web browsers.

mod array;
mod colormap;
mod dct;
mod drawer;
mod simulator;
mod tdm;

use wasm_bindgen::prelude::*;

/// Stores all states.
#[wasm_bindgen]
pub struct WebNSSolver {
    /// Stores and updates a flow field.
    field: crate::simulator::Field,
    /// Draws flow field to HTML canvas element.
    drawer: crate::drawer::Drawer,
}

#[wasm_bindgen]
impl WebNSSolver {
    /// Constructor of [`WebNSSolver`] struct.
    ///
    /// # Arguments
    /// * `ra` : Rayleigh number.
    /// * `pr` : Prandtl number.
    ///
    /// # Returns
    /// A new initialised array.
    pub fn new(ra: f64, pr: f64) -> WebNSSolver {
        let field: simulator::Field = simulator::Field::new(ra, pr);
        let drawer = crate::drawer::Drawer::new(simulator::NPOINTS);
        return WebNSSolver { field, drawer };
    }

    /// Integrates the NS equation in time for a desired time and display the resulting flow field.
    pub fn update(&mut self) -> () {
        // loop until desired time is reached
        const RATE: f64 = 5e-2;
        let mut time: f64 = 0.;
        loop {
            let dt: f64 = simulator::evolve(RATE, &mut self.field);
            time += dt;
            if RATE < time {
                break;
            }
        }
        // draw
        self.drawer.draw(&self.field);
    }

    /// Switches the type of a scalar field to be drawn.
    pub fn change_field(&mut self) -> () {
        self.drawer.change_field();
    }
}

#[wasm_bindgen(start)]
/// Handles things when the WASM is loaded.
///
/// For now nothing to do.
pub fn init() -> () {}
